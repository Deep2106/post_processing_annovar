#!/usr/bin/env python3
"""
Post-annotation Excel generator - batch PED-driven processing.

Reads ALL families from a PED file and processes each one.
For each family, detects proband (A), father (B), mother (C),
finds the corresponding *hg38_multianno.csv files, and writes
one Excel workbook per proband.

Required arguments:
  -f / --ped          PED file (tab-separated, no header)
  -d / --directory    Directory containing *hg38_multianno.csv files
  -o / --output_dir   Output directory for Excel files

Optional arguments:
  -m / --omim         OMIM summary CSV
  -p / --hpo          HPO genes-to-phenotype TSV
  -r / --reference    Indexed reference FASTA (hg38) for exact indel linkouts

PED format (tab-separated, no header):
  FamilyID  SampleID  FatherID  MotherID  Sex  Phenotype

Role assignment (two supported conventions, auto-detected per sample):
  Legacy   : SampleID ending in A  ->  proband
             SampleID ending in B  ->  father
             SampleID ending in C  ->  mother
  Explicit : SampleID ending in P  ->  proband
             SampleID ending in F  ->  father  (also accepts D as alias)
             SampleID ending in D  ->  father  (alias for F)
             SampleID ending in M  ->  mother
  Both conventions may coexist across families in the same PED file.

Usage:
  python processing_annotated_v5.py \\
      -f cohort.ped \\
      -d /data/annovar \\
      -o /data/output \\
      -m /dbs/OMIM_Summary_File.csv \\
      -p /dbs/hpo_genes_to_phenotype.txt

Changes v8 -> v9 (this file):
  - Quad / multi-child family support added in parse_ped().
    Families with >1 proband sharing the same FamilyID are automatically
    split into one trio per proband, each inheriting the shared parents.
    Single-proband families are unaffected (fully backward-compatible).

Changes v7 -> v8 (this file):
  - DeNovo_Het sheet added to Excel output for duo/trio families only.
    Filter: Inheritance == 'De Novo', GT == 'het', BA1 == False, intronic == False.
    Singletons are unaffected (Inheritance is always 'Singleton' for them).

Changes v6 -> v7 (this file):
  - GEL panel annotation added (-g / --gel argument).
    Adds a GEL_annotation column by matching Gene.refGene against Gene.Symbol
    in the supplied CSV. Unmatched genes are labelled "Not_Present".
    Column appears in output after HPO_phenotypes. Optional — script runs
    without it (all rows labelled Not_Present).

Changes v5 -> v6:
  - Role detection extended: P/F/M/D suffixes now supported in addition to
    the original A/B/C convention. D is treated as an alias for F (father).
    Both conventions can coexist across families in the same PED file.

Changes v4 -> v5:
  - Fixed gnomAD linkout for indels: ANNOVAR strips anchor base (Ref="-" or
    Alt="-") so direct variant URLs break. Now uses:
      SNPs/MNVs  → direct variant link  (ANNOVAR coords = VCF coords)
      Indels     → /variant/{rsID} disambiguation page (if avsnp151 available)
                    Shows all alleles at the site so clinician can select the
                    correct one matching the patient's Alt column
      Indels     → region link ±25 bp   (fallback when no rsID available)
"""

import pandas as pd
import numpy as np
import re
import sys
from glob import glob
from pathlib import Path
from datetime import datetime
import argparse

# Optional: pysam for reference genome lookups (anchor base reconstruction)
try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False

# Module-level reference FASTA handle (opened once in main, shared by all linkouts)
_FASTA = None


def _get_anchor_base(chrom: str, pos_1based: int) -> str:
    """
    Fetch a single reference base from the indexed FASTA.
    Returns empty string if FASTA is not loaded or lookup fails.

    Args:
        chrom:      Chromosome (with or without 'chr' prefix — auto-adjusted)
        pos_1based: 1-based genomic position of the base to fetch
    """
    if _FASTA is None:
        return ""
    try:
        # Ensure chr prefix matches the FASTA contig naming
        chrom_str = str(chrom)
        if not chrom_str.lower().startswith("chr"):
            chrom_str = "chr" + chrom_str

        # pysam uses 0-based half-open coordinates
        base = _FASTA.fetch(chrom_str, pos_1based - 1, pos_1based).upper()
        return base if len(base) == 1 else ""
    except Exception:
        return ""


def _annovar_to_vcf(chrom: str, start: int, end: int, ref: str, alt: str):
    """
    Reconstruct VCF-format alleles from ANNOVAR indel coordinates.
    Returns (vcf_pos, vcf_ref, vcf_alt) or None if anchor base unavailable.

    ANNOVAR insertion (Ref="-"):
      ANNOVAR: Start=227387887, Ref=-, Alt=CCACCGCAG
      Anchor:  ref base at Start = "A"
      VCF:     pos=227387887, Ref=A, Alt=ACCACCGCAG

    ANNOVAR deletion (Alt="-"):
      ANNOVAR: Start=101, End=103, Ref=GCT, Alt=-
      Anchor:  ref base at (Start - 1) = "A"
      VCF:     pos=100, Ref=AGCT, Alt=A
    """
    if ref == "-":
        # Insertion: anchor is the ref base at ANNOVAR Start position
        anchor = _get_anchor_base(chrom, start)
        if not anchor:
            return None
        return (start, anchor, anchor + alt)

    if alt == "-":
        # Deletion: anchor is the ref base BEFORE the deleted region
        anchor_pos = start - 1
        anchor = _get_anchor_base(chrom, anchor_pos)
        if not anchor:
            return None
        return (anchor_pos, anchor + ref, anchor)

    return None


# ─────────────────────────────────────────────────────────────────────────────
# TeeLogger: mirrors all stdout AND stderr to a log file in real time.
# Activated ONLY after the user confirms at the interactive prompt,
# so the confirmation question itself is never affected.
# ─────────────────────────────────────────────────────────────────────────────
class TeeLogger:
    """
    Replaces sys.stdout and sys.stderr with objects that write to both
    the original stream and a log file simultaneously (line-buffered).
    Call .stop() to restore original streams and flush/close the file.
    """
    def __init__(self, log_path: str):
        self.log_path     = log_path
        self._log_file    = open(log_path, "w", buffering=1)   # line-buffered
        self._orig_stdout = sys.stdout
        self._orig_stderr = sys.stderr
        sys.stdout = self._TeeStream(self._orig_stdout, self._log_file)
        sys.stderr = self._TeeStream(self._orig_stderr, self._log_file)

    class _TeeStream:
        def __init__(self, stream, log_file):
            self._stream   = stream
            self._log_file = log_file

        def write(self, data):
            self._stream.write(data)
            self._log_file.write(data)

        def flush(self):
            self._stream.flush()
            self._log_file.flush()

        def __getattr__(self, attr):
            return getattr(self._stream, attr)

    def stop(self):
        sys.stdout = self._orig_stdout
        sys.stderr = self._orig_stderr
        self._log_file.flush()
        self._log_file.close()
        # Print to restored original stdout so the path itself is also visible
        print(f"\n  Log saved -> {self.log_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Desired output column order
# Columns not in this list are appended at the end unchanged.
# Script-generated columns sit between genomicSuperDups and GNOMAD_linkout.
# ─────────────────────────────────────────────────────────────────────────────
DESIRED_COLUMN_ORDER = [
    # Variant coordinates & call info
    "Chr", "Start", "End", "Ref", "Alt", "GT", "DeNovo", "Caller",
    # Gene annotation
    "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
    "ExonicFunc.refGene", "AAChange.refGene",
    "cytoBand", "genomicSuperDups",
    # Script-generated analysis columns
    "OMIM", "HPO_phenotypes", "GEL_annotation", "Inheritance", "LOF", "intronic", "BA1", "variant",
    # External linkouts
    "GNOMAD_linkout", "UCSC_linkout", "Decipher_linkout",
    # Top-tier pathogenicity scores
    "REVEL_score",
    "ClinPred_score", "ClinPred_rankscore", "ClinPred_pred",
    "SpliceAI_gene", "SpliceAI_DS_AG", "SpliceAI_DS_AL",
    "SpliceAI_DS_DG", "SpliceAI_DS_DL", "SpliceAI_max", "SpliceAI_raw",
    "AlphaMissense_score", "AlphaMissense_rankscore", "AlphaMissense_pred",
    # ClinVar
    "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG",
    # Population frequencies - ExAC
    "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS",
    "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS",
    # Population frequencies - gnomAD 4.1 exome
    "gnomad41_exome_AF", "gnomad41_exome_AF_raw",
    "gnomad41_exome_AF_XX", "gnomad41_exome_AF_XY",
    "gnomad41_exome_AF_grpmax", "gnomad41_exome_faf95", "gnomad41_exome_faf99",
    "gnomad41_exome_fafmax_faf95_max", "gnomad41_exome_fafmax_faf99_max",
    "gnomad41_exome_AF_afr", "gnomad41_exome_AF_amr", "gnomad41_exome_AF_asj",
    "gnomad41_exome_AF_eas", "gnomad41_exome_AF_fin", "gnomad41_exome_AF_mid",
    "gnomad41_exome_AF_nfe", "gnomad41_exome_AF_remaining", "gnomad41_exome_AF_sas",
    # Population frequencies - gnomAD 4.1 genome
    "gnomad41_genome_AF", "gnomad41_genome_AF_raw",
    "gnomad41_genome_AF_XX", "gnomad41_genome_AF_XY",
    "gnomad41_genome_AF_grpmax", "gnomad41_genome_faf95", "gnomad41_genome_faf99",
    "gnomad41_genome_fafmax_faf95_max", "gnomad41_genome_fafmax_faf99_max",
    "gnomad41_genome_AF_afr", "gnomad41_genome_AF_ami",
    "gnomad41_genome_AF_amr", "gnomad41_genome_AF_asj", "gnomad41_genome_AF_eas",
    "gnomad41_genome_AF_fin", "gnomad41_genome_AF_mid", "gnomad41_genome_AF_nfe",
    "gnomad41_genome_AF_remaining", "gnomad41_genome_AF_sas",
    # dbSNP
    "avsnp151",
    # Functional predictors
    "SIFT_score", "SIFT_converted_rankscore", "SIFT_pred",
    "SIFT4G_score", "SIFT4G_converted_rankscore", "SIFT4G_pred",
    "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore", "Polyphen2_HDIV_pred",
    "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore", "Polyphen2_HVAR_pred",
    "LRT_score", "LRT_converted_rankscore", "LRT_pred", "LRT_Omega",
    "MutationTaster_score", "MutationTaster_converted_rankscore", "MutationTaster_pred",
    "MutationAssessor_score", "MutationAssessor_rankscore", "MutationAssessor_pred",
    "FATHMM_score", "FATHMM_converted_rankscore", "FATHMM_pred",
    "PROVEAN_score", "PROVEAN_converted_rankscore", "PROVEAN_pred",
    "VEST4_score", "VEST4_rankscore",
    "MetaSVM_score", "MetaSVM_rankscore", "MetaSVM_pred",
    "MetaLR_score", "MetaLR_rankscore", "MetaLR_pred",
    "Reliability_index",
    "MetaRNN_score", "MetaRNN_rankscore", "MetaRNN_pred",
    "M-CAP_score", "M-CAP_rankscore", "M-CAP_pred",
    "REVEL_rankscore",
    "MutPred_score", "MutPred_rankscore",
    "MVP_score", "MVP_rankscore",
    "gMVP_score", "gMVP_rankscore",
    "MPC_score", "MPC_rankscore",
    "PrimateAI_score", "PrimateAI_rankscore", "PrimateAI_pred",
    "DEOGEN2_score", "DEOGEN2_rankscore", "DEOGEN2_pred",
    "BayesDel_addAF_score", "BayesDel_addAF_rankscore", "BayesDel_addAF_pred",
    "BayesDel_noAF_score", "BayesDel_noAF_rankscore", "BayesDel_noAF_pred",
    "LIST-S2_score", "LIST-S2_rankscore", "LIST-S2_pred",
    "VARITY_R_score", "VARITY_R_rankscore",
    "VARITY_ER_score", "VARITY_ER_rankscore",
    "VARITY_R_LOO_score", "VARITY_R_LOO_rankscore",
    "VARITY_ER_LOO_score", "VARITY_ER_LOO_rankscore",
    "ESM1b_score", "ESM1b_rankscore", "ESM1b_pred",
    "EVE_score", "EVE_rankscore",
    "Aloft_pred", "Aloft_Confidence",
    "CADD_raw", "CADD_raw_rankscore", "CADD_phred",
    "DANN_score", "DANN_rankscore",
    "fathmm-MKL_coding_score", "fathmm-MKL_coding_rankscore",
    "fathmm-MKL_coding_pred", "fathmm-MKL_coding_group",
    "fathmm-XF_coding_score", "fathmm-XF_coding_rankscore", "fathmm-XF_coding_pred",
    "Eigen-raw_coding", "Eigen-raw_coding_rankscore", "Eigen-phred_coding",
    "Eigen-PC-raw_coding", "Eigen-PC-raw_coding_rankscore", "Eigen-PC-phred_coding",
    "GenoCanyon_score", "GenoCanyon_rankscore",
    "integrated_fitCons_score", "integrated_fitCons_rankscore", "integrated_confidence_value",
    "GM12878_fitCons_score", "GM12878_fitCons_rankscore", "GM12878_confidence_value",
    "H1-hESC_fitCons_score", "H1-hESC_fitCons_rankscore", "H1-hESC_confidence_value",
    "HUVEC_fitCons_score", "HUVEC_fitCons_rankscore", "HUVEC_confidence_value",
    "LINSIGHT", "LINSIGHT_rankscore",
    "GERP++_NR", "GERP++_RS", "GERP++_RS_rankscore",
    "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore",
    "phyloP470way_mammalian", "phyloP470way_mammalian_rankscore",
    "phyloP17way_primate", "phyloP17way_primate_rankscore",
    "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore",
    "phastCons470way_mammalian", "phastCons470way_mammalian_rankscore",
    "phastCons17way_primate", "phastCons17way_primate_rankscore",
    "SiPhy_29way_pi", "SiPhy_29way_logOdds", "SiPhy_29way_logOdds_rankscore",
    "bStatistic", "bStatistic_converted_rankscore",
    "Interpro_domain",
    "GTEx_V8_eQTL_gene", "GTEx_V8_eQTL_tissue",
    "GTEx_V8_sQTL_gene", "GTEx_V8_sQTL_tissue",
    "eQTLGen_snp_id",
    "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE",
    "MCAP",
    # ONC / SCI
    "ONCDN", "ONCDISDB", "ONCREVSTAT", "ONC",
    "SCIDN", "SCIDISDB", "SCIREVSTAT", "SCI",
    # REVEL standalone + study metadata
    "REVEL",
    "DN ID", "Patient ID", "Phenotype", "Platform", "Study", "Pubmed ID",
]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_base_id(sample_id: str) -> str:
    """
    Split SampleID on the last underscore and return (family_part, role_letter).

    The PED SampleID format is always:  {family_id}_{ROLE}
    where ROLE is a single letter as the last underscore-delimited token.

    Supported role letters:
      P  - proband       F / D  - father       M  - mother
      A  - proband (legacy)     B - father (legacy)   C - mother (legacy)

    Examples:
      CINDI013_EPI013SB_P  ->  family=CINDI013_EPI013SB,  role=P
      CINDI013_EPI013SB_F  ->  family=CINDI013_EPI013SB,  role=F
      CINDI013_EPI013SB_M  ->  family=CINDI013_EPI013SB,  role=M
      LH0001A              ->  family=LH0001,             role=A  (legacy, no underscore)
    """
    s = str(sample_id)
    if "_" in s:
        family_part, role_letter = s.rsplit("_", 1)
        return family_part, role_letter.upper()
    # Legacy format with no trailing underscore: role letter is the last character
    return s[:-1], s[-1].upper()


def family_type_label(n: int) -> str:
    return {1: "singleton", 2: "duo", 3: "trio"}.get(n, "unknown")


def find_multianno_files(annovar_dir: str, base_id: str) -> list:
    patterns = [
        f"{annovar_dir}/*{base_id}*.hg38_multianno.csv",
        f"{annovar_dir}/{base_id}*.hg38_multianno.csv",
        f"{annovar_dir}/{base_id}.hg38_multianno.csv",
    ]
    found = []
    for pat in patterns:
        found.extend(glob(pat))
    return sorted(set(found))


# ─────────────────────────────────────────────────────────────────────────────
# PED parsing  ->  returns dict of { family_id: {role: base_id, ...}, ... }
# ─────────────────────────────────────────────────────────────────────────────
def parse_ped(ped_path: str) -> dict:
    """
    Parse entire PED file and return all families.

    Returns:
        {
          "LH0001": {"proband": "LH0001A", "father": "LH0001B", "mother": "LH0001C"},  # A/B/C convention
          "FAM02":  {"proband": "FAM02P",  "father": "FAM02F",  "mother": "FAM02M"},   # P/F/M convention
          "FAM03":  {"proband": "FAM03P",  "father": "FAM03D",  "mother": "FAM03M"},   # D as father alias
          "LH0002": {"proband": "LH0002A"},
          ...
        }
    """
    try:
        ped = pd.read_csv(
            ped_path, sep="\t", header=None,
            names=["FamilyID", "SampleID", "FatherID", "MotherID", "Sex", "Phenotype"],
            comment="#",
            dtype=str,
        )
    except Exception as e:
        print(f"ERROR: Could not read PED file '{ped_path}': {e}", file=sys.stderr)
        sys.exit(1)

    # Drop founder/missing rows (SampleID == 0 or empty)
    ped = ped[~ped["SampleID"].isin(["0", "nan", ""])].copy()

    # Legacy convention:   A=proband, B=father, C=mother
    # Explicit convention: P=proband, F=father, M=mother
    # D is an alias for father (equivalent to B/F)
    role_map = {
        "A": "proband", "B": "father", "C": "mother",   # legacy A/B/C
        "P": "proband", "F": "father", "M": "mother",   # explicit P/F/M
        "D": "father",                                   # D alias for father
    }

    # Use _probands list to accumulate all probands per family before exploding.
    # This handles quad/multi-child families where >1 proband shares a family ID.
    raw_families = {}
    for _, row in ped.iterrows():
        fid  = str(row["FamilyID"]).strip()
        sid  = str(row["SampleID"]).strip()
        family_part, letter = get_base_id(sid)
        role = role_map.get(letter)

        if role is None:
            print(f"  Warning: SampleID '{sid}' (family {fid}) — last token "
                  f"'{letter}' is not a recognised role letter (P/F/M/D or A/B/C) — skipped.",
                  file=sys.stderr)
            continue

        if role == "proband":
            # Collect into a list — there may be >1 proband (quad families)
            raw_families.setdefault(fid, {}).setdefault("_probands", []).append(sid)
        else:
            raw_families.setdefault(fid, {})[role] = sid

    # Explode multi-proband families into one entry per proband, each sharing
    # the same parents. Single-proband families keep their family ID as key.
    families = {}
    for fid, members in raw_families.items():
        probands = members.pop("_probands", [])
        parents  = dict(members)  # father / mother (if present)

        if not probands:
            print(f"  Warning: Family '{fid}' has no proband — skipped.",
                  file=sys.stderr)
            continue

        if len(probands) == 1:
            # Standard singleton/duo/trio — backward-compatible
            families[fid] = {"proband": probands[0], **parents}
        else:
            # Quad / multi-child: split into one trio per proband
            print(f"  [INFO] Family '{fid}': {len(probands)} probands detected "
                  f"— splitting into {len(probands)} separate trios.")
            for p_sid in probands:
                families[p_sid] = {"proband": p_sid, **parents}

    if not families:
        print("ERROR: No valid families found in PED file.", file=sys.stderr)
        sys.exit(1)

    return families


# ─────────────────────────────────────────────────────────────────────────────
# Reference databases
# ─────────────────────────────────────────────────────────────────────────────
def get_omim(omim_path: str) -> pd.DataFrame:
    try:
        omim = pd.read_csv(omim_path)
        omim["Approved Symbol"] = omim["Approved Symbol"].astype(str).str.upper()
        return omim
    except Exception as e:
        print(f"Warning: Could not load OMIM file: {e}", file=sys.stderr)
        return pd.DataFrame()


def get_hpo(hpo_path: str) -> dict:
    try:
        hpo = pd.read_csv(hpo_path, sep="\t")
        return (
            hpo[["gene_symbol", "hpo_name"]]
            .groupby("gene_symbol")
            .agg(" | ".join)
            .to_dict(orient="index")
        )
    except Exception as e:
        print(f"Warning: Could not load HPO file: {e}", file=sys.stderr)
        return {}


def get_gel_panel(gel_path: str) -> dict:
    """
    Load GEL panel annotation file and return a gene → annotation lookup dict.

    Expected columns: Gene.Symbol, GEL_annotation
    Gene symbols are upper-cased for case-insensitive matching.

    Returns:
        {"AARS": "Green_202", "ABAT": "Green_202", ...}
        Empty dict if file is not provided or cannot be loaded.
    """
    if not gel_path:
        return {}
    try:
        gel = pd.read_csv(gel_path)
        # Accept common column name variants
        col_map = {c.strip().lower(): c for c in gel.columns}
        sym_col = col_map.get("gene.symbol") or col_map.get("gene_symbol") or col_map.get("genesymbol")
        ann_col = col_map.get("gel_annotation") or col_map.get("gel annotation")
        if sym_col is None or ann_col is None:
            print(f"  Warning: GEL panel file must have 'Gene.Symbol' and "
                  f"'GEL_annotation' columns — found: {list(gel.columns)}",
                  file=sys.stderr)
            return {}
        gel[sym_col] = gel[sym_col].astype(str).str.upper()
        return dict(zip(gel[sym_col], gel[ann_col]))
    except Exception as e:
        print(f"  Warning: Could not load GEL panel file: {e}", file=sys.stderr)
        return {}


# ─────────────────────────────────────────────────────────────────────────────
# Linkout generators (Excel HYPERLINK formulas)
# ─────────────────────────────────────────────────────────────────────────────
def _is_valid_rsid(val) -> bool:
    """Check if a value is a usable rsID (not missing/empty/dot)."""
    if pd.isna(val):
        return False
    s = str(val).strip()
    return bool(s) and s != "." and s.lower().startswith("rs")


def make_gnomad_link(row: pd.Series) -> str:
    """
    Build gnomAD linkout handling ANNOVAR's indel coordinate conventions.

    Strategy (prioritised):
      1. SNPs/MNVs               → direct variant link (ANNOVAR coords = VCF)
      2. Indels + ref genome      → direct variant link (anchor base reconstructed)
      3. Indels + rsID (no ref)   → /variant/{rsID} disambiguation page
      4. Indels (no ref, no rsID) → region link ±25 bp

    Excel label:
      "gnomAD"          direct variant link
      "gnomAD:rsXXX"    rsID disambiguation page
      "gnomAD:region"   region browser fallback
    """
    try:
        chrom = str(row["Chr"])
        start = int(row["Start"])
        end   = int(row["End"])
        ref   = str(row["Ref"])
        alt   = str(row["Alt"])

        # Strip chr prefix — gnomAD URLs use bare chromosome numbers
        chrom_clean = re.sub(r"(?i)^chr", "", chrom)

        # ── Case 1: SNPs and MNVs ──────────────────────────────────────────
        if ref != "-" and alt != "-":
            url = (
                f"https://gnomad.broadinstitute.org/variant/"
                f"{chrom_clean}-{start}-{ref}-{alt}?dataset=gnomad_r4"
            )
            return f'=HYPERLINK("{url}","gnomAD")'

        # ── Case 2: Indels with reference genome → direct variant link ─────
        vcf = _annovar_to_vcf(chrom, start, end, ref, alt)
        if vcf is not None:
            vcf_pos, vcf_ref, vcf_alt = vcf
            url = (
                f"https://gnomad.broadinstitute.org/variant/"
                f"{chrom_clean}-{vcf_pos}-{vcf_ref}-{vcf_alt}?dataset=gnomad_r4"
            )
            return f'=HYPERLINK("{url}","gnomAD")'

        # ── Case 3: Indels with rsID → disambiguation page ─────────────────
        rsid = row.get("avsnp151", None)
        if _is_valid_rsid(rsid):
            rsid_str = str(rsid).strip()
            label = f"gnomAD:{rsid_str[:15]}" if len(rsid_str) > 15 else f"gnomAD:{rsid_str}"
            url = (
                f"https://gnomad.broadinstitute.org/variant/"
                f"{rsid_str}?dataset=gnomad_r4"
            )
            return f'=HYPERLINK("{url}","{label}")'

        # ── Case 4: Indels without rsID → region fallback ──────────────────
        region_start = max(1, start - 25)
        region_end   = end + 25
        url = (
            f"https://gnomad.broadinstitute.org/region/"
            f"{chrom_clean}-{region_start}-{region_end}?dataset=gnomad_r4"
        )
        return f'=HYPERLINK("{url}","gnomAD:region")'

    except Exception:
        return ""


def make_ucsc_link(row: pd.Series) -> str:
    try:
        chrom = str(row["Chr"])
        if not chrom.lower().startswith("chr"):
            chrom = "chr" + chrom
        start = int(row["Start"])
        end   = int(row["End"])
        url = (f"https://genome.ucsc.edu/cgi-bin/hgTracks"
               f"?db=hg38&position={chrom}:{max(1, start - 50)}-{end + 50}")
        return f'=HYPERLINK("{url}","UCSC")'
    except Exception:
        return ""


def make_decipher_link(row: pd.Series) -> str:
    """
    Build DECIPHER linkout handling ANNOVAR's indel coordinate conventions.

    Strategy (prioritised):
      1. SNPs/MNVs               → chrom:pos REF>ALT (exact variant page)
      2. Indels + ref genome      → chrom:pos REF>ALT (anchor base reconstructed)
      3. Deletions (no ref)       → SPDI chrom:pos:ref: (ANNOVAR preserves deleted bases)
      4. Insertions (no ref)      → position window ±10bp (fallback)

    Excel label:
      "DECIPHER"          exact variant page
      "DECIPHER:region"   position window fallback (insertions without ref genome)
    """
    try:
        chrom = re.sub(r"(?i)^chr", "", str(row["Chr"]))
        start = int(row["Start"])
        end   = int(row["End"])
        ref   = str(row["Ref"])
        alt   = str(row["Alt"])

        # ── Case 1: SNPs and MNVs ──────────────────────────────────────────
        if ref != "-" and alt != "-":
            query = f"{chrom}:{start} {ref}>{alt}"
            url = (
                f"https://www.deciphergenomics.org/search"
                f"?q={query.replace(' ', '%20')}"
            )
            return f'=HYPERLINK("{url}","DECIPHER")'

        # ── Case 2: Indels with reference genome → exact variant page ──────
        vcf = _annovar_to_vcf(row["Chr"], start, end, ref, alt)
        if vcf is not None:
            vcf_pos, vcf_ref, vcf_alt = vcf
            query = f"{chrom}:{vcf_pos} {vcf_ref}>{vcf_alt}"
            url = (
                f"https://www.deciphergenomics.org/search"
                f"?q={query.replace(' ', '%20')}"
            )
            return f'=HYPERLINK("{url}","DECIPHER")'

        # ── Case 3: Deletions without ref → SPDI (deleted bases preserved) ─
        if alt == "-":
            spdi_pos = start - 1
            url = (
                f"https://www.deciphergenomics.org/search"
                f"?q={chrom}:{spdi_pos}:{ref}:"
            )
            return f'=HYPERLINK("{url}","DECIPHER")'

        # ── Case 4: Insertions without ref → position window fallback ──────
        if ref == "-":
            pad_start = max(1, start - 10)
            pad_end   = end + 10
            url = (
                f"https://www.deciphergenomics.org/search"
                f"?q={chrom}:{pad_start}-{pad_end}"
            )
            return f'=HYPERLINK("{url}","DECIPHER:region")'

        return ""

    except Exception:
        return ""


# ─────────────────────────────────────────────────────────────────────────────
# Load and prepare one sample's multianno CSV
# ─────────────────────────────────────────────────────────────────────────────
def prepare_data(annovar_dir: str, base_id: str, omim: pd.DataFrame, gel: dict = None) -> pd.DataFrame:
    files = find_multianno_files(annovar_dir, base_id)
    if not files:
        print(f"  ERROR: No multianno file found for '{base_id}'", file=sys.stderr)
        return pd.DataFrame()

    dfs = []
    for path in files:
        print(f"    -> {Path(path).name}")
        try:
            dfs.append(pd.read_csv(path, sep=",", low_memory=False))
        except Exception as e:
            print(f"    Warning: Failed to load {Path(path).name}: {e}", file=sys.stderr)

    if not dfs:
        return pd.DataFrame()

    df = pd.concat(dfs, ignore_index=True)

    required = ["Chr", "Start", "End", "Ref", "Alt", "GT", "Gene.refGene"]
    missing  = [c for c in required if c not in df.columns]
    if missing:
        print(f"  ERROR: Missing required columns for {base_id}: {', '.join(missing)}",
              file=sys.stderr)
        return pd.DataFrame()

    df["Gene.refGene"] = df["Gene.refGene"].astype(str).str.upper()

    # OMIM lookup
    if not omim.empty:
        omim_lookup = (
            omim[["Approved Symbol", "Phenotypes"]]
            .groupby("Approved Symbol")
            .agg(" | ".join)
            .to_dict(orient="index")
        )
        df["OMIM"] = df["Gene.refGene"].map(omim_lookup).apply(
            lambda x: (
                x["Phenotypes"]
                 .replace("0 |", "").replace("| 0", "").replace("0 | 0", "#N/A").strip()
            ) if isinstance(x, dict) else "#N/A"
        )
    else:
        df["OMIM"] = "#N/A"

    # GEL panel annotation
    if gel:
        df["GEL_annotation"] = df["Gene.refGene"].map(gel).fillna("Not_Present")
    else:
        df["GEL_annotation"] = "Not_Present"

    # ExonicFunc fill
    if "ExonicFunc.refGene" in df.columns:
        df["ExonicFunc.refGene"] = (
            df["ExonicFunc.refGene"].fillna("missing_value").replace(".", "missing_value")
        )
    else:
        df["ExonicFunc.refGene"] = "missing_value"

    # BA1 flag (any gnomAD / ExAC AF >= 0.01)
    af_cols = [c for c in df.columns if "gnomad" in c.lower() or "exac" in c.lower()]
    if af_cols:
        df[af_cols] = (
            df[af_cols].fillna(0).replace(".", 0)
            .apply(pd.to_numeric, errors="coerce").fillna(0)
        )
        df["BA1"] = df[af_cols].ge(0.01).any(axis=1)
    else:
        df["BA1"] = False

    # Variant key
    df["variant"] = (
        df["Chr"].astype(str) + ":" + df["Start"].astype(str) + ":" +
        df["Ref"].astype(str) + ":" + df["Alt"].astype(str)
    )

    return df


# ─────────────────────────────────────────────────────────────────────────────
# Column reordering
# ─────────────────────────────────────────────────────────────────────────────
def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    present  = [c for c in DESIRED_COLUMN_ORDER if c in df.columns]
    leftover = [c for c in df.columns if c not in set(DESIRED_COLUMN_ORDER)]
    if leftover:
        print(f"    Note: {len(leftover)} unrecognised column(s) appended at end.")
    return df[present + leftover]


# ─────────────────────────────────────────────────────────────────────────────
# Analysis flags + linkouts + Excel writing (one family)
# ─────────────────────────────────────────────────────────────────────────────
def analyze_and_write(
    proband_df: pd.DataFrame,
    family_type: str,
    output_dir: str,
    proband_id: str,
) -> None:

    # LOF flag
    lof_consequences = {
        "frameshift deletion", "frameshift insertion",
        "stopgain", "startloss", "stoploss",
    }
    func_col = proband_df.get("Func.refGene", pd.Series(dtype=str))
    proband_df["LOF"] = (
        proband_df["ExonicFunc.refGene"].isin(lof_consequences) |
        func_col.str.contains("splicing", na=False)
    )

    # Intronic flag
    proband_df["intronic"] = (func_col == "intronic")

    # Linkouts — gnomAD uses row-level function for indel-safe handling
    proband_df["GNOMAD_linkout"]   = proband_df.apply(make_gnomad_link, axis=1)
    proband_df["UCSC_linkout"]     = proband_df.apply(make_ucsc_link, axis=1)
    proband_df["Decipher_linkout"] = proband_df.apply(make_decipher_link, axis=1)

    # Reorder columns
    proband_df = reorder_columns(proband_df)

    # Build sheets
    sheets = {}

    sheets["Raw_Data"] = proband_df.copy()

    # DeNovo_Het: only meaningful for duo/trio families where Inheritance is
    # assigned. Singletons always have Inheritance='Singleton' so this sheet
    # is intentionally omitted for them (mirrors post_annotation.py main_trio logic).
    if family_type in ("duo", "trio"):
        sheets["DeNovo_Het"] = proband_df.query(
            "Inheritance == 'De Novo' and GT == 'het' and BA1 == False and intronic == False"
        ).copy()

    sheets["rare_Hom"] = proband_df.query(
        "GT == 'hom' and BA1 == False and intronic == False"
    ).copy()

    het_clean   = proband_df.query("GT == 'het' and BA1 == False and intronic == False")
    gene_counts = het_clean["Gene.refGene"].value_counts()
    comp_genes  = gene_counts[gene_counts > 1].index
    sheets["Comp_Het"] = (
        het_clean[het_clean["Gene.refGene"].isin(comp_genes)]
        .sort_values("Gene.refGene").copy()
    )

    sheets["Het"] = het_clean.copy()

    sheets["Synonymous"] = proband_df.query(
        "`ExonicFunc.refGene` == 'synonymous SNV' and BA1 == False and intronic == False"
    ).copy()

    # Write Excel
    outfile = Path(output_dir) / f"{proband_id}_analysis.xlsx"
    with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
        for name, df_sheet in sheets.items():
            if df_sheet.empty and name != "Raw_Data":
                print(f"    Skipping empty sheet: {name}")
                continue
            df_sheet.to_excel(writer, sheet_name=name, index=False)
            print(f"    Sheet '{name}': {len(df_sheet):,} rows")

    print(f"  -> Written: {outfile}")


# ─────────────────────────────────────────────────────────────────────────────
# Process one family
# ─────────────────────────────────────────────────────────────────────────────
def process_family(
    family_id: str,
    members: dict,
    annovar_dir: str,
    output_dir: str,
    omim: pd.DataFrame,
    hpo: dict,
    gel: dict = None,
) -> bool:
    """Returns True on success, False if proband data is unavailable."""
    ftype      = family_type_label(len(members))
    proband_id = members["proband"]

    print(f"\n{'─' * 70}")
    print(f"  Family : {family_id}  |  Type : {ftype.upper()}  |  Proband : {proband_id}")
    print(f"{'─' * 70}")

    # Load proband
    print(f"  Loading proband ({proband_id}):")
    proband_df = prepare_data(annovar_dir, proband_id, omim, gel)
    if proband_df.empty:
        print(f"  SKIPPING family {family_id} - proband data unavailable.", file=sys.stderr)
        return False

    # HPO annotation
    proband_df["HPO_phenotypes"] = proband_df["Gene.refGene"].map(hpo).apply(
        lambda x: (
            x.get("hpo_name", "#N/A")
             .replace("0 |", "").replace("| 0", "").replace("0 | 0", "#N/A").strip()
        ) if isinstance(x, dict) else "#N/A"
    )

    # Inheritance
    if ftype in ("duo", "trio"):
        father_df = pd.DataFrame()
        mother_df = pd.DataFrame()

        if "father" in members:
            print(f"  Loading father  ({members['father']}):")
            father_df = prepare_data(annovar_dir, members["father"], omim, gel)

        if "mother" in members:
            print(f"  Loading mother  ({members['mother']}):")
            mother_df = prepare_data(annovar_dir, members["mother"], omim, gel)

        in_father = (
            proband_df["variant"].isin(father_df["variant"])
            if not father_df.empty
            else pd.Series(False, index=proband_df.index)
        )
        in_mother = (
            proband_df["variant"].isin(mother_df["variant"])
            if not mother_df.empty
            else pd.Series(False, index=proband_df.index)
        )

        conditions = [in_father & in_mother, in_mother & ~in_father, in_father & ~in_mother]
        choices    = ["Familial", "Maternal", "Paternal"]
        proband_df["Inheritance"] = np.select(conditions, choices, default="De Novo")
    else:
        proband_df["Inheritance"] = "Singleton"

    analyze_and_write(proband_df, ftype, output_dir, proband_id)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Batch confirmation prompt
# ─────────────────────────────────────────────────────────────────────────────
def print_batch_summary_and_ask(
    families: dict,
    annovar_dir: str,
    output_dir: str,
    omim_path: str,
    hpo_path: str,
    ped_path: str,
    ref_path: str = None,
) -> bool:
    print("\n" + "=" * 80)
    print(" " * 22 + "ANNOVAR BATCH POST-PROCESSING SUMMARY")
    print("=" * 80)
    print(f"  PED file       : {ped_path}")
    print(f"  Annovar dir    : {annovar_dir}")
    print(f"  Output dir     : {output_dir}")
    print(f"  OMIM db        : {omim_path}")
    print(f"  HPO db         : {hpo_path}")
    print(f"  Reference      : {ref_path or 'not provided (indel linkouts use fallbacks)'}")
    print(f"  Families found : {len(families)}")
    print()

    print(f"  {'Family':<12} {'Type':<12} {'Proband':<16} {'Father':<16} {'Mother':<16} Files")
    print(f"  {'-'*12} {'-'*12} {'-'*16} {'-'*16} {'-'*16} {'-'*6}")

    all_files_present = True
    for fid, members in sorted(families.items()):
        ftype   = family_type_label(len(members))
        proband = members.get("proband", "-")
        father  = members.get("father",  "-")
        mother  = members.get("mother",  "-")
        n_files = sum(len(find_multianno_files(annovar_dir, mid)) for mid in members.values())
        flag    = "  *** MISSING ***" if n_files < len(members) else ""
        if flag:
            all_files_present = False
        print(f"  {fid:<12} {ftype:<12} {proband:<16} {father:<16} {mother:<16} {n_files}{flag}")

    if not all_files_present:
        print("\n  WARNING: Some multianno files are missing (see *** above).")
        print("           Those samples will be skipped during processing.")

    print()
    print("  Proceed with all families listed above?")
    while True:
        ans = input("  yes/no: ").strip().lower()
        if ans in ("yes", "y"):
            print("\n  Starting batch processing...\n")
            return True
        if ans in ("no", "n"):
            print("  Aborted by user.")
            return False
        print("  Please type yes or no.")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main(ped_path, annovar_dir, output_dir, omim_path, hpo_path, gel_path=None, ref_path=None):
    global _FASTA

    # Validate inputs
    if not Path(ped_path).is_file():
        print(f"ERROR: PED file not found: {ped_path}", file=sys.stderr)
        sys.exit(1)
    if not Path(annovar_dir).is_dir():
        print(f"ERROR: Annovar directory not found: {annovar_dir}", file=sys.stderr)
        sys.exit(1)

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Parse all families from PED
    families = parse_ped(ped_path)
    print(f"\n  Parsed {len(families)} valid family/families from PED.")

    # Confirmation (interactive - before log starts)
    if not print_batch_summary_and_ask(
        families, annovar_dir, output_dir, omim_path, hpo_path, ped_path, ref_path
    ):
        sys.exit(0)

    # Start logging AFTER user confirms
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path  = Path(output_dir) / f"run_{timestamp}.log"
    logger    = TeeLogger(str(log_path))

    try:
        # Write run parameters to log header
        print("=" * 80)
        print(f"  Run started   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  PED file      : {ped_path}")
        print(f"  Annovar dir   : {annovar_dir}")
        print(f"  Output dir    : {output_dir}")
        print(f"  OMIM db       : {omim_path}")
        print(f"  HPO db        : {hpo_path}")
        print(f"  Reference     : {ref_path or 'not provided'}")
        print(f"  GEL panel     : {gel_path or 'not provided (GEL_annotation = Not_Present for all)'}")
        print(f"  Families      : {', '.join(sorted(families.keys()))}")
        print("=" * 80 + "\n")

        # Load reference DBs once - shared across all families
        print("  Loading reference databases...")
        omim = get_omim(omim_path)
        hpo  = get_hpo(hpo_path)
        gel  = get_gel_panel(gel_path)
        print(f"  GEL panel     : {len(gel):,} genes loaded" if gel else "  GEL panel     : not loaded")

        # Open reference FASTA for anchor base reconstruction (indel linkouts)
        if ref_path and Path(ref_path).is_file():
            if not _PYSAM_AVAILABLE:
                print("  WARNING: pysam not installed — indel linkouts will use fallbacks.")
                print("           Install with: pip install pysam")
            else:
                fai_path = ref_path + ".fai"
                if not Path(fai_path).is_file():
                    print(f"  WARNING: FASTA index not found ({fai_path})")
                    print("           Run: samtools faidx {ref_path}")
                    print("           Indel linkouts will use fallbacks.")
                else:
                    _FASTA = pysam.FastaFile(ref_path)
                    print(f"  Reference FASTA loaded: {ref_path}")
                    print(f"    -> Indel linkouts will use exact variant URLs")
        elif ref_path:
            print(f"  WARNING: Reference FASTA not found: {ref_path}")
            print("           Indel linkouts will use fallbacks.")
        else:
            print("  No reference FASTA provided (-r) — indel linkouts use rsID/region fallbacks.")

        print("  Done.\n")

        # Batch loop
        succeeded = []
        failed    = []

        for family_id, members in sorted(families.items()):
            ok = process_family(family_id, members, annovar_dir, output_dir, omim, hpo, gel)
            (succeeded if ok else failed).append(family_id)

        # Final report
        print("\n" + "=" * 80)
        print(" " * 30 + "BATCH COMPLETE")
        print("=" * 80)
        print(f"  Succeeded : {len(succeeded)}  ->  {', '.join(succeeded) or '-'}")
        print(f"  Failed    : {len(failed)}  ->  {', '.join(failed) or '-'}")
        print(f"  Output    : {output_dir}")
        print(f"  Finished  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80 + "\n")

    finally:
        # Close reference FASTA if opened
        if _FASTA is not None:
            _FASTA.close()
            _FASTA = None
        # Always restore streams and close log, even on unexpected error
        logger.stop()


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch post-annotation Excel generator - full PED-driven processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Minimal - process every family in the PED
  python processing_annotated_v5.py \\
      -f cohort.ped \\
      -d /data/annovar \\
      -o /data/output

  # With explicit DB paths
  python processing_annotated_v5.py \\
      -f cohort.ped \\
      -d /data/annovar \\
      -o /data/output \\
      -m /dbs/OMIM_Summary_File.csv \\
      -p /dbs/hpo_genes_to_phenotype.txt
        """
    )
    parser.add_argument(
        "-f", "--ped", type=str, required=True,
        help="Path to PED file (tab-separated, no header)"
    )
    parser.add_argument(
        "-d", "--directory", type=str, required=True,
        help="Directory containing *hg38_multianno.csv files"
    )
    parser.add_argument(
        "-o", "--output_dir", type=str, required=True,
        help="Output directory for Excel files"
    )
    parser.add_argument(
        "-m", "--omim", type=str,
        default="/path/to/dbs/OMIM_Summary_File.csv",
        help="Path to OMIM summary CSV (default: %(default)s)"
    )
    parser.add_argument(
        "-p", "--hpo", type=str,
        default="/path/to/dbs/hpo_genes_to_phenotype.txt",
        help="Path to HPO genes-to-phenotype TSV (default: %(default)s)"
    )
    parser.add_argument(
        "-g", "--gel", type=str, default=None,
        help="Path to GEL panel CSV with columns Gene.Symbol and GEL_annotation. "
             "Genes not in the panel are labelled 'Not_Present'. (default: not used)"
    )
    parser.add_argument(
        "-r", "--reference", type=str, default="/path/to/reference/hg38.fasta",
        help="Path to indexed reference FASTA (hg38). Enables exact variant "
             "URLs for indels in gnomAD and DECIPHER linkouts. "
             "Requires pysam and a .fai index file. (default: not used)"
    )

    args = parser.parse_args()

    main(
        ped_path=args.ped,
        annovar_dir=args.directory,
        output_dir=args.output_dir,
        omim_path=args.omim,
        hpo_path=args.hpo,
        gel_path=args.gel,
        ref_path=args.reference,
    )
