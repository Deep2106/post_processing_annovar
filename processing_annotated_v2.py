#!/usr/bin/env python3
"""
Post-annotation Excel generator for trio/duo/singleton
Filenames expected: LH0007A.hg38_multianno.csv etc.
Generates one Excel per run: {proband}_analysis.xlsx
"""

import pandas as pd
import numpy as np
import os
import sys
from glob import glob
from collections import defaultdict
import argparse
from pathlib import Path

def get_hpo(hpo_path):
    try:
        hpo = pd.read_csv(hpo_path, sep="\t")
        return hpo[["gene_symbol", "hpo_name"]].groupby("gene_symbol").agg(" | ".join).to_dict(orient="index")
    except Exception as e:
        print(f"Warning: Could not load HPO file: {e}", file=sys.stderr)
        return {}

def get_omim(omim_path):
    try:
        omim = pd.read_csv(omim_path)
        omim["Approved Symbol"] = omim["Approved Symbol"].astype(str).str.upper()
        return omim
    except Exception as e:
        print(f"Warning: Could not load OMIM file: {e}", file=sys.stderr)
        return pd.DataFrame()

def detect_family_members(annovar_dir, family_code, project_name=None):
    annovar_path = Path(annovar_dir).resolve()
    if not annovar_path.is_dir():
        print(f"ERROR: Directory does not exist: {annovar_path}", file=sys.stderr)
        sys.exit(1)

    patterns = []
    if project_name and project_name.strip():
        patterns.append(f"{annovar_dir}/{project_name}*{family_code}[ABC]*.hg38_multianno.csv")

    patterns.extend([
        f"{annovar_dir}/*{family_code}[ABC]*.hg38_multianno.csv",
        f"{annovar_dir}/{family_code}[ABC]*.hg38_multianno.csv",
        f"{annovar_dir}/{family_code}[ABC].hg38_multianno.csv",
    ])

    found_files = []
    for pat in patterns:
        found_files.extend(glob(pat))
    found_files = sorted(set(found_files))

    if not found_files:
        print(f"\nERROR: No multianno files found for {family_code}", file=sys.stderr)
        print("Tried patterns:")
        for p in patterns:
            print(f"  • {p}")
        sys.exit(1)

    members = defaultdict(list)
    for f in found_files:
        basename = Path(f).name
        if f"{family_code}A" in basename:
            members["A"].append(f)
        elif f"{family_code}B" in basename:
            members["B"].append(f)
        elif f"{family_code}C" in basename:
            members["C"].append(f)

    detected = {}
    if members["A"]:
        detected["proband"] = f"{family_code}A"
    if members["B"]:
        detected["father"] = f"{family_code}B"
    if members["C"]:
        detected["mother"] = f"{family_code}C"

    if not detected:
        print("Found files but no A/B/C identified:", file=sys.stderr)
        for f in found_files:
            print(f"  • {Path(f).name}", file=sys.stderr)
        sys.exit(1)

    count = len(detected)
    family_type = "singleton" if count == 1 else "duo" if count == 2 else "trio" if count == 3 else "unknown"

    return detected, family_type, found_files

def print_confirmation_and_ask(annovar_dir, output_dir, omim_path, hpo_path, project_name,
                               family_code, detected, family_type, found_files):
    print("\n" + "="*80)
    print(" " * 25 + "ANNOVAR POST-PROCESSING SUMMARY")
    print("="*80)
    print(f"Family code:          {family_code}")
    print(f"Family type:          ** {family_type.upper()} **")
    print(f"Project prefix:       {project_name or 'None (direct naming)'}")
    print(f"Annovar directory:    {annovar_dir}")
    print(f"Output directory:     {output_dir}")
    print()

    print("Detected samples:")
    for role, sid in detected.items():
        print(f"  • {role.capitalize():<8}: {sid}")

    print(f"\nFound files ({len(found_files)}):")
    for f in sorted(found_files)[:8]:  # limit to avoid flooding terminal
        print(f"  • {Path(f).name}")
    if len(found_files) > 8:
        print(f"  ... and {len(found_files)-8} more")

    print("\nIs everything correct?")
    while True:
        ans = input("yes/no: ").strip().lower()
        if ans in ["yes", "y"]:
            print("\nStarting analysis...\n")
            return True
        if ans in ["no", "n"]:
            print("Aborted by user.")
            return False
        print("Please type yes or no.")

def prepare_data(annovar_dir, sample_id, omim, project_name=None):
    patterns = []
    if project_name:
        patterns.append(f"{annovar_dir}/{project_name}*{sample_id}*.hg38_multianno.csv")
    patterns.extend([
        f"{annovar_dir}/*{sample_id}*.hg38_multianno.csv",
        f"{annovar_dir}/{sample_id}*.hg38_multianno.csv",
        f"{annovar_dir}/{sample_id}.hg38_multianno.csv",
    ])

    candidates = []
    for pat in patterns:
        candidates.extend(glob(pat))
    candidates = sorted(set(candidates))

    if not candidates:
        print(f"ERROR: No multianno file found for {sample_id}", file=sys.stderr)
        sys.exit(1)

    dfs = []
    for path in candidates:
        print(f"  → {Path(path).name}")
        try:
            df = pd.read_csv(path, sep=",", low_memory=False)
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Failed to load {Path(path).name}: {e}", file=sys.stderr)

    if not dfs:
        sys.exit(1)

    df = pd.concat(dfs, ignore_index=True)

    required = ["Chr", "Start", "End", "Ref", "Alt", "GT", "Gene.refGene"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: Missing required columns: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    df["Gene.refGene"] = df["Gene.refGene"].astype(str).str.upper()

    gamze = omim[["Approved Symbol", "Phenotypes"]].groupby("Approved Symbol").agg(" | ".join).to_dict(orient="index")
    df["OMIM"] = df["Gene.refGene"].map(gamze).apply(
        lambda x: x["Phenotypes"].replace("0 |","").replace("| 0","").replace("0 | 0","#N/A")
        if isinstance(x, dict) else "#N/A"
    )

    df["ExonicFunc.refGene"] = df["ExonicFunc.refGene"].fillna("missing_value").replace(".", "missing_value")

    gnomad_cols = [c for c in df.columns if "gnomad" in c.lower() or "exac" in c.lower()]
    if gnomad_cols:
        df[gnomad_cols] = df[gnomad_cols].fillna(0).replace(".", 0).apply(pd.to_numeric, errors='coerce').fillna(0)
        df["BA1"] = df[gnomad_cols].ge(0.01).any(axis=1)
    else:
        df["BA1"] = False

    df["variant"] = df.apply(lambda r: f"{r['Chr']}:{r['Start']}:{r['Ref']}:{r['Alt']}", axis=1)

    return df

def analyze_and_write(proband_df, family_type, output_dir, proband_id):
    # Flags
    proband_df["LOF"] = (
        proband_df["ExonicFunc.refGene"].isin(["frameshift deletion", "frameshift insertion", "stopgain", "startloss", "stoploss"]) |
        proband_df["Func.refGene"].str.contains("splicing", na=False)
    )

    proband_df["intronic"] = proband_df["Func.refGene"] == "intronic"

    # GNOMAD_linkout column
    def make_gnomad_link(variant):
        if pd.isna(variant) or not isinstance(variant, str) or ':' not in variant:
            return ""
        try:
            chrom, pos, ref, alt = variant.split(':')
            formatted = f"{chrom}-{pos}-{ref}-{alt}"
            url = f"https://gnomad.broadinstitute.org/variant/{formatted}?dataset=gnomad_r4"
            return f'=HYPERLINK("{url}", "GNOMELINK")'
        except:
            return ""

    if "variant" in proband_df.columns:
        proband_df["GNOMAD_linkout"] = proband_df["variant"].apply(make_gnomad_link)
    else:
        proband_df["GNOMAD_linkout"] = ""

    # Reorder columns
    cols = list(proband_df.columns)

    # 1. Group after GT: OMIM, BA1, variant, CLNSIG, REVEL
    try:
        gt_pos = cols.index("GT") + 1
    except ValueError:
        gt_pos = 0

    to_insert_after_gt = ["OMIM", "BA1", "variant", "CLNSIG", "REVEL"]
    for col in to_insert_after_gt:
        if col in cols:
            cols.remove(col)

    for col in reversed(to_insert_after_gt):
        if col in proband_df.columns:
            cols.insert(gt_pos, col)

    # 2. GNOMAD_linkout after AAChange.refGene
    try:
        aachange_pos = cols.index("AAChange.refGene") + 1
    except ValueError:
        aachange_pos = len(cols)

    if "GNOMAD_linkout" in cols:
        cols.remove("GNOMAD_linkout")
    cols.insert(aachange_pos, "GNOMAD_linkout")

    # 3. Analysis columns after Caller
    try:
        caller_pos = cols.index("Caller") + 1
    except ValueError:
        caller_pos = len(cols)

    analysis_cols = ["Inheritance", "LOF", "intronic", "HPO_phenotypes"]
    for col in reversed([c for c in analysis_cols if c in cols]):
        cols.remove(col)
        cols.insert(caller_pos, col)

    proband_df = proband_df[cols]

    # ── Sheets ────────────────────────────────────────────────────────
    sheets = {}

    sheets["Raw_Data"] = proband_df.copy()

    #if family_type in ("duo", "trio"):
    #    sheets["DeNovo_Het"] = proband_df.query(
     #       "Inheritance == 'De Novo' and GT == 'het' and BA1 == False and intronic == False"
     #   )

    sheets["rare_Hom"] = proband_df.query(
        "GT == 'hom' and BA1 == False and intronic == False"
    )

    het_clean = proband_df.query("GT == 'het' and BA1 == False and intronic == False")
    gene_counts = het_clean["Gene.refGene"].value_counts()
    comp_genes = gene_counts[gene_counts > 1].index
    sheets["Comp_Het"] = het_clean[het_clean["Gene.refGene"].isin(comp_genes)].sort_values("Gene.refGene")

    sheets["Het"] = proband_df.query("GT == 'het' and BA1 == False and intronic == False")

    sheets["Synonymous"] = proband_df.query(
        "`ExonicFunc.refGene` == 'synonymous SNV' and BA1 == False and intronic == False"
    )

    outfile = Path(output_dir) / f"{proband_id}_analysis.xlsx"
    print(f"\nWriting Excel → {outfile}")

    with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
        for name, df_sheet in sheets.items():
            if df_sheet.empty and name != "Raw_Data":
                print(f"  Skipping empty sheet: {name}")
                continue
            df_sheet.to_excel(writer, sheet_name=name, index=False)

    print("Done!")

def main(family_code, annovar_dir, omim_path, hpo_path, output_dir, project_name):
    detected, family_type, found_files = detect_family_members(annovar_dir, family_code, project_name)

    if not print_confirmation_and_ask(
        annovar_dir, output_dir, omim_path, hpo_path, project_name,
        family_code, detected, family_type, found_files
    ):
        sys.exit(0)

    omim = get_omim(omim_path)
    gizem = get_hpo(hpo_path)

    proband_id = detected["proband"]
    proband_df = prepare_data(annovar_dir, proband_id, omim, project_name)

    proband_df["HPO_phenotypes"] = proband_df["Gene.refGene"].map(gizem).apply(
        lambda x: x.get("hpo_name", "#N/A").replace("0 |","").replace("| 0","").replace("0 | 0","#N/A")
        if isinstance(x, dict) else "#N/A"
    )

    if family_type in ("duo", "trio"):
        mother_df = prepare_data(annovar_dir, detected.get("mother", ""), omim, project_name) if "mother" in detected else pd.DataFrame()
        father_df = prepare_data(annovar_dir, detected.get("father", ""), omim, project_name) if "father" in detected else pd.DataFrame()

        cond = [
            (proband_df["variant"].isin(mother_df["variant"]) & proband_df["variant"].isin(father_df["variant"])),
            proband_df["variant"].isin(mother_df["variant"]),
            proband_df["variant"].isin(father_df["variant"]),
        ]
        proband_df["Inheritance"] = np.select(cond, ["Familial", "Maternal", "Paternal"], default="De Novo")
    else:
        proband_df["Inheritance"] = "Singleton"

    analyze_and_write(proband_df, family_type, output_dir, proband_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Post-annotation Excel generator – A/B/C family logic",
        epilog="""
Examples:
  python this_script.py -s LH0007
  python this_script.py -s LH0007 -P SomePrefix
  python this_script.py -s LH0007 -d /custom/annovar/path -o /output -m /omim/omim.csv -p /hpo/hpo.txt -P prj01
        """
    )
    parser.add_argument("-s", "--sample_code", type=str, required=True, help="e.g. LH0007")
    parser.add_argument("-d", "--directory", type=str, default=None)
    parser.add_argument("-o", "--output_dir", type=str, default=None)
    parser.add_argument("-m", "--omim", type=str, default=None)
    parser.add_argument("-p", "--hpo", type=str, default=None)
    parser.add_argument("-P", "--project_name", type=str, default=None,
                        help="Optional prefix (leave empty if files are like LH0007A.hg38_multianno.csv)")

    args = parser.parse_args()

    main(
        family_code=args.sample_code,
        annovar_dir=args.directory,
        omim_path=args.omim,
        hpo_path=args.hpo,
        output_dir=args.output_dir,
        project_name=args.project_name
    )
