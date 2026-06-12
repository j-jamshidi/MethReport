#!/usr/bin/env python3
"""Convert NanoImprint Excel control files to the TSV format expected by MethReport.

Usage:
    python scripts/convert_controls.py \\
        --t2t  path/to/ctrls_2.xlsx \\
        --hg38 path/to/ctrls_hg38.xlsx \\
        --out  methreport/data/

The Excel files are from carolinehey/NanoImprint.
"""

import argparse
from pathlib import Path

import pandas as pd


REGION_COLUMN_MAP = {
    # Each sheet in the Excel corresponds to a region.
    # Keys are sheet name substrings → region_name values used in MethReport.
    "PLAGL1": "PLAGL1_alt_TSS_DMR",
    "GRB10": "GRB10_alt_TSS_DMR",
    "MEST": "MEST_alt_TSS_DMR",
    "H19": "H19_IGF2_IG_DMR",
    "KCNQ1": "KCNQ1OT1_TSS_DMR",
    "MEG3": "MEG3_TSS_DMR",
    "MEG8": "MEG8_Int2_DMR",
    "MAGEL2": "MAGEL2_TSS_DMR",
    "SNURF": "SNURF_TSS_DMR",
    "PEG3": "PEG3_TSS_DMR",
    "NESP": "GNAS_NESP_TSS_DMR",
    "AS1": "GNAS_AS1_TSS_DMR",
    "GNASXL": "GNASXL_Ex1_DMR",
    "AB": "GNAS_AB_TSS_DMR",
}


def sheet_to_region(sheet_name: str) -> str | None:
    for key, region in REGION_COLUMN_MAP.items():
        if key.upper() in sheet_name.upper():
            return region
    return None


def convert_excel(xlsx_path: Path, genome: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    xl = pd.ExcelFile(xlsx_path)
    frames = []

    for sheet in xl.sheet_names:
        region_name = sheet_to_region(sheet)
        if region_name is None:
            print(f"  Skipping sheet '{sheet}' (no region match)")
            continue

        df = xl.parse(sheet)
        # Expected columns in NanoImprint controls: position, then one column per sample
        # Melt into long format
        if "position" not in [c.lower() for c in df.columns]:
            print(f"  Skipping sheet '{sheet}' (no 'position' column)")
            continue

        pos_col = next(c for c in df.columns if c.lower() == "position")
        sample_cols = [c for c in df.columns if c != pos_col]

        melted = df.melt(id_vars=pos_col, value_vars=sample_cols,
                         var_name="sample_id", value_name="methylation_pct")
        melted = melted.rename(columns={pos_col: "position"})
        melted["region_name"] = region_name
        melted = melted.dropna(subset=["methylation_pct"])
        melted = melted[["region_name", "position", "sample_id", "methylation_pct"]]
        frames.append(melted)
        print(f"  {sheet} → {region_name}: {len(melted)} rows")

    if not frames:
        print(f"No valid sheets found in {xlsx_path}")
        return

    combined = pd.concat(frames, ignore_index=True)
    out_path = out_dir / f"controls_{genome}.tsv"
    combined.to_csv(out_path, sep="\t", index=False)
    print(f"Written: {out_path}  ({len(combined)} rows total)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--t2t", type=Path, help="ctrls_2.xlsx (T2T genome)")
    parser.add_argument("--hg38", type=Path, help="ctrls_hg38.xlsx")
    parser.add_argument("--out", type=Path, default=Path("methreport/data"), help="Output directory")
    args = parser.parse_args()

    if args.t2t:
        print(f"Converting T2T controls: {args.t2t}")
        convert_excel(args.t2t, "t2t", args.out)
    if args.hg38:
        print(f"Converting hg38 controls: {args.hg38}")
        convert_excel(args.hg38, "hg38", args.out)

    if not args.t2t and not args.hg38:
        parser.print_help()


if __name__ == "__main__":
    main()
