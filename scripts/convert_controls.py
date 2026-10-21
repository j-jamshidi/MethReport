#!/usr/bin/env python3
"""Convert the NanoImprint control spreadsheets into MethReport control TSVs.

The published NanoImprint control files (ctrls_2.xlsx for T2T-CHM13v2.0,
ctrls_hg38.xlsx for hg38, from carolinehey/NanoImprint) hold one sheet of
already-aggregated statistics:

    <index>  start_position  methylation_frequency  sd

There is no per-sample column, no chromosome and no control count, so the
per-sample long format MethReport also accepts cannot be reconstructed from
them, and this script does not try: it writes the aggregated schema instead,

    region_name  position  mean  sd  n_controls

leaving n_controls empty because the source does not record it. Positions are
assigned to a DMR by coordinate, using MethReport's own region definitions so
the two cannot drift apart. A position that falls in no region is reported and
dropped; a position that falls in more than one is an error, since the regions
are meant to be disjoint.

Usage:
    python scripts/convert_controls.py \\
        --t2t  path/to/ctrls_2.xlsx \\
        --hg38 path/to/ctrls_hg38.xlsx \\
        --out  methreport/data/

The control data is the intellectual property of the NanoImprint authors and
carries no licence in that repository. Convert it for local use; do not
redistribute the result without their permission.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from methreport.regions import get_regions

# Column names accepted in the source sheet, lowercased.
_POSITION_COLUMNS = ("start_position", "position", "start")
_MEAN_COLUMNS = ("methylation_frequency", "mean", "methylation_pct")
_SD_COLUMNS = ("sd", "std", "stdev")


def _pick(columns: list[str], candidates: tuple[str, ...], label: str) -> str:
    lowered = {str(c).strip().lower(): c for c in columns}
    for name in candidates:
        if name in lowered:
            return lowered[name]
    raise SystemExit(
        f"error: no {label} column found. Looked for {list(candidates)}, "
        f"sheet has {columns}"
    )


def convert_excel(xlsx_path: Path, genome: str, out_dir: Path) -> Path:
    regions = get_regions(genome)
    frame = pd.read_excel(xlsx_path)

    pos_col = _pick(list(frame.columns), _POSITION_COLUMNS, "position")
    mean_col = _pick(list(frame.columns), _MEAN_COLUMNS, "methylation")
    sd_col = _pick(list(frame.columns), _SD_COLUMNS, "standard deviation")

    frame = frame[[pos_col, mean_col, sd_col]].copy()
    frame.columns = ["position", "mean", "sd"]
    frame = frame.dropna(subset=["position", "mean"])
    frame["position"] = frame["position"].astype(float).astype(int)

    def assign(position: int) -> str | None:
        hits = [r.name for r in regions if r.start <= position <= r.end]
        if len(hits) > 1:
            raise SystemExit(
                f"error: position {position} falls in {len(hits)} regions "
                f"({', '.join(hits)}); the region definitions overlap"
            )
        return hits[0] if hits else None

    frame["region_name"] = frame["position"].map(assign)
    outside = frame["region_name"].isna().sum()
    frame = frame.dropna(subset=["region_name"])
    frame["n_controls"] = pd.NA
    frame = frame[["region_name", "position", "mean", "sd", "n_controls"]]
    frame = frame.sort_values(["region_name", "position"]).reset_index(drop=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"controls_{genome}.tsv"
    frame.to_csv(out_path, sep="\t", index=False)

    print(f"{xlsx_path.name} -> {out_path}")
    print(f"  {len(frame)} positions assigned across "
          f"{frame['region_name'].nunique()} of {len(regions)} regions")
    if outside:
        print(f"  {outside} positions fell outside every MethReport region and "
              f"were dropped")
    for region in regions:
        n = int((frame["region_name"] == region.name).sum())
        if n == 0:
            print(f"  WARNING: no control positions for {region.label}")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--t2t", type=Path, help="ctrls_2.xlsx (T2T-CHM13v2.0)")
    parser.add_argument("--hg38", type=Path, help="ctrls_hg38.xlsx")
    parser.add_argument("--out", type=Path, default=Path("methreport/data"),
                        help="output directory")
    args = parser.parse_args()

    if not args.t2t and not args.hg38:
        parser.print_help()
        sys.exit(1)
    if args.t2t:
        convert_excel(args.t2t, "t2t", args.out)
    if args.hg38:
        convert_excel(args.hg38, "hg38", args.out)


if __name__ == "__main__":
    main()
