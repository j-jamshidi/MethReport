"""Read methylation data from modkit BEDMethyl files.

This is the recommended primary input for MethReport. Modkit (ONT's official
tool) correctly handles:
  - Strand-aware CpG merging (C on + strand and complement C on - strand
    are merged into a single CpG call)
  - Probability calibration for Dorado basecaller output
  - MM/ML encoding variants across basecaller versions

BEDMethyl format (modkit pileup output):
  chrom  start  end  modtype  score  strand  thick_start  thick_end
  color  N_valid_cov  pct_modified

Column indices (0-based):
  [0] chrom
  [1] start (0-based)
  [3] modtype — we select "m" (5mC)
  [9] N_valid_cov
  [10] pct_modified (0–100)

Phased analysis requires three separate files: unphased, hp1, hp2.
These are produced by running modkit with --partition-tag HP:
  modkit pileup sample.bam unphased.bed --ref genome.fa
  modkit pileup sample.bam hp1.bed --ref genome.fa --partition-tag HP --filter-tag HP:1
  modkit pileup sample.bam hp2.bed --ref genome.fa --partition-tag HP --filter-tag HP:2

Or from wf-human-variation (--mod --phased flags), which outputs:
  sample.methyl.cpg.acc.bed        (unphased)
  sample.methyl.cpg.hp1.acc.bed    (HP1)
  sample.methyl.cpg.hp2.acc.bed    (HP2)
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from methreport.reader import CpGSite, RegionMethylation
from methreport.regions import DMRegion

log = logging.getLogger(__name__)

_MOD_TYPE_5MC = "m"
_COL_CHROM   = 0
_COL_START   = 1
_COL_MODTYPE = 3
_COL_COV     = 9
_COL_PCT     = 10


def _load_bed(path: Path) -> pd.DataFrame:
    """Load a BEDMethyl file into a DataFrame, filtering to 5mC rows only."""
    rows = []
    n_skipped = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= _COL_PCT:
                n_skipped += 1
                continue
            if parts[_COL_MODTYPE] != _MOD_TYPE_5MC:
                continue
            try:
                rows.append((
                    parts[_COL_CHROM],
                    int(parts[_COL_START]),
                    int(parts[_COL_COV]),
                    float(parts[_COL_PCT]),
                ))
            except (ValueError, IndexError):
                n_skipped += 1

    if n_skipped:
        log.debug("%s: skipped %d malformed lines", path.name, n_skipped)

    if not rows:
        return pd.DataFrame(columns=["chrom", "position", "n_valid_cov", "pct_modified"])

    df = pd.DataFrame(rows, columns=["chrom", "position", "n_valid_cov", "pct_modified"])
    log.debug("%s: loaded %d 5mC sites", path.name, len(df))
    return df


def _chrom_variants(chrom: str) -> list[str]:
    """Return both 'chr1' and '1' forms to tolerate mixed conventions."""
    if chrom.startswith("chr"):
        return [chrom, chrom[3:]]
    return [chrom, f"chr{chrom}"]


def extract_region_from_bed(
    bed_df: pd.DataFrame,
    region: DMRegion,
    haplotype: str,
    min_coverage: int,
) -> RegionMethylation:
    """Extract CpG sites for one DMR from a pre-loaded BED DataFrame."""
    chroms = _chrom_variants(region.chrom)
    mask = (
        bed_df["chrom"].isin(chroms)
        & (bed_df["position"] >= region.start)
        & (bed_df["position"] < region.end)
        & (bed_df["n_valid_cov"] >= min_coverage)
    )
    subset = bed_df[mask].copy()

    sites = []
    for _, row in subset.iterrows():
        cov = int(row["n_valid_cov"])
        pct = float(row["pct_modified"])
        n_methyl = round(pct * cov / 100.0)
        sites.append(CpGSite(
            chrom=region.chrom,
            position=int(row["position"]),
            n_total=cov,
            n_methyl=n_methyl,
        ))

    return RegionMethylation(region=region, haplotype=haplotype, sites=sites)


class BedMethylReader:
    """Loads BEDMethyl files once and serves per-region queries efficiently."""

    def __init__(
        self,
        unphased_path: Path | str,
        hp1_path: Path | str | None = None,
        hp2_path: Path | str | None = None,
    ) -> None:
        self._dfs: dict[str, pd.DataFrame] = {}

        unphased_path = Path(unphased_path)
        if not unphased_path.exists():
            raise FileNotFoundError(f"Unphased BED not found: {unphased_path}")
        log.info("Loading BEDMethyl: %s", unphased_path)
        self._dfs["unphased"] = _load_bed(unphased_path)

        for hp_key, path in [("hp1", hp1_path), ("hp2", hp2_path)]:
            if path is not None:
                path = Path(path)
                if not path.exists():
                    raise FileNotFoundError(f"BED not found ({hp_key}): {path}")
                log.info("Loading BEDMethyl (%s): %s", hp_key, path)
                self._dfs[hp_key] = _load_bed(path)
            else:
                self._dfs[hp_key] = pd.DataFrame(
                    columns=["chrom", "position", "n_valid_cov", "pct_modified"]
                )

    def extract(
        self,
        region: DMRegion,
        min_coverage: int = 5,
    ) -> dict[str, RegionMethylation]:
        """Return {unphased, hp1, hp2} RegionMethylation for a DMR.

        HP tracks use a lower coverage threshold (ceil(min_coverage / 2)) because
        phasing splits reads across two haplotypes, roughly halving per-HP depth.
        """
        min_cov_hp = max(1, (min_coverage + 1) // 2)
        thresholds = {"unphased": min_coverage, "hp1": min_cov_hp, "hp2": min_cov_hp}
        return {
            hp: extract_region_from_bed(self._dfs[hp], region, hp, thresholds[hp])
            for hp in ("unphased", "hp1", "hp2")
        }

    @property
    def is_phased(self) -> bool:
        return len(self._dfs.get("hp1", pd.DataFrame())) > 0
