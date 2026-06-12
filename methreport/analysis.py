"""Orchestrates per-sample methylation analysis across all DMRs.

Ties together reader → controls → summary for one sample BAM.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from methreport.controls import compute_reference_ranges, flag_result, merge_controls
from methreport.reader import (
    DEFAULT_CALL_THRESHOLD,
    DEFAULT_MIN_COVERAGE,
    RegionMethylation,
    extract_region_methylation,
)
from methreport.regions import DMRegion, get_regions

log = logging.getLogger(__name__)


@dataclass
class RegionResult:
    region: DMRegion
    unphased: RegionMethylation
    hp1: RegionMethylation
    hp2: RegionMethylation
    reference: pd.DataFrame   # per-position control stats
    flag: str = "NA"          # "NORMAL" | "LOW" | "HIGH" | "NA"

    @property
    def mean_methylation(self) -> float:
        return self.unphased.mean_methylation_pct

    @property
    def mean_methylation_hp1(self) -> float:
        return self.hp1.mean_methylation_pct

    @property
    def mean_methylation_hp2(self) -> float:
        return self.hp2.mean_methylation_pct

    @property
    def is_phased(self) -> bool:
        return self.hp1.n_cpg > 0 or self.hp2.n_cpg > 0

    def summary_row(self) -> dict:
        return {
            "region": self.region.label,
            "disorder": self.region.disorder,
            "chrom": self.region.chrom,
            "start": self.region.start,
            "end": self.region.end,
            "n_cpg": self.unphased.n_cpg,
            "mean_coverage": round(self.unphased.mean_coverage, 1),
            "mean_coverage_hp1": round(self.hp1.mean_coverage, 1),
            "mean_coverage_hp2": round(self.hp2.mean_coverage, 1),
            "methylation_pct": round(self.mean_methylation, 1) if not np.isnan(self.mean_methylation) else None,
            "methylation_pct_hp1": round(self.mean_methylation_hp1, 1) if not np.isnan(self.mean_methylation_hp1) else None,
            "methylation_pct_hp2": round(self.mean_methylation_hp2, 1) if not np.isnan(self.mean_methylation_hp2) else None,
            "flag": self.flag,
        }


@dataclass
class SampleAnalysis:
    sample_id: str
    bam_path: Path
    genome: str
    results: list[RegionResult] = field(default_factory=list)

    def summary_table(self) -> pd.DataFrame:
        rows = [r.summary_row() for r in self.results]
        return pd.DataFrame(rows)

    @property
    def flagged_regions(self) -> list[RegionResult]:
        return [r for r in self.results if r.flag in ("LOW", "HIGH")]


def run_analysis(
    bam_path: Path | str,
    genome: str = "t2t",
    sample_id: str | None = None,
    user_controls_path: Path | str | None = None,
    replace_controls: bool = False,
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    min_coverage: int = DEFAULT_MIN_COVERAGE,
) -> SampleAnalysis:
    """Run full methylation analysis for a sample BAM.

    Parameters
    ----------
    bam_path:
        Path to the indexed modbam file.
    genome:
        Reference genome used ('t2t' or 'hg38').
    sample_id:
        Name shown in the report. Defaults to BAM filename stem.
    user_controls_path:
        Optional user-supplied controls TSV/XLSX.
    replace_controls:
        If True, replace bundled controls with user-supplied ones.
    call_threshold:
        Probability threshold to call a base as methylated (0–1).
    min_coverage:
        Minimum read depth at a CpG to include it.
    """
    bam_path = Path(bam_path)
    if sample_id is None:
        sample_id = bam_path.stem

    log.info("Starting analysis: sample=%s  genome=%s  bam=%s", sample_id, genome, bam_path)

    controls = merge_controls(genome, user_controls_path, replace=replace_controls)
    regions = get_regions(genome)

    analysis = SampleAnalysis(sample_id=sample_id, bam_path=bam_path, genome=genome)

    for region in regions:
        log.info("  Processing region: %s", region.label)
        meth = extract_region_methylation(
            bam_path=bam_path,
            region=region,
            call_threshold=call_threshold,
            min_coverage=min_coverage,
        )
        ref = compute_reference_ranges(controls, region.name)
        flag = flag_result(meth["unphased"].mean_methylation_pct)

        result = RegionResult(
            region=region,
            unphased=meth["unphased"],
            hp1=meth["hp1"],
            hp2=meth["hp2"],
            reference=ref,
            flag=flag,
        )
        analysis.results.append(result)
        log.info(
            "    → %.1f%% methylation (%s), %d CpGs, coverage %.1fx",
            result.mean_methylation if not np.isnan(result.mean_methylation) else 0,
            flag,
            result.unphased.n_cpg,
            result.unphased.mean_coverage,
        )

    n_flagged = len(analysis.flagged_regions)
    log.info("Analysis complete. %d/%d regions flagged.", n_flagged, len(regions))
    return analysis
