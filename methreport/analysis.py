"""Orchestrates per-sample methylation analysis across all DMRs.

Supports two input modes:
  - BEDMethyl (recommended): modkit pileup output — strand-aware, calibrated
  - modbam: direct BAM reading via MM/ML tag parser (fallback)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from methreport.controls import compute_reference_ranges, flag_result, merge_controls, z_score_flag
from methreport.reader import (
    DEFAULT_CALL_THRESHOLD,
    DEFAULT_MIN_COVERAGE,
    RegionMethylation,
)
from methreport.regions import DMRegion, get_regions

log = logging.getLogger(__name__)


@dataclass
class RegionResult:
    region: DMRegion
    unphased: RegionMethylation
    hp1: RegionMethylation
    hp2: RegionMethylation
    reference: pd.DataFrame     # per-position control stats (informative CpGs only)
    flag: str = "NA"            # "NORMAL" | "LOW" | "HIGH" | "NA"
    z_score: float = float("nan")   # mean z-score across informative CpGs
    n_informative: int = 0      # number of CpGs used for z-score (subset of n_cpg)

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
        meth = self.mean_methylation
        hp1 = self.mean_methylation_hp1
        hp2 = self.mean_methylation_hp2
        return {
            "region": self.region.label,
            "disorder": self.region.disorder,
            "chrom": self.region.chrom,
            "start": self.region.start,
            "end": self.region.end,
            "n_cpg": self.unphased.n_cpg,
            "n_informative_cpg": self.n_informative,
            "mean_coverage": round(self.unphased.mean_coverage, 1),
            "mean_coverage_hp1": round(self.hp1.mean_coverage, 1),
            "mean_coverage_hp2": round(self.hp2.mean_coverage, 1),
            "methylation_pct": round(meth, 1) if not np.isnan(meth) else None,
            "methylation_pct_hp1": round(hp1, 1) if not np.isnan(hp1) else None,
            "methylation_pct_hp2": round(hp2, 1) if not np.isnan(hp2) else None,
            "z_score": round(self.z_score, 2) if not np.isnan(self.z_score) else None,
            "flag": self.flag,
            "unreliable": self.region.unreliable,
            "unreliable_reason": self.region.unreliable_reason if self.region.unreliable else "",
        }


@dataclass
class SampleAnalysis:
    sample_id: str
    input_path: Path        # BAM or unphased BED
    genome: str
    input_mode: str = "bam"  # "bam" | "bed"
    results: list[RegionResult] = field(default_factory=list)

    def summary_table(self) -> pd.DataFrame:
        rows = [r.summary_row() for r in self.results]
        return pd.DataFrame(rows)

    @property
    def flagged_regions(self) -> list[RegionResult]:
        return [r for r in self.results if r.flag in ("LOW", "HIGH") and not r.region.unreliable]

    @property
    def unreliable_regions(self) -> list[RegionResult]:
        return [r for r in self.results if r.region.unreliable]


_MIN_PHASED_CPG = 5   # per-haplotype CpG minimum for phased flagging
_HP_HIGH_THRESH = 65.0  # clearly methylated allele
_HP_LOW_THRESH  = 35.0  # clearly unmethylated allele


def _compute_flag(
    unphased: RegionMethylation,
    hp1: RegionMethylation,
    hp2: RegionMethylation,
    reference: pd.DataFrame,
) -> tuple[str, float, int]:
    """Compute flag using the best available evidence, in priority order:

    1. z-score against per-position controls  (most sensitive, avoids fixed cutoffs)
    2. Phased allelic symmetry (HP1 vs HP2)   (reliable when phased, no controls needed)
    3. Fixed threshold on unphased mean        (last resort — low confidence)

    For imprinted regions, the expected pattern is strong allelic asymmetry:
    one haplotype ~75-100% methylated, the other ~0-25%. Any sample where BOTH
    haplotypes are high (biallelic methylation) or BOTH are low (biallelic
    unmethylation) is flagged. If one is clearly high and the other clearly low,
    the region is NORMAL regardless of the unphased mean.
    """
    # ── Strategy 1: z-score vs controls ────────────────────────────────────
    if not reference.empty:
        sample_df = unphased.to_dataframe()
        if not sample_df.empty:
            return z_score_flag(sample_df, reference)

    # ── Strategy 2: phased allelic symmetry ────────────────────────────────
    if hp1.n_cpg >= _MIN_PHASED_CPG and hp2.n_cpg >= _MIN_PHASED_CPG:
        m1 = hp1.mean_methylation_pct
        m2 = hp2.mean_methylation_pct
        if not (np.isnan(m1) or np.isnan(m2)):
            hp_high = max(m1, m2)
            hp_low  = min(m1, m2)

            if hp_high >= _HP_HIGH_THRESH and hp_low <= _HP_LOW_THRESH:
                # Canonical imprinting pattern: one allele methylated, one not
                return "NORMAL", float("nan"), 0

            if hp_low >= _HP_HIGH_THRESH:
                # Both alleles heavily methylated
                return "HIGH", float("nan"), 0

            if hp_high <= _HP_LOW_THRESH:
                # Both alleles unmethylated
                return "LOW", float("nan"), 0
            # Ambiguous phasing — fall through to fixed threshold

    # ── Strategy 3: fixed threshold (last resort) ───────────────────────────
    mean_meth = unphased.mean_methylation_pct
    if np.isnan(mean_meth) or unphased.n_cpg == 0:
        return "NA", float("nan"), 0
    return flag_result(mean_meth), float("nan"), 0


def run_analysis(
    genome: str = "t2t",
    sample_id: str | None = None,
    # modbam input
    bam_path: Path | str | None = None,
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    # BEDMethyl input
    bed_unphased: Path | str | None = None,
    bed_hp1: Path | str | None = None,
    bed_hp2: Path | str | None = None,
    # shared
    min_coverage: int = DEFAULT_MIN_COVERAGE,
    user_controls_path: Path | str | None = None,
    replace_controls: bool = False,
) -> SampleAnalysis:
    """Run full methylation analysis for a sample.

    Accepts either a modbam file (bam_path) or pre-processed modkit BEDMethyl
    files (bed_unphased + optionally bed_hp1/bed_hp2). BEDMethyl is strongly
    recommended as it avoids strand-merging artefacts from raw MM/ML parsing.

    Parameters
    ----------
    genome:
        Reference genome ('t2t' or 'hg38').
    sample_id:
        Label used in the report. Defaults to filename stem of the input.
    bam_path:
        Path to an indexed modbam file. Mutually exclusive with bed_unphased.
    call_threshold:
        MM/ML probability threshold for methylation calls (modbam mode only).
    bed_unphased:
        Path to modkit BEDMethyl file for all reads. Enables BED input mode.
    bed_hp1 / bed_hp2:
        Optional modkit BEDMethyl files for HP1 and HP2 (BED mode only).
    min_coverage:
        Minimum read depth at a CpG site to include it.
    user_controls_path:
        Optional user-supplied controls TSV/XLSX to supplement or replace bundled data.
    replace_controls:
        If True, bundled controls are replaced rather than supplemented.
    """
    if bed_unphased is not None and bam_path is not None:
        raise ValueError("Provide either bam_path or bed_unphased, not both.")
    if bed_unphased is None and bam_path is None:
        raise ValueError("Provide either bam_path or bed_unphased.")

    # Determine input mode and label
    use_bed = bed_unphased is not None
    input_path = Path(bed_unphased) if use_bed else Path(bam_path)
    if sample_id is None:
        sample_id = input_path.stem

    log.info(
        "Starting analysis: sample=%s  genome=%s  mode=%s  input=%s",
        sample_id, genome, "bed" if use_bed else "bam", input_path,
    )

    controls = merge_controls(genome, user_controls_path, replace=replace_controls)
    regions = get_regions(genome)

    analysis = SampleAnalysis(
        sample_id=sample_id,
        input_path=input_path,
        genome=genome,
        input_mode="bed" if use_bed else "bam",
    )

    # Pre-extract all regions at once (more efficient, enables strand-aware BAM pileup)
    bed_reader = None
    bam_all_regions: dict | None = None

    if use_bed:
        from methreport.bed_reader import BedMethylReader
        bed_reader = BedMethylReader(
            unphased_path=bed_unphased,
            hp1_path=bed_hp1,
            hp2_path=bed_hp2,
        )
    else:
        from methreport.modbam_pileup import extract_from_bam
        log.info("Extracting methylation from BAM (strand-normalised, HP-aware)...")
        bam_all_regions = extract_from_bam(
            bam_path=bam_path,
            genome=genome,
            min_coverage=min_coverage,
            call_threshold=call_threshold,
        )

    for region in regions:
        log.info("  Processing region: %s", region.label)

        # --- Extract methylation ---
        if use_bed:
            meth = bed_reader.extract(region, min_coverage=min_coverage)
        else:
            meth = bam_all_regions[region.name]

        # --- Reference ranges (informative CpGs only) ---
        ref = compute_reference_ranges(controls, region.name)

        # --- Flagging ---
        if region.unreliable:
            flag, z_score, n_informative = "UNRELIABLE", float("nan"), 0
        else:
            flag, z_score, n_informative = _compute_flag(
                meth["unphased"], meth["hp1"], meth["hp2"], ref
            )

        result = RegionResult(
            region=region,
            unphased=meth["unphased"],
            hp1=meth["hp1"],
            hp2=meth["hp2"],
            reference=ref,
            flag=flag,
            z_score=z_score,
            n_informative=n_informative,
        )
        analysis.results.append(result)

        meth_val = result.mean_methylation
        log.info(
            "    → %.1f%% methylation  flag=%s  z=%.2f  informative_cpg=%d/%d  cov=%.1f×",
            meth_val if not np.isnan(meth_val) else 0,
            flag,
            z_score if not np.isnan(z_score) else 0,
            n_informative,
            result.unphased.n_cpg,
            result.unphased.mean_coverage,
        )

    n_flagged = len(analysis.flagged_regions)
    log.info("Analysis complete: %d/%d regions flagged.", n_flagged, len(regions))
    return analysis
