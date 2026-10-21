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


@dataclass(frozen=True)
class FlagDecision:
    """What both methods concluded, and what they amount to together."""
    flag: str
    z_score: float
    n_informative: int
    flag_controls: str
    flag_asymmetry: str
    hp_delta: float
    confidence: str
    evidence: str


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
    # Dual-mode record: what each method concluded on its own, and how much the
    # two together are worth. Kept separately from `flag` so that a combined
    # call can never hide the evidence it was built from.
    flag_controls: str = "NA"   # z-score against controls, alone
    flag_asymmetry: str = "NA"  # HP1 vs HP2 allelic asymmetry, alone
    hp_delta: float = float("nan")  # |HP1 - HP2| mean methylation, percentage points
    confidence: str = "none"    # see methreport.analysis.CONFIDENCE_* below
    evidence: str = ""          # one line naming the methods that contributed

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
            "flag_controls": self.flag_controls,
            "flag_asymmetry": self.flag_asymmetry,
            "hp_delta": round(self.hp_delta, 1) if not np.isnan(self.hp_delta) else None,
            "confidence": self.confidence,
            "evidence": self.evidence,
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

# How much the combined call is worth. These are recorded on every region so a
# reader can tell a call both methods agree on from one that rests on a single
# method, or on the fixed-threshold fallback.
CONFIDENCE_CONCORDANT = "concordant"       # both methods available and agreeing
CONFIDENCE_DISCORDANT = "discordant"       # both available, disagreeing
CONFIDENCE_CONTROLS = "controls only"      # asymmetry unavailable
CONFIDENCE_ASYMMETRY = "asymmetry only"    # controls unavailable
CONFIDENCE_FALLBACK = "fixed threshold"    # neither available; low confidence
CONFIDENCE_NONE = "none"                   # nothing measurable


def _asymmetry_flag(
    hp1: RegionMethylation,
    hp2: RegionMethylation,
) -> tuple[str, float]:
    """Flag a region from the symmetry of its two haplotypes.

    An imprinted DMR is expected to be strongly asymmetric: one haplotype
    methylated, the other not. Both high is hypermethylation, both low is
    hypomethylation, and one of each is the normal pattern whatever the
    unphased mean happens to be. Returns ("NA", delta) when the region is not
    phased deeply enough to judge, or when the two haplotypes sit between the
    thresholds so that neither asymmetry nor its loss can be asserted.
    """
    m1, m2 = hp1.mean_methylation_pct, hp2.mean_methylation_pct
    if np.isnan(m1) or np.isnan(m2):
        return "NA", float("nan")

    delta = abs(m1 - m2)
    if hp1.n_cpg < _MIN_PHASED_CPG or hp2.n_cpg < _MIN_PHASED_CPG:
        return "NA", delta

    hp_high, hp_low = max(m1, m2), min(m1, m2)
    if hp_high >= _HP_HIGH_THRESH and hp_low <= _HP_LOW_THRESH:
        return "NORMAL", delta          # canonical imprinting pattern
    if hp_low >= _HP_HIGH_THRESH:
        return "HIGH", delta            # both alleles methylated
    if hp_high <= _HP_LOW_THRESH:
        return "LOW", delta             # both alleles unmethylated
    return "NA", delta                  # between thresholds: no assertion


def combine_flags(flag_controls: str, flag_asymmetry: str) -> tuple[str, str, str]:
    """Decide a region from both methods at once.

    Returns (flag, confidence, evidence). The rule where the two disagree is
    deliberately asymmetric: an abnormal result from either method is carried
    through rather than averaged away, because a missed epimutation costs more
    than a region sent for review. The disagreement is recorded in the
    confidence field so that such a call is never mistaken for a settled one.
    """
    have_c = flag_controls in ("NORMAL", "LOW", "HIGH")
    have_a = flag_asymmetry in ("NORMAL", "LOW", "HIGH")

    if have_c and have_a:
        if flag_controls == flag_asymmetry:
            return (flag_controls, CONFIDENCE_CONCORDANT,
                    "z-score against controls and allelic asymmetry agree")
        if flag_controls == "NORMAL":
            return (flag_asymmetry, CONFIDENCE_DISCORDANT,
                    "allelic asymmetry reports "
                    f"{flag_asymmetry.lower()} where the z-score against "
                    "controls is normal")
        if flag_asymmetry == "NORMAL":
            return (flag_controls, CONFIDENCE_DISCORDANT,
                    "the z-score against controls reports "
                    f"{flag_controls.lower()} where allelic asymmetry is normal")
        # Both abnormal, in opposite directions. The z-score is positional and
        # quantitative, so it sets the direction, but the contradiction stands.
        return (flag_controls, CONFIDENCE_DISCORDANT,
                f"the z-score against controls reports {flag_controls.lower()} "
                f"and allelic asymmetry reports {flag_asymmetry.lower()}; the "
                "two contradict each other")

    if have_c:
        return (flag_controls, CONFIDENCE_CONTROLS,
                "z-score against controls; the region is not phased deeply "
                "enough for allelic asymmetry")
    if have_a:
        return (flag_asymmetry, CONFIDENCE_ASYMMETRY,
                "allelic asymmetry; no control data for this region")
    return ("NA", CONFIDENCE_NONE, "")


def _compute_flag(
    unphased: RegionMethylation,
    hp1: RegionMethylation,
    hp2: RegionMethylation,
    reference: pd.DataFrame,
) -> FlagDecision:
    """Flag a region from both independent lines of evidence.

    The two methods answer the same question in different ways and fail in
    different places. A z-score against per-position controls is quantitative
    and sensitive, but needs control data for the region and says nothing about
    which allele carries the methylation. Allelic asymmetry needs no controls
    and speaks directly to the imprinting mechanism, but needs the region to be
    phased and cannot judge haplotypes that sit between its thresholds.

    Both are therefore computed whenever they can be, and `combine_flags`
    decides. Only when neither is available does the fixed threshold on the
    unphased mean apply, which the caller records as low confidence.
    """
    # ── Evidence 1: z-score against per-position controls ──────────────────
    flag_controls, z_score, n_informative = "NA", float("nan"), 0
    if not reference.empty:
        sample_df = unphased.to_dataframe()
        if not sample_df.empty:
            flag_controls, z_score, n_informative = z_score_flag(sample_df, reference)

    # ── Evidence 2: allelic asymmetry between the two haplotypes ───────────
    flag_asymmetry, hp_delta = _asymmetry_flag(hp1, hp2)

    # ── Decide ─────────────────────────────────────────────────────────────
    flag, confidence, evidence = combine_flags(flag_controls, flag_asymmetry)

    # ── Fallback: neither method could speak ───────────────────────────────
    if flag == "NA":
        mean_meth = unphased.mean_methylation_pct
        if np.isnan(mean_meth) or unphased.n_cpg == 0:
            return FlagDecision("NA", z_score, n_informative, flag_controls,
                                flag_asymmetry, hp_delta, CONFIDENCE_NONE,
                                "no methylation measurable in this region")
        flag = flag_result(mean_meth)
        confidence = CONFIDENCE_FALLBACK
        evidence = ("fixed threshold on the unphased mean; neither control "
                    "data nor sufficient phasing was available")

    return FlagDecision(flag, z_score, n_informative, flag_controls,
                        flag_asymmetry, hp_delta, confidence, evidence)


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
            # Still measured and reported, but never flagged: the region is
            # known to mismap, so a departure from the controls would say more
            # about the alignment than about the sample.
            decision = FlagDecision(
                "UNRELIABLE", float("nan"), 0, "NA", "NA", float("nan"),
                CONFIDENCE_NONE, region.unreliable_reason)
        else:
            decision = _compute_flag(
                meth["unphased"], meth["hp1"], meth["hp2"], ref
            )

        result = RegionResult(
            region=region,
            unphased=meth["unphased"],
            hp1=meth["hp1"],
            hp2=meth["hp2"],
            reference=ref,
            flag=decision.flag,
            z_score=decision.z_score,
            n_informative=decision.n_informative,
            flag_controls=decision.flag_controls,
            flag_asymmetry=decision.flag_asymmetry,
            hp_delta=decision.hp_delta,
            confidence=decision.confidence,
            evidence=decision.evidence,
        )
        analysis.results.append(result)

        meth_val = result.mean_methylation
        log.info(
            "    → %.1f%% methylation  flag=%s (%s)  controls=%s  asymmetry=%s"
            "  z=%.2f  informative_cpg=%d/%d  cov=%.1f×",
            meth_val if not np.isnan(meth_val) else 0,
            decision.flag,
            decision.confidence,
            decision.flag_controls,
            decision.flag_asymmetry,
            decision.z_score if not np.isnan(decision.z_score) else 0,
            decision.n_informative,
            result.unphased.n_cpg,
            result.unphased.mean_coverage,
        )

    n_flagged = len(analysis.flagged_regions)
    log.info("Analysis complete: %d/%d regions flagged.", n_flagged, len(regions))
    return analysis
