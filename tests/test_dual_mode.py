"""Tests for the dual-mode flag decision.

Two independent lines of evidence are combined: a z-score against per-position
controls, and the allelic asymmetry between the two haplotypes. These tests pin
the combination rule, which is deliberately asymmetric — an abnormal result from
either method survives the combination rather than being averaged away — and
pin the asymmetry call itself.
"""

import numpy as np
import pytest

from methreport.analysis import (
    CONFIDENCE_ASYMMETRY,
    CONFIDENCE_CONCORDANT,
    CONFIDENCE_CONTROLS,
    CONFIDENCE_DISCORDANT,
    CONFIDENCE_NONE,
    _MIN_PHASED_CPG,
    _asymmetry_flag,
    combine_flags,
)
from methreport.reader import CpGSite, RegionMethylation
from methreport.regions import get_regions


def _methylation(pcts, haplotype="hp1", coverage=100):
    """A RegionMethylation carrying one CpG per given percentage.

    Coverage is 100 so that a requested percentage is reachable exactly:
    CpGSite counts reads, so at 30x the nearest achievable value to 95%
    is 93.3%, which would make the expected deltas below approximate for
    no good reason.
    """
    region = get_regions("hg38")[0]
    sites = [
        CpGSite(
            chrom=region.chrom,
            position=region.start + i,
            n_total=coverage,
            n_methyl=int(round(coverage * pct / 100.0)),
        )
        for i, pct in enumerate(pcts)
    ]
    return RegionMethylation(region=region, haplotype=haplotype, sites=sites)


# ─── the combination rule ──────────────────────────────────────────────────

@pytest.mark.parametrize("flag", ["NORMAL", "LOW", "HIGH"])
def test_agreement_is_concordant(flag):
    combined, confidence, evidence = combine_flags(flag, flag)
    assert combined == flag
    assert confidence == CONFIDENCE_CONCORDANT
    assert evidence


@pytest.mark.parametrize("abnormal", ["LOW", "HIGH"])
def test_abnormal_asymmetry_survives_a_normal_zscore(abnormal):
    """A normal z-score must not mask loss of allelic asymmetry."""
    combined, confidence, _ = combine_flags("NORMAL", abnormal)
    assert combined == abnormal
    assert confidence == CONFIDENCE_DISCORDANT


@pytest.mark.parametrize("abnormal", ["LOW", "HIGH"])
def test_abnormal_zscore_survives_normal_asymmetry(abnormal):
    """Symmetric haplotypes must not mask a departure from the controls."""
    combined, confidence, _ = combine_flags(abnormal, "NORMAL")
    assert combined == abnormal
    assert confidence == CONFIDENCE_DISCORDANT


def test_opposite_directions_take_the_zscore_and_are_marked_discordant():
    combined, confidence, evidence = combine_flags("LOW", "HIGH")
    assert combined == "LOW"
    assert confidence == CONFIDENCE_DISCORDANT
    assert "contradict" in evidence


def test_single_method_is_recorded_as_such():
    assert combine_flags("NORMAL", "NA")[1] == CONFIDENCE_CONTROLS
    assert combine_flags("NA", "NORMAL")[1] == CONFIDENCE_ASYMMETRY
    assert combine_flags("LOW", "NA")[0] == "LOW"
    assert combine_flags("NA", "HIGH")[0] == "HIGH"


def test_no_evidence_at_all():
    combined, confidence, evidence = combine_flags("NA", "NA")
    assert combined == "NA"
    assert confidence == CONFIDENCE_NONE
    assert evidence == ""


# ─── the asymmetry call ────────────────────────────────────────────────────

def test_canonical_imprinting_is_normal():
    n = _MIN_PHASED_CPG
    flag, delta = _asymmetry_flag(_methylation([95.0] * n), _methylation([3.0] * n))
    assert flag == "NORMAL"
    assert delta == pytest.approx(92.0, abs=1.0)


def test_both_haplotypes_methylated_is_high():
    n = _MIN_PHASED_CPG
    flag, delta = _asymmetry_flag(_methylation([90.0] * n), _methylation([85.0] * n))
    assert flag == "HIGH"
    assert delta < 10.0


def test_both_haplotypes_unmethylated_is_low():
    n = _MIN_PHASED_CPG
    flag, _ = _asymmetry_flag(_methylation([4.0] * n), _methylation([8.0] * n))
    assert flag == "LOW"


def test_intermediate_haplotypes_assert_nothing():
    """A normal unphased mean can hide two half-methylated alleles."""
    n = _MIN_PHASED_CPG
    flag, delta = _asymmetry_flag(_methylation([52.0] * n), _methylation([48.0] * n))
    assert flag == "NA"
    assert delta == pytest.approx(4.0, abs=1.0)


def test_shallow_phasing_asserts_nothing_but_still_reports_delta():
    flag, delta = _asymmetry_flag(
        _methylation([95.0] * (_MIN_PHASED_CPG - 1)),
        _methylation([3.0] * (_MIN_PHASED_CPG - 1)),
    )
    assert flag == "NA"
    assert not np.isnan(delta)


def test_unphased_region_asserts_nothing():
    flag, delta = _asymmetry_flag(_methylation([]), _methylation([]))
    assert flag == "NA"
    assert np.isnan(delta)
