"""Tests for z-score based flagging."""

import numpy as np
import pandas as pd
import pytest

from methreport.controls import z_score_flag, MIN_INFORMATIVE_CPG


def _make_sample(positions, methylation_pcts):
    return pd.DataFrame({"position": positions, "methylation_pct": methylation_pcts})


def _make_reference(positions, means, sds):
    return pd.DataFrame({"position": positions, "mean": means, "sd": sds})


def test_normal_sample_not_flagged():
    # Sample matches control mean exactly → z=0
    positions = list(range(1000, 1010))
    ref = _make_reference(positions, [50.0] * 10, [5.0] * 10)
    sample = _make_sample(positions, [50.0] * 10)
    flag, z, n = z_score_flag(sample, ref)
    assert flag == "NORMAL"
    assert abs(z) < 0.1
    assert n == 10


def test_low_methylation_flagged():
    # Sample at 10% where controls expect 50% → z ≈ -8
    positions = list(range(1000, 1010))
    ref = _make_reference(positions, [50.0] * 10, [5.0] * 10)
    sample = _make_sample(positions, [10.0] * 10)
    flag, z, n = z_score_flag(sample, ref)
    assert flag == "LOW"
    assert z < -2.0


def test_high_methylation_flagged():
    # Sample at 90% where controls expect 50% → z ≈ +8
    positions = list(range(1000, 1010))
    ref = _make_reference(positions, [50.0] * 10, [5.0] * 10)
    sample = _make_sample(positions, [90.0] * 10)
    flag, z, n = z_score_flag(sample, ref)
    assert flag == "HIGH"
    assert z > 2.0


def test_border_z2_not_flagged():
    # Mean z-score exactly at 2.0 is not flagged (strictly > 2.0 triggers flag)
    positions = list(range(1000, 1010))
    ref = _make_reference(positions, [50.0] * 10, [5.0] * 10)
    # z = (60 - 50) / 5 = 2.0 exactly — should be NORMAL (not > 2.0)
    sample = _make_sample(positions, [60.0] * 10)
    flag, z, n = z_score_flag(sample, ref)
    assert flag == "NORMAL"
    assert abs(z - 2.0) < 0.01


def test_no_overlap_returns_na():
    # Sample positions don't match reference positions
    ref = _make_reference([1000, 1001], [50.0, 50.0], [5.0, 5.0])
    sample = _make_sample([2000, 2001], [30.0, 30.0])
    flag, z, n = z_score_flag(sample, ref)
    assert flag == "NA"
    assert n == 0


def test_empty_reference_returns_na():
    sample = _make_sample([1000, 1001], [30.0, 30.0])
    ref = pd.DataFrame(columns=["position", "mean", "sd"])
    flag, z, n = z_score_flag(sample, ref)
    assert flag == "NA"


def test_too_few_informative_falls_back_to_fixed_threshold():
    # Only 1 overlapping position (< MIN_INFORMATIVE_CPG) and that value is LOW
    ref = _make_reference([1000], [50.0], [5.0])
    sample = _make_sample([1000], [10.0])
    flag, z, n = z_score_flag(sample, ref, min_informative=MIN_INFORMATIVE_CPG)
    # Falls back to fixed threshold: 10% < 30% → LOW
    assert flag == "LOW"
    assert n == 1
    assert np.isnan(z)


def test_zero_sd_clipped_to_min():
    # All controls at exactly 50% → SD=0. Should not divide by zero.
    positions = list(range(1000, 1010))
    ref = _make_reference(positions, [50.0] * 10, [0.0] * 10)
    sample = _make_sample(positions, [10.0] * 10)
    flag, z, n = z_score_flag(sample, ref)
    # SD is clipped to _MIN_SD=2.0, so z = (10-50)/2 = -20 → LOW
    assert flag == "LOW"
    assert z < -2.0
