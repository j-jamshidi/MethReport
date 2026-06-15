"""Tests for control data utilities."""

import pandas as pd
import numpy as np
import pytest
from methreport.controls import compute_reference_ranges, flag_result


def _make_controls(region: str, n_samples: int = 10, mean: float = 50.0, sd: float = 5.0) -> pd.DataFrame:
    rng = np.random.default_rng(42)
    positions = list(range(1000, 1010))
    rows = []
    for pos in positions:
        for i in range(n_samples):
            rows.append({
                "region_name": region,
                "position": pos,
                "sample_id": f"ctrl_{i}",
                "methylation_pct": float(np.clip(rng.normal(mean, sd), 0, 100)),
            })
    return pd.DataFrame(rows)


def test_compute_reference_ranges_basic():
    df = _make_controls("TEST_REGION", mean=50.0, sd=3.0)
    ref = compute_reference_ranges(df, "TEST_REGION")
    assert len(ref) == 10
    assert "mean" in ref.columns
    assert "sd" in ref.columns
    assert all(ref["mean"].between(40, 60))


def test_compute_reference_ranges_filters_non_imprinted():
    # Controls with mean ~80% should be excluded (not imprinted range)
    df = _make_controls("TEST_REGION_HI", mean=80.0, sd=3.0)
    ref = compute_reference_ranges(df, "TEST_REGION_HI")
    assert len(ref) == 0


def test_compute_reference_ranges_empty_region():
    df = _make_controls("OTHER")
    ref = compute_reference_ranges(df, "MISSING_REGION")
    assert len(ref) == 0


def test_flag_result():
    assert flag_result(50.0) == "NORMAL"
    assert flag_result(10.0) == "LOW"
    assert flag_result(90.0) == "HIGH"
    assert flag_result(float("nan")) == "NA"
    # Boundaries: FLAG_LOW=20.0, FLAG_HIGH=80.0
    assert flag_result(20.0) == "NORMAL"   # 20.0 is not < 20.0
    assert flag_result(19.9) == "LOW"
    assert flag_result(80.0) == "NORMAL"   # 80.0 is not > 80.0
    assert flag_result(80.1) == "HIGH"
    # Values previously flagged at old 30% threshold should now be NORMAL
    assert flag_result(28.0) == "NORMAL"
    assert flag_result(29.9) == "NORMAL"
