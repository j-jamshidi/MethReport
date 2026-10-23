"""Tests for building an in-house control cohort from per-sample output."""

import pandas as pd
import pytest

from methreport.controls import compute_reference_ranges, normalise_controls
from methreport.controls_build import build_controls, provenance
from methreport.regions import get_regions

REGION = get_regions("hg38")[0]


def _tsv(tmp_path, sample, pcts, haplotype="unphased", region=None, n_total=30):
    rows = [
        (region or REGION.label, haplotype, REGION.chrom, REGION.start + i, n_total,
         int(round(n_total * p / 100)), p)
        for i, p in enumerate(pcts)
    ]
    df = pd.DataFrame(rows, columns=["region", "haplotype", "chrom", "position",
                                     "n_total", "n_methyl", "methylation_pct"])
    path = tmp_path / f"{sample}_cpg_methylation.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_pools_unphased_rows_into_long_schema(tmp_path):
    a = _tsv(tmp_path, "ctrl_01", [50, 48])
    b = _tsv(tmp_path, "ctrl_02", [52, 50])
    out = build_controls([a, b], "hg38")
    assert list(out.columns) == ["region_name", "position", "sample_id", "methylation_pct"]
    assert set(out["sample_id"]) == {"ctrl_01", "ctrl_02"}
    assert len(out) == 4
    assert (out["region_name"] == REGION.name).all()


def test_phased_rows_and_thin_coverage_are_excluded(tmp_path):
    p = tmp_path / "ctrl_01_cpg_methylation.tsv"
    keep = _tsv(tmp_path, "keep", [50])
    df = pd.read_csv(keep, sep="\t")
    hp = df.assign(haplotype="hp1")
    thin = df.assign(position=df["position"] + 1, n_total=2)
    pd.concat([df, hp, thin]).to_csv(p, sep="\t", index=False)
    out = build_controls([p], "hg38", min_cov=5)
    assert len(out) == 1


def test_yields_mean_sd_and_n_per_position(tmp_path):
    paths = [_tsv(tmp_path, f"ctrl_0{i}", [50 + i]) for i in range(1, 5)]
    out = build_controls(paths, "hg38")
    agg = normalise_controls(out)
    assert int(agg["n_controls"].iloc[0]) == 4
    ranges = compute_reference_ranges(out, REGION.name)
    assert int(ranges["n_controls"].iloc[0]) == 4


def test_accession_shaped_sample_id_is_refused(tmp_path):
    path = _tsv(tmp_path, "26R279664X", [50])
    with pytest.raises(ValueError, match="accession"):
        build_controls([path], "hg38")


def test_region_from_another_build_is_refused(tmp_path):
    path = _tsv(tmp_path, "ctrl_01", [50], region="NOT A REGION")
    with pytest.raises(ValueError, match="not defined"):
        build_controls([path], "hg38")


def test_provenance_records_n_per_position(tmp_path):
    paths = [_tsv(tmp_path, f"ctrl_0{i}", [50, 51]) for i in range(1, 4)]
    out = build_controls(paths, "hg38")
    rec = provenance(out, "hg38", 5, "5mC", "model-x", "v0")
    assert rec["n_samples"] == 3
    assert rec["per_region"][REGION.name]["n_min"] == 3
