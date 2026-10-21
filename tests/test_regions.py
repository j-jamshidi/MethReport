"""Basic tests for region definitions."""

import pytest
from methreport.regions import get_regions, get_region, REGIONS


def test_both_genomes_have_15_regions():
    for genome in ("t2t", "hg38"):
        regions = get_regions(genome)
        assert len(regions) == 15, f"{genome} should have 15 regions"


def test_region_names_consistent_across_genomes():
    t2t_names = {r.name for r in get_regions("t2t")}
    hg38_names = {r.name for r in get_regions("hg38")}
    assert t2t_names == hg38_names


def test_coordinates_differ_between_genomes():
    t2t = {r.name: (r.chrom, r.start, r.end) for r in get_regions("t2t")}
    hg38 = {r.name: (r.chrom, r.start, r.end) for r in get_regions("hg38")}
    different = sum(1 for name in t2t if t2t[name] != hg38[name])
    assert different > 0, "T2T and hg38 coordinates should differ"


def test_get_region_raises_for_unknown():
    with pytest.raises(KeyError):
        get_region("t2t", "NONEXISTENT")


def test_region_fields():
    r = get_region("t2t", "H19_IGF2_IG_DMR")
    assert r.chrom == "chr11"
    assert r.start < r.end
    assert r.disorder == "BWS/SRS"


def test_invalid_genome():
    from methreport.regions import get_regions
    with pytest.raises(ValueError):
        get_regions("hg19")


@pytest.mark.parametrize("genome,start,end", [
    ("hg38", 100809335, 100811397),
    ("t2t", 95043978, 95046609),
])
def test_ig_dmr_present_in_both_builds(genome, start, end):
    r = get_region(genome, "DLK1_MEG3_IG_DMR")
    assert (r.chrom, r.start, r.end) == ("chr14", start, end)
    assert r.disorder == "Kagami-Ogata/Temple"
    assert r.expected_methylation == "paternal"
    assert not r.unreliable


@pytest.mark.parametrize("genome", ["t2t", "hg38"])
def test_regions_are_disjoint(genome):
    """Controls are assigned to a region by coordinate, so overlap is an error."""
    regions = sorted(get_regions(genome), key=lambda r: (r.chrom, r.start))
    for a, b in zip(regions, regions[1:]):
        if a.chrom == b.chrom:
            assert a.end < b.start, f"{a.name} overlaps {b.name}"
