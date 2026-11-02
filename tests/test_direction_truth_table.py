"""Check the clinical wording and parental-allele labels against the truth table.

docs/direction_truth_table.tsv records, for every DMR and direction (LOW/HIGH),
which parental allele is methylated and which condition that direction points
to, with the supporting references. These tests stop the prose in
methreport/disorders.py and the labels in methreport/regions.py drifting away
from that table (an inverted Prader-Willi / Angelman direction went unnoticed
in 0.3.0 because nothing compared them).
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import pytest

from methreport.disorders import DISORDER_INFO, get_disorder_info
from methreport.regions import DISORDER_GROUPS, REGIONS, get_regions

TABLE = Path(__file__).resolve().parents[1] / "docs" / "direction_truth_table.tsv"

NO_CONDITION = "none established"
COLUMNS = [
    "disorder_key",
    "region_name",
    "methylated_parental_allele",
    "direction",
    "condition",
    "mechanisms",
    "source",
    "confidence",
]

# Regular expressions that identify each condition in the interpretation text.
# Abbreviations are matched case-sensitively and on word boundaries so that
# "AS" in "GNAS-AS1" is never read as Angelman syndrome.
CONDITION_PATTERNS = {
    "Prader-Willi syndrome": re.compile(r"Prader-Willi|\bPWS\b"),
    "Angelman syndrome": re.compile(r"Angelman"),
    "Silver-Russell syndrome": re.compile(r"Silver-Russell|\bSRS\b"),
    "Beckwith-Wiedemann syndrome": re.compile(r"Beckwith-Wiedemann|\bBWS\b"),
    "Kagami-Ogata syndrome": re.compile(r"Kagami-Ogata"),
    "Temple syndrome": re.compile(r"Temple"),
    "TNDM": re.compile(r"TNDM|neonatal diabetes"),
    "MLID": re.compile(r"MLID|[Mm]ulti-locus imprinting"),
    "PHP type 1B": re.compile(r"PHP\s*type\s*1[bB]|PHP1B|pseudohypoparathyroidism type 1[bB]"),
}


def _load_table() -> list[dict[str, str]]:
    with TABLE.open(encoding="utf8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        assert reader.fieldnames == COLUMNS
        return list(reader)


ROWS = _load_table()


def _text(info, direction: str) -> str:
    return info.interpretation_low if direction == "LOW" else info.interpretation_high


def _named(condition: str, text: str) -> bool:
    return bool(CONDITION_PATTERNS[condition].search(text))


def _conditions(disorder: str, direction: str | None = None) -> set[str]:
    return {
        r["condition"]
        for r in ROWS
        if r["disorder_key"] == disorder
        and r["condition"] != NO_CONDITION
        and (direction is None or r["direction"] == direction)
    }


# ---------------------------------------------------------------------------
# Table integrity
# ---------------------------------------------------------------------------
def test_table_covers_every_region_and_direction_once():
    seen = defaultdict(list)
    for r in ROWS:
        seen[(r["region_name"], r["direction"])].append(r)
    names = {r.name for r in get_regions("hg38")}
    expected = {(n, d) for n in names for d in ("LOW", "HIGH")}
    assert set(seen) == expected
    assert all(len(v) == 1 for v in seen.values())


def test_table_values_are_valid():
    for r in ROWS:
        assert r["direction"] in {"LOW", "HIGH"}
        assert r["methylated_parental_allele"] in {"maternal", "paternal"}
        assert r["confidence"] in {"confirmed", "probable", "unresolved"}
        assert r["source"].strip() and r["mechanisms"].strip()
        assert r["condition"] == NO_CONDITION or r["condition"] in CONDITION_PATTERNS
        for pmid in re.findall(r"PMID:(\S+?)(?=[;\s)]|$)", r["source"]):
            assert pmid.isdigit(), r


def test_allele_is_consistent_within_each_region():
    alleles = defaultdict(set)
    for r in ROWS:
        alleles[r["region_name"]].add(r["methylated_parental_allele"])
    assert all(len(v) == 1 for v in alleles.values())


def test_table_disorder_keys_match_regions_and_disorder_info():
    by_name = {r.name: r.disorder for r in get_regions("hg38")}
    for r in ROWS:
        assert by_name[r["region_name"]] == r["disorder_key"]
        assert r["disorder_key"] in DISORDER_INFO
    assert set(DISORDER_INFO) == {r["disorder_key"] for r in ROWS}


def test_disorder_groups_are_consistent_with_regions():
    names = {r.name for r in get_regions("hg38")}
    for members in DISORDER_GROUPS.values():
        assert set(members) <= names


# ---------------------------------------------------------------------------
# (a) interpretation text names the right condition for each direction
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("disorder", sorted(DISORDER_INFO))
@pytest.mark.parametrize("direction", ["LOW", "HIGH"])
def test_interpretation_names_the_condition_the_table_assigns(disorder, direction):
    info = get_disorder_info(disorder)
    text = _text(info, direction)
    required = _conditions(disorder, direction)
    for condition in required:
        assert _named(condition, text), (
            f"{disorder} {direction} interpretation does not name {condition}: {text!r}"
        )
    # Where a disorder entry covers two conditions, a direction must not name
    # a condition that the table assigns only to the opposite direction.
    all_conditions = _conditions(disorder)
    if len(all_conditions) >= 2:
        for condition in all_conditions - required:
            assert not _named(condition, text), (
                f"{disorder} {direction} interpretation names {condition}, which the table "
                f"assigns to the opposite direction: {text!r}"
            )


# ---------------------------------------------------------------------------
# (b) region labels equal the table
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("genome", sorted(REGIONS))
def test_expected_methylation_matches_table(genome):
    table = {r["region_name"]: r["methylated_parental_allele"] for r in ROWS}
    for region in get_regions(genome):
        assert region.expected_methylation == table[region.name], region.name


# ---------------------------------------------------------------------------
# (c) Prader-Willi / Angelman regression (0.3.0 had the direction inverted)
# ---------------------------------------------------------------------------
PWAS = "Prader-Willi/Angelman"


@pytest.mark.parametrize("region", ["SNURF_TSS_DMR", "MAGEL2_TSS_DMR"])
def test_pwas_regions_are_maternally_methylated(region):
    for genome in REGIONS:
        r = next(x for x in get_regions(genome) if x.name == region)
        assert r.disorder == PWAS
        assert r.expected_methylation == "maternal"


@pytest.mark.parametrize("region", ["SNURF_TSS_DMR", "MAGEL2_TSS_DMR"])
def test_pwas_table_rows(region):
    rows = {r["direction"]: r for r in ROWS if r["region_name"] == region}
    assert rows["HIGH"]["condition"] == "Prader-Willi syndrome"
    assert rows["LOW"]["condition"] == "Angelman syndrome"


def test_pwas_high_is_prader_willi_and_low_is_angelman():
    info = get_disorder_info(PWAS)
    high, low = info.interpretation_high, info.interpretation_low
    assert "Prader-Willi" in high and "Angelman" not in high
    assert "Angelman" in low and "Prader-Willi" not in low and "PWS" not in low
    # Mechanisms sit on the correct side of the pair.
    assert "upd(15)mat" in high and "upd(15)pat" not in high
    assert "upd(15)pat" in low and "upd(15)mat" not in low
    # Normal methylation does not exclude Angelman syndrome (UBE3A variants).
    assert "UBE3A" in low and "UBE3A" not in high


def test_pwas_key_features_and_description_match_direction():
    info = get_disorder_info(PWAS)
    by_direction = {f.split()[0]: f for f in info.key_features}
    assert "Prader-Willi" in by_direction["HIGH"] and "Angelman" not in by_direction["HIGH"]
    assert "Angelman" in by_direction["LOW"] and "Prader-Willi" not in by_direction["LOW"]
    desc = info.description
    # Prader-Willi is introduced with the paternal loss and hypermethylation,
    # Angelman with the maternal loss and hypomethylation.
    pw_at = desc.index("Prader-Willi syndrome is caused")
    as_at = desc.index("Angelman syndrome is caused")
    pw, asyn = desc[pw_at:as_at], desc[as_at:]
    assert "hypermethylation" in pw and "hypomethylation" not in pw
    assert "hypomethylation" in asyn and "hypermethylation" not in asyn


# ---------------------------------------------------------------------------
# Other corrected statements
# ---------------------------------------------------------------------------
def test_srs_upd7mat_is_hypermethylation_at_grb10_and_mest():
    info = get_disorder_info("SRS")
    assert "upd(7)mat" in info.interpretation_high
    assert "upd(7)mat" not in info.interpretation_low
    for genome in REGIONS:
        for r in get_regions(genome):
            if r.name in {"GRB10_alt_TSS_DMR", "MEST_alt_TSS_DMR"}:
                assert r.expected_methylation == "maternal"


def test_meg8_int2_dmr_moves_opposite_to_ig_and_meg3_dmrs():
    info = get_disorder_info("Kagami-Ogata/Temple")
    for genome in REGIONS:
        regs = {r.name: r for r in get_regions(genome)}
        assert regs["MEG8_Int2_DMR"].expected_methylation == "maternal"
        assert regs["DLK1_MEG3_IG_DMR"].expected_methylation == "paternal"
        assert regs["MEG3_TSS_DMR"].expected_methylation == "paternal"

    def sentence(text: str, key: str) -> str:
        return next(s for s in re.split(r"(?<=\.)\s+", text) if key in s)

    low_meg8 = sentence(info.interpretation_low, "MEG8:Int2-DMR")
    high_meg8 = sentence(info.interpretation_high, "MEG8:Int2-DMR")
    assert "Kagami-Ogata" in low_meg8 and "Temple" not in low_meg8
    assert "Temple" in high_meg8 and "Kagami-Ogata" not in high_meg8
    # The IG-DMR / MEG3:TSS-DMR sentence keeps the original direction.
    assert "Temple" in info.interpretation_low.split("Hypomethylation at the MEG8")[0]
    assert "Kagami-Ogata" in info.interpretation_high.split("Hypermethylation at the MEG8")[0]
    assert "MEG8" in info.description and "maternal allele" in info.description


def test_gnas_maternal_dmrs_and_upd20_direction():
    for genome in REGIONS:
        regs = {r.name: r.expected_methylation for r in get_regions(genome)}
        assert regs["GNAS_NESP_TSS_DMR"] == "paternal"
        assert regs["GNAS_AB_TSS_DMR"] == "maternal"
        assert regs["GNAS_AS1_TSS_DMR"] == "maternal"
        assert regs["GNASXL_Ex1_DMR"] == "maternal"
    info = get_disorder_info("PHP")
    text = " ".join([info.description, info.interpretation_low, info.interpretation_high])
    # Paternal (not maternal) UPD of 20q causes the PHP1B methylation pattern.
    assert "upd(20)mat" not in text
    assert "uniparental disomy of 20q" in info.interpretation_low
    assert "supports maternal-allele loss" not in text
    # Gain of NESP55 methylation, not loss, accompanies PHP1B.
    assert "gain of methylation" in info.interpretation_high
    assert "GNAS-NESP:TSS-DMR" in info.interpretation_high


def test_public_dataclass_fields_unchanged():
    info = get_disorder_info(PWAS)
    for field in (
        "full_name", "subtitle", "omim", "description", "key_features",
        "interpretation_low", "interpretation_high", "inheritance", "gene_locus",
    ):
        assert hasattr(info, field)
