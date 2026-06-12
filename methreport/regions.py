"""DMR region definitions for T2T-CHM13v2.0 and hg38 reference genomes.

Coordinates derived from NanoImprint (carolinehey/NanoImprint) covering the
14 differentially methylated regions associated with known imprinting disorders.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class DMRegion:
    name: str
    label: str          # display name in reports
    chrom: str
    start: int
    end: int
    disorder: str
    expected_methylation: str   # "maternal" | "paternal" | "biallelic"
    # For imprinted DMRs: which allele is methylated in normal controls
    # maternal = methylated on maternal allele (~50% in unphased)
    # paternal = methylated on paternal allele (~50% in unphased)


# ---------------------------------------------------------------------------
# T2T-CHM13v2.0 regions
# ---------------------------------------------------------------------------
T2T_REGIONS: list[DMRegion] = [
    DMRegion(
        name="PLAGL1_alt_TSS_DMR",
        label="PLAGL1:alt-TSS-DMR",
        chrom="chr6",
        start=145200200,
        end=145201432,
        disorder="TNDM",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GRB10_alt_TSS_DMR",
        label="GRB10:alt-TSS-DMR",
        chrom="chr7",
        start=50943355,
        end=50944182,
        disorder="SRS",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="MEST_alt_TSS_DMR",
        label="MEST:alt-TSS-DMR",
        chrom="chr7",
        start=131804587,
        end=131807000,
        disorder="SRS",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="H19_IGF2_IG_DMR",
        label="H19/IGF2:IG-DMR",
        chrom="chr11",
        start=2085438,
        end=2090906,
        disorder="BWS/SRS",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="KCNQ1OT1_TSS_DMR",
        label="KCNQ1OT1:TSS-DMR",
        chrom="chr11",
        start=2788106,
        end=2790240,
        disorder="BWS",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="MEG3_TSS_DMR",
        label="MEG3:TSS-DMR",
        chrom="chr14",
        start=95059075,
        end=95062529,
        disorder="Kagami-Ogata/Temple",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="MEG8_Int2_DMR",
        label="MEG8:Int2-DMR",
        chrom="chr14",
        start=95140000,
        end=95140384,
        disorder="Kagami-Ogata/Temple",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="MAGEL2_TSS_DMR",
        label="MAGEL2:TSS-DMR",
        chrom="chr15",
        start=21381784,
        end=21383056,
        disorder="Prader-Willi/Angelman",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="SNURF_TSS_DMR",
        label="SNURF:TSS-DMR",
        chrom="chr15",
        start=22691523,
        end=22693493,
        disorder="Prader-Willi/Angelman",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="PEG3_TSS_DMR",
        label="PEG3:TSS-DMR",
        chrom="chr19",
        start=59932193,
        end=59936667,
        disorder="MLID",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GNAS_NESP_TSS_DMR",
        label="GNAS-NESP:TSS-DMR",
        chrom="chr20",
        start=60622124,
        end=60626697,
        disorder="PHP",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="GNAS_AS1_TSS_DMR",
        label="GNAS-AS1:TSS-DMR",
        chrom="chr20",
        start=60633713,
        end=60636097,
        disorder="PHP",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GNASXL_Ex1_DMR",
        label="GNASXL:Ex1-DMR",
        chrom="chr20",
        start=60637080,
        end=60639527,
        disorder="PHP",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GNAS_AB_TSS_DMR",
        label="GNAS A/B:TSS-DMR",
        chrom="chr20",
        start=60671509,
        end=60673270,
        disorder="PHP",
        expected_methylation="paternal",
    ),
]

# ---------------------------------------------------------------------------
# hg38 regions
# ---------------------------------------------------------------------------
HG38_REGIONS: list[DMRegion] = [
    DMRegion(
        name="PLAGL1_alt_TSS_DMR",
        label="PLAGL1:alt-TSS-DMR",
        chrom="chr6",
        start=144007487,
        end=144008719,
        disorder="TNDM",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GRB10_alt_TSS_DMR",
        label="GRB10:alt-TSS-DMR",
        chrom="chr7",
        start=50782178,
        end=50783005,
        disorder="SRS",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="MEST_alt_TSS_DMR",
        label="MEST:alt-TSS-DMR",
        chrom="chr7",
        start=130490866,
        end=130493279,
        disorder="SRS",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="H19_IGF2_IG_DMR",
        label="H19/IGF2:IG-DMR",
        chrom="chr11",
        start=1997761,
        end=2003229,
        disorder="BWS/SRS",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="KCNQ1OT1_TSS_DMR",
        label="KCNQ1OT1:TSS-DMR",
        chrom="chr11",
        start=2698767,
        end=2700902,
        disorder="BWS",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="MEG3_TSS_DMR",
        label="MEG3:TSS-DMR",
        chrom="chr14",
        start=100824187,
        end=100827641,
        disorder="Kagami-Ogata/Temple",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="MEG8_Int2_DMR",
        label="MEG8:Int2-DMR",
        chrom="chr14",
        start=100904568,
        end=100904952,
        disorder="Kagami-Ogata/Temple",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="MAGEL2_TSS_DMR",
        label="MAGEL2:TSS-DMR",
        chrom="chr15",
        start=23647365,
        end=23648635,
        disorder="Prader-Willi/Angelman",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="SNURF_TSS_DMR",
        label="SNURF:TSS-DMR",
        chrom="chr15",
        start=24954857,
        end=24956829,
        disorder="Prader-Willi/Angelman",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="PEG3_TSS_DMR",
        label="PEG3:TSS-DMR",
        chrom="chr19",
        start=56837125,
        end=56841599,
        disorder="MLID",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GNAS_NESP_TSS_DMR",
        label="GNAS-NESP:TSS-DMR",
        chrom="chr20",
        start=58838984,
        end=58843557,
        disorder="PHP",
        expected_methylation="paternal",
    ),
    DMRegion(
        name="GNAS_AS1_TSS_DMR",
        label="GNAS-AS1:TSS-DMR",
        chrom="chr20",
        start=58850594,
        end=58852978,
        disorder="PHP",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GNASXL_Ex1_DMR",
        label="GNASXL:Ex1-DMR",
        chrom="chr20",
        start=58853961,
        end=58856408,
        disorder="PHP",
        expected_methylation="maternal",
    ),
    DMRegion(
        name="GNAS_AB_TSS_DMR",
        label="GNAS A/B:TSS-DMR",
        chrom="chr20",
        start=58888390,
        end=58890146,
        disorder="PHP",
        expected_methylation="paternal",
    ),
]

REGIONS: dict[str, list[DMRegion]] = {
    "t2t": T2T_REGIONS,
    "hg38": HG38_REGIONS,
}

DISORDER_GROUPS: dict[str, list[str]] = {
    "TNDM": ["PLAGL1_alt_TSS_DMR"],
    "SRS": ["GRB10_alt_TSS_DMR", "MEST_alt_TSS_DMR", "H19_IGF2_IG_DMR"],
    "BWS": ["H19_IGF2_IG_DMR", "KCNQ1OT1_TSS_DMR"],
    "Kagami-Ogata/Temple": ["MEG3_TSS_DMR", "MEG8_Int2_DMR"],
    "Prader-Willi/Angelman": ["MAGEL2_TSS_DMR", "SNURF_TSS_DMR"],
    "MLID": ["PEG3_TSS_DMR"],
    "PHP": ["GNAS_NESP_TSS_DMR", "GNAS_AS1_TSS_DMR", "GNASXL_Ex1_DMR", "GNAS_AB_TSS_DMR"],
}


def get_regions(genome: str) -> list[DMRegion]:
    genome = genome.lower()
    if genome not in REGIONS:
        raise ValueError(f"Unknown genome '{genome}'. Choose from: {list(REGIONS)}")
    return REGIONS[genome]


def get_region(genome: str, name: str) -> DMRegion:
    for r in get_regions(genome):
        if r.name == name:
            return r
    raise KeyError(f"Region '{name}' not found for genome '{genome}'")
