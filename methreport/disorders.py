"""Clinical disorder information and OMIM references for imprinting disorders.

Each entry maps a disorder key (matching DMRegion.disorder) to a DisorderInfo
describing the clinical phenotype, OMIM accessions, and direction-specific
methylation interpretation.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class OmimEntry:
    accession: str      # numeric OMIM ID, no '#'
    label: str

    @property
    def url(self) -> str:
        return f"https://www.omim.org/entry/{self.accession}"

    @property
    def badge(self) -> str:
        return f'<a class="omim-badge" href="{self.url}" target="_blank" rel="noopener">OMIM #{self.accession}</a>'


@dataclass(frozen=True)
class DisorderInfo:
    full_name: str
    subtitle: str
    omim: list[OmimEntry]
    description: str
    key_features: list[str]
    # Direction-specific interpretation: what LOW / HIGH methylation means
    interpretation_low: str
    interpretation_high: str
    inheritance: str
    gene_locus: str


DISORDER_INFO: dict[str, DisorderInfo] = {
    "TNDM": DisorderInfo(
        full_name="Transient Neonatal Diabetes Mellitus",
        subtitle="6q24-related imprinting disorder",
        omim=[OmimEntry("601410", "TNDM, 6q24-related")],
        description=(
            "A rare form of neonatal diabetes presenting within the first weeks of life due to "
            "overexpression of the imprinted genes PLAGL1 and HYMAI on 6q24. Diabetes typically "
            "resolves within 12–18 months but recurs in approximately 50% of patients during "
            "adolescence or early adulthood. Caused by paternal uniparental disomy (upd(6)pat), "
            "paternal duplication of 6q24, or loss of methylation at the maternal PLAGL1:alt-TSS-DMR."
        ),
        key_features=[
            "Neonatal hyperglycaemia (onset &lt;6 months)",
            "Intrauterine growth restriction",
            "Macroglossia",
            "Umbilical hernia",
            "Transient course (resolves by 18 months)",
            "Risk of relapse in adolescence (~50%)",
        ],
        interpretation_low=(
            "Hypomethylation at PLAGL1:alt-TSS-DMR indicates loss of the normally-methylated "
            "maternal allele, leading to biallelic PLAGL1/HYMAI expression. This is consistent "
            "with 6q24-related TNDM via epimutation."
        ),
        interpretation_high=(
            "Hypermethylation at PLAGL1:alt-TSS-DMR is not a recognised cause of TNDM and may "
            "represent a technical artefact or benign variant. Clinical correlation advised."
        ),
        inheritance="Paternal overexpression (loss of maternal imprint)",
        gene_locus="PLAGL1/HYMAI — chr6q24",
    ),

    "SRS": DisorderInfo(
        full_name="Silver-Russell Syndrome",
        subtitle="11p15 / chromosome 7 imprinting disorder",
        omim=[OmimEntry("180860", "Silver-Russell Syndrome")],
        description=(
            "A growth restriction syndrome characterised by severe pre- and post-natal growth "
            "retardation, relative macrocephaly at birth, and variable body asymmetry. Two main "
            "molecular causes: hypomethylation of H19/IGF2:IG-DMR (ICR1) on 11p15 (~40%) and "
            "maternal uniparental disomy of chromosome 7 (~10%), reflected in GRB10 and MEST "
            "methylation changes."
        ),
        key_features=[
            "Pre- and postnatal growth restriction (height &lt;-2 SD)",
            "Relative macrocephaly at birth",
            "Prominent forehead",
            "Body asymmetry (hemihypotrophy)",
            "Feeding difficulties and failure to thrive",
            "Clinodactyly of 5th finger",
        ],
        interpretation_low=(
            "Hypomethylation suggests loss of paternal imprint, reducing IGF2 expression. "
            "At H19/IGF2:IG-DMR this is the key finding in ~40% of SRS. "
            "At GRB10/MEST, reduced methylation may reflect upd(7)mat."
        ),
        interpretation_high=(
            "Hypermethylation at these loci is not a recognised mechanism for SRS. "
            "ICR1 hypermethylation is instead associated with Beckwith-Wiedemann syndrome."
        ),
        inheritance="Paternal imprint (biallelic IGF2 loss → undergrowth)",
        gene_locus="H19/IGF2 — chr11p15.5 · GRB10, MEST — chr7",
    ),

    "BWS": DisorderInfo(
        full_name="Beckwith-Wiedemann Syndrome",
        subtitle="11p15 imprinting disorder",
        omim=[OmimEntry("130650", "Beckwith-Wiedemann Syndrome")],
        description=(
            "An overgrowth and tumour predisposition syndrome caused by dysregulation of "
            "imprinted genes on chromosome 11p15.5. The two imprinting centres involved are "
            "ICR1 (H19/IGF2:IG-DMR) and ICR2 (KCNQ1OT1:TSS-DMR). Loss of methylation at ICR2 "
            "is the most common molecular cause (~50%), while gain of methylation at ICR1 "
            "accounts for ~5–10% of cases."
        ),
        key_features=[
            "Macrosomia and overgrowth",
            "Macroglossia",
            "Anterior abdominal wall defects (omphalocele, umbilical hernia)",
            "Neonatal hypoglycaemia",
            "Ear creases / pits",
            "Hemihyperplasia",
            "Embryonal tumour predisposition (Wilms tumour, hepatoblastoma) — 7–10%",
        ],
        interpretation_low=(
            "Hypomethylation at KCNQ1OT1:TSS-DMR (ICR2) is the most common molecular cause of BWS "
            "(∼50% of cases), resulting in paternal expression of KCNQ1OT1 from both alleles "
            "and silencing of maternally expressed growth-suppressor genes in the domain."
        ),
        interpretation_high=(
            "Hypermethylation at KCNQ1OT1:TSS-DMR is not a recognised BWS mechanism. "
            "Consider H19/IGF2:IG-DMR (ICR1) results for a complete 11p15 assessment."
        ),
        inheritance="Maternal imprint (loss of ICR2 maternal methylation)",
        gene_locus="KCNQ1OT1, CDKN1C, KCNQ1 — chr11p15.5",
    ),

    "BWS/SRS": DisorderInfo(
        full_name="Beckwith-Wiedemann Syndrome / Silver-Russell Syndrome",
        subtitle="H19/IGF2 IG-DMR (ICR1) — 11p15.5",
        omim=[
            OmimEntry("130650", "Beckwith-Wiedemann Syndrome"),
            OmimEntry("180860", "Silver-Russell Syndrome"),
        ],
        description=(
            "The H19/IGF2 Imprinting Control Region 1 (ICR1) on 11p15.5 is normally methylated "
            "on the paternal allele, enabling monoallelic IGF2 expression (paternal) and H19 "
            "expression (maternal). Methylation defects at this locus cause opposite growth "
            "phenotypes: hypomethylation → biallelic H19 / reduced IGF2 → Silver-Russell syndrome; "
            "hypermethylation → biallelic IGF2 / silenced H19 → Beckwith-Wiedemann syndrome."
        ),
        key_features=[
            "LOW methylation → SRS: growth restriction, relative macrocephaly, feeding difficulties",
            "HIGH methylation → BWS: overgrowth, macroglossia, tumour predisposition",
        ],
        interpretation_low=(
            "Hypomethylation at H19/IGF2:IG-DMR indicates loss of normal paternal methylation, "
            "resulting in biallelic H19 expression and reduced IGF2 levels. This pattern is found "
            "in approximately 40% of patients with Silver-Russell syndrome."
        ),
        interpretation_high=(
            "Hypermethylation at H19/IGF2:IG-DMR indicates gain of methylation on the normally "
            "unmethylated maternal allele, causing biallelic IGF2 overexpression. This is found "
            "in approximately 5–10% of Beckwith-Wiedemann syndrome cases."
        ),
        inheritance="Paternal imprint — ICR1",
        gene_locus="H19 / IGF2 — chr11p15.5",
    ),

    "Kagami-Ogata/Temple": DisorderInfo(
        full_name="Kagami-Ogata Syndrome / Temple Syndrome",
        subtitle="14q32 imprinting disorder (MEG3/DLK1 domain)",
        omim=[
            OmimEntry("608149", "Kagami-Ogata Syndrome"),
            OmimEntry("616222", "Temple Syndrome"),
        ],
        description=(
            "The DLK1-MEG3 imprinting domain on chromosome 14q32 is controlled by a "
            "germline-derived differentially methylated region (IG-DMR) and a "
            "somatic DMR at MEG3:TSS-DMR. The paternal allele expresses DLK1 and MEG8; "
            "the maternal allele expresses MEG3 (GTL2), RTL1as, and other non-coding RNAs. "
            "Methylation defects cause opposite phenotypes: "
            "gain-of-methylation silences maternal transcripts (Kagami-Ogata); "
            "loss-of-methylation silences paternal DLK1 expression (Temple syndrome)."
        ),
        key_features=[
            "HIGH methylation → Kagami-Ogata: bell-shaped thorax, placentomegaly, coat-hanger ribs",
            "LOW methylation → Temple: growth restriction, truncal obesity, early puberty",
            "Both associated with upd(14) or 14q32 epimutations",
        ],
        interpretation_low=(
            "Hypomethylation at MEG3:TSS-DMR and/or MEG8:Int2-DMR suggests loss of the normally "
            "methylated paternal allele, causing reduced DLK1 expression. This is the hallmark "
            "finding in Temple syndrome (upd(14)mat or epimutation)."
        ),
        interpretation_high=(
            "Hypermethylation at MEG3:TSS-DMR and/or MEG8:Int2-DMR indicates silencing of the "
            "maternally expressed MEG3 transcript, consistent with Kagami-Ogata syndrome "
            "(upd(14)pat or paternal duplication of 14q32)."
        ),
        inheritance="Paternal imprint — IG-DMR / MEG3:TSS-DMR",
        gene_locus="MEG3, MEG8, DLK1 — chr14q32",
    ),

    "Prader-Willi/Angelman": DisorderInfo(
        full_name="Prader-Willi Syndrome / Angelman Syndrome",
        subtitle="15q11-q13 imprinting disorder (SNURF-SNRPN domain)",
        omim=[
            OmimEntry("176270", "Prader-Willi Syndrome"),
            OmimEntry("105830", "Angelman Syndrome"),
        ],
        description=(
            "The SNURF-SNRPN imprinting centre on 15q11-q13 controls expression of a large "
            "cluster of paternally expressed genes (SNURF/SNRPN, MAGEL2, NDN, snoRNAs) and "
            "the maternally expressed UBE3A. The DMRs at SNURF:TSS and MAGEL2:TSS are "
            "methylated on the maternal allele. Loss of paternal expression (hypomethylation) "
            "causes Prader-Willi syndrome; loss of maternal UBE3A expression causes Angelman syndrome."
        ),
        key_features=[
            "LOW methylation → Prader-Willi: neonatal hypotonia, hyperphagia, obesity, intellectual disability, hypogonadism",
            "HIGH methylation → Angelman: absent/minimal speech, seizures, ataxia, happy demeanour, intellectual disability",
        ],
        interpretation_low=(
            "Hypomethylation at SNURF:TSS-DMR and/or MAGEL2:TSS-DMR indicates loss of the normally "
            "methylated maternal allele at these loci, consistent with Prader-Willi syndrome "
            "(upd(15)mat, 15q11-q13 deletion of paternal allele, or imprinting centre defect). "
            "DNA methylation analysis has >99% sensitivity for PW syndrome."
        ),
        interpretation_high=(
            "Hypermethylation at SNURF:TSS-DMR may indicate loss of the normally unmethylated "
            "paternal allele, as seen in Angelman syndrome caused by paternal uniparental disomy "
            "(upd(15)pat) or imprinting centre defects. Note: UBE3A mutation and 15q11-q13 "
            "deletion AS cases may show normal or only mildly abnormal methylation."
        ),
        inheritance="Maternal imprint — SNURF/SNRPN imprinting centre",
        gene_locus="SNURF, SNRPN, MAGEL2, NDN, UBE3A — chr15q11-q13",
    ),

    "MLID": DisorderInfo(
        full_name="Multi-Locus Imprinting Disturbance",
        subtitle="Imprinting defects at multiple loci",
        omim=[OmimEntry("617107", "MLID")],
        description=(
            "A condition in which loss of methylation (or rarely, gain) is present at multiple "
            "imprinting loci simultaneously. MLID most commonly presents with a Beckwith-Wiedemann "
            "or TNDM-like phenotype, but with additional DMR abnormalities. PEG3:TSS-DMR "
            "hypomethylation is frequently found in MLID. Pathogenic variants in imprinting "
            "maintenance factors (ZFP57, NLRP5, KCNQ1OT1, DPPA3) can cause MLID."
        ),
        key_features=[
            "Hypomethylation at PEG3:TSS-DMR",
            "Often co-occurs with BWS or TNDM methylation defects",
            "Variable clinical phenotype",
            "May be caused by ZFP57 or NLRP5 pathogenic variants",
            "Maternal-effect mutations (NLRP genes) may cause recurrence",
        ],
        interpretation_low=(
            "Hypomethylation at PEG3:TSS-DMR on 19q13 is a recognised marker of MLID. "
            "Consider testing other imprinting loci (particularly PLAGL1, KCNQ1OT1) and "
            "sequencing of ZFP57 and NLRP5 if multiple loci are affected."
        ),
        interpretation_high=(
            "Hypermethylation at PEG3:TSS-DMR is not a well-characterised MLID finding. "
            "Clinical correlation and repeat testing advised."
        ),
        inheritance="Variable (trans-acting maternal-effect mutations or sporadic)",
        gene_locus="PEG3 — chr19q13",
    ),

    "PHP": DisorderInfo(
        full_name="Pseudohypoparathyroidism — GNAS-related Imprinting Disorder",
        subtitle="20q13 GNAS locus imprinting disorder",
        omim=[
            OmimEntry("103580", "Pseudohypoparathyroidism type Ia / Albright hereditary osteodystrophy"),
            OmimEntry("603233", "Pseudohypoparathyroidism type Ib"),
        ],
        description=(
            "A group of disorders characterised by end-organ resistance to PTH (and other "
            "hormones using Gsα signalling) due to dysregulation of the complex GNAS imprinting "
            "cluster on 20q13. The locus encodes multiple transcripts with distinct imprinting "
            "patterns: GNAS exon A/B and XL are paternally expressed; NESP55 is maternally "
            "expressed; AS1 is expressed from both alleles with paternal predominance. "
            "PHP type 1b is most commonly caused by methylation defects at GNAS A/B:TSS-DMR "
            "and/or GNAS-AS1:TSS-DMR."
        ),
        key_features=[
            "Hypocalcaemia (PTH resistance)",
            "Hyperphosphataemia",
            "Elevated PTH",
            "TSH resistance (PHP type 1a)",
            "Albright hereditary osteodystrophy features (PHP1a): short stature, obesity, brachydactyly",
            "Ossifications (progressive osseous heteroplasia — POH)",
        ],
        interpretation_low=(
            "Hypomethylation at GNAS A/B:TSS-DMR and/or AS1/XL DMRs indicates loss of "
            "normal methylation at maternally methylated elements, consistent with PHP type 1b "
            "(AD-PHP1b or sporadic epimutation). This pattern is associated with upd(20)mat "
            "or STX16/NESP55 deletions. Reduced NESP55:TSS-DMR methylation supports "
            "maternal-allele loss."
        ),
        interpretation_high=(
            "Hypermethylation at GNAS loci is less common but may indicate gain of methylation "
            "on the paternal allele. Clinical correlation with PTH and calcium levels is essential."
        ),
        inheritance="Complex — tissue-specific maternal/paternal imprinting",
        gene_locus="GNAS, NESP55, XLαs, GNASAS — chr20q13",
    ),
}


def get_disorder_info(disorder_key: str) -> DisorderInfo | None:
    """Return disorder info for a disorder key, or None if not found."""
    return DISORDER_INFO.get(disorder_key)
