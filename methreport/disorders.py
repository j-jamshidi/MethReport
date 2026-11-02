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
            "maternal uniparental disomy of chromosome 7 (~10%). The GRB10:alt-TSS-DMR (7p12) and "
            "MEST:alt-TSS-DMR (7q32) are methylated on the maternal allele and unmethylated on "
            "the paternal allele, so upd(7)mat gives two methylated maternal copies and "
            "increased methylation at these loci."
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
            "Hypomethylation at GRB10:alt-TSS-DMR and/or MEST:alt-TSS-DMR indicates loss of "
            "methylation on the normally methylated maternal allele. This is not the pattern of "
            "maternal uniparental disomy of chromosome 7 (which increases methylation at these "
            "loci) and is not an established mechanism of SRS. Loss of methylation at maternally "
            "methylated DMRs is described in multi-locus imprinting disturbance (MLID); consider "
            "the other maternally methylated DMRs. Hypomethylation of ICR1 (the commonest "
            "molecular finding in SRS) is reported separately under H19/IGF2:IG-DMR."
        ),
        interpretation_high=(
            "Hypermethylation at GRB10:alt-TSS-DMR and/or MEST:alt-TSS-DMR indicates absence of "
            "the normally unmethylated paternal allele, leaving only the methylated maternal "
            "pattern. This is the pattern expected with maternal uniparental disomy of "
            "chromosome 7 (upd(7)mat), found in approximately 5-10% of SRS. Confirmation "
            "requires microsatellite or SNP-array analysis of the proband and parents."
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
            "somatic DMR at MEG3:TSS-DMR; an isolated IG-DMR defect can occur without upd(14). The paternal allele expresses DLK1 and RTL1; "
            "the maternal allele expresses MEG3 (GTL2), RTL1as, MEG8 and other non-coding RNAs. "
            "The IG-DMR and MEG3:TSS-DMR are methylated on the paternal allele; the MEG8:Int2-DMR "
            "is a secondary DMR methylated on the maternal allele, so it changes in the opposite "
            "direction to the IG-DMR and MEG3:TSS-DMR. "
            "Methylation defects cause opposite phenotypes: "
            "gain-of-methylation at the IG-DMR and MEG3:TSS-DMR silences maternal transcripts (Kagami-Ogata); "
            "loss-of-methylation at these DMRs silences paternal DLK1 and RTL1 expression (Temple syndrome)."
        ),
        key_features=[
            "HIGH methylation at IG-DMR / MEG3:TSS-DMR → Kagami-Ogata: bell-shaped thorax, placentomegaly, coat-hanger ribs",
            "LOW methylation at IG-DMR / MEG3:TSS-DMR → Temple: growth restriction, truncal obesity, early puberty",
            "MEG8:Int2-DMR (maternally methylated) moves in the opposite direction: LOW in Kagami-Ogata, HIGH in Temple",
            "Kagami-Ogata: upd(14)pat, maternal deletion or epimutation; Temple: upd(14)mat, paternal deletion or epimutation",
        ],
        interpretation_low=(
            "Hypomethylation at the DLK1/MEG3:IG-DMR and/or MEG3:TSS-DMR (normally methylated on the "
            "paternal allele) indicates loss of the paternal methylation pattern, causing reduced "
            "DLK1 and RTL1 expression. This is the hallmark finding in Temple syndrome "
            "(upd(14)mat, paternal deletion or epimutation). Hypomethylation at the MEG8:Int2-DMR "
            "(normally methylated on the maternal allele) is the reverse: it is the expected "
            "finding in Kagami-Ogata syndrome (upd(14)pat, maternal deletion or epimutation) and "
            "should be read together with the IG-DMR and MEG3:TSS-DMR."
        ),
        interpretation_high=(
            "Hypermethylation at the DLK1/MEG3:IG-DMR and/or MEG3:TSS-DMR (normally methylated on "
            "the paternal allele) indicates silencing of the maternally expressed MEG3 transcript, "
            "consistent with Kagami-Ogata syndrome (upd(14)pat, maternal deletion or epimutation). "
            "Hypermethylation at the MEG8:Int2-DMR (normally methylated on the maternal allele) is "
            "the reverse: it is the expected finding in Temple syndrome (upd(14)mat, paternal "
            "deletion or epimutation) and should be read together with the IG-DMR and "
            "MEG3:TSS-DMR."
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
            "methylated on the maternal allele and unmethylated on the paternal allele. "
            "Prader-Willi syndrome is caused by loss of the paternally expressed genes (paternal "
            "15q11-q13 deletion, maternal uniparental disomy 15 or an imprinting defect), which "
            "leaves only the methylated maternal pattern (hypermethylation). Angelman syndrome is "
            "caused by loss of maternal UBE3A expression (maternal 15q11-q13 deletion, paternal "
            "uniparental disomy 15 or an imprinting defect), which leaves only the unmethylated "
            "paternal pattern (hypomethylation). Pathogenic UBE3A variants do not change methylation."
        ),
        key_features=[
            "HIGH methylation → Prader-Willi: neonatal hypotonia, hyperphagia, obesity, intellectual disability, hypogonadism",
            "LOW methylation → Angelman: absent/minimal speech, seizures, ataxia, happy demeanour, intellectual disability",
        ],
        interpretation_low=(
            "Hypomethylation at SNURF:TSS-DMR and/or MAGEL2:TSS-DMR indicates absence of the "
            "normally methylated maternal allele, leaving only the unmethylated paternal pattern. "
            "This is the methylation finding in Angelman syndrome caused by a maternal 15q11-q13 "
            "deletion, paternal uniparental disomy 15 (upd(15)pat) or an imprinting centre defect. "
            "Methylation alone does not distinguish these mechanisms; dosage and uniparental "
            "disomy analysis are needed. Normal methylation does not exclude Angelman syndrome: "
            "UBE3A pathogenic variants and cases without an identified cause have normal "
            "methylation (methylation testing detects approximately 80% of Angelman syndrome)."
        ),
        interpretation_high=(
            "Hypermethylation at SNURF:TSS-DMR and/or MAGEL2:TSS-DMR indicates absence of the "
            "normally unmethylated paternal allele, leaving only the methylated maternal pattern. "
            "This is the methylation finding in Prader-Willi syndrome caused by a paternal "
            "15q11-q13 deletion, maternal uniparental disomy 15 (upd(15)mat) or an imprinting "
            "centre defect, and is detected in more than 99% of cases. Methylation alone does "
            "not distinguish these mechanisms; dosage and uniparental disomy analysis are needed."
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
            "patterns: GNAS exon A/B, XL and the antisense transcript (AS1) are paternally expressed; "
            "NESP55 is maternally expressed. The GNAS A/B:TSS-DMR, GNASXL:Ex1-DMR and GNAS-AS1:TSS-DMR "
            "are methylated on the maternal allele; the GNAS-NESP:TSS-DMR is methylated on the "
            "paternal allele. PHP type 1B is characterised by loss of methylation at GNAS "
            "A/B:TSS-DMR, either alone (familial disease with a maternal STX16 deletion) or "
            "together with loss of methylation at GNASXL:Ex1-DMR and GNAS-AS1:TSS-DMR and gain of "
            "methylation at GNAS-NESP:TSS-DMR (sporadic disease, paternal uniparental disomy of 20q, "
            "or maternal NESP/AS deletions)."
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
            "Hypomethylation at GNAS A/B:TSS-DMR, GNASXL:Ex1-DMR and/or GNAS-AS1:TSS-DMR "
            "(all normally methylated on the maternal allele) indicates loss of the maternal "
            "methylation imprint, consistent with PHP type 1B. Isolated loss at GNAS A/B:TSS-DMR "
            "is typical of familial disease with a maternal STX16 deletion; loss at two or more of "
            "these DMRs is typical of sporadic PHP1B, maternal NESP/AS deletions and paternal "
            "uniparental disomy of 20q. Normal methylation does not exclude PHP: GNAS coding "
            "variants (PHP type 1A) do not change methylation. Hypomethylation at "
            "GNAS-NESP:TSS-DMR (normally methylated on the paternal allele) has no established "
            "interpretation here; clinical correlation advised."
        ),
        interpretation_high=(
            "Hypermethylation at GNAS-NESP:TSS-DMR (normally methylated on the paternal allele) "
            "indicates gain of methylation on the maternal allele. Together with hypomethylation "
            "at the maternally methylated DMRs this is the pattern of sporadic PHP type 1B, "
            "including PHP1B due to paternal uniparental disomy of 20q. Hypermethylation at "
            "GNAS A/B:TSS-DMR, GNASXL:Ex1-DMR or GNAS-AS1:TSS-DMR (normally methylated on the "
            "maternal allele) is not an established pattern of PHP1B and has no established "
            "interpretation here. Clinical correlation with PTH and calcium levels is essential."
        ),
        inheritance="Complex — tissue-specific maternal/paternal imprinting",
        gene_locus="GNAS, NESP55, XLαs, GNASAS — chr20q13",
    ),
}


def get_disorder_info(disorder_key: str) -> DisorderInfo | None:
    """Return disorder info for a disorder key, or None if not found."""
    return DISORDER_INFO.get(disorder_key)
