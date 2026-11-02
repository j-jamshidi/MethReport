# Review of methylation direction statements (release 0.3.1)

Scope: every clinical statement in `methreport/disorders.py` that says which
methylation direction (LOW or HIGH) points to which condition, and every
`expected_methylation` label in `methreport/regions.py`. The flagging logic,
coordinates, region names and dataclass fields are unchanged.

The machine-readable record is `docs/direction_truth_table.tsv` (one row per DMR
and direction, 30 rows). `tests/test_direction_truth_table.py` compares the
prose and the labels with that table.

Confidence terms used in the table:

* confirmed: a source I opened states the allele assignment and the direction
  (or the mechanism and its methylation result).
* probable: the allele assignment is confirmed and the direction follows from it
  by the mechanism, or the absence of a finding rests on its absence from the
  consensus statements consulted. No source states the exact result.
* unresolved: I could not establish it. The text does not make a claim.

## 1. Summary of the defect

In 0.3.0 the Prader-Willi / Angelman entry had the direction inverted, and
contradicted itself. SNURF:TSS-DMR and MAGEL2:TSS-DMR are methylated on the
maternal allele. Loss of the paternal allele (paternal deletion, upd(15)mat,
imprinting-centre defect) leaves only the methylated maternal pattern, so
hypermethylation indicates Prader-Willi syndrome. Loss of the maternal allele
(maternal deletion, upd(15)pat, imprinting-centre defect) leaves only the
unmethylated paternal pattern, so hypomethylation indicates Angelman syndrome.
The numeric LOW/HIGH flags were correct; the explanatory text attached to them
was not.

The check of the remaining entries found further errors: SRS (GRB10/MEST),
Kagami-Ogata / Temple (MEG8:Int2-DMR), PHP (upd(20) and NESP55), and three
`expected_methylation` labels. Each is listed below.

## 2. Statements changed

### 2.1 Prader-Willi / Angelman (`disorders.py`)

| Field | Old (0.3.0) | New (0.3.1) | Sources |
|---|---|---|---|
| description | "Loss of paternal expression (hypomethylation) causes Prader-Willi syndrome; loss of maternal UBE3A expression causes Angelman syndrome." | Prader-Willi syndrome leaves only the methylated maternal pattern (hypermethylation); Angelman syndrome leaves only the unmethylated paternal pattern (hypomethylation); UBE3A variants do not change methylation. | PMID:31235867, 20459762, 26583054 |
| key_features[0] | "LOW methylation → Prader-Willi: ..." | "HIGH methylation → Prader-Willi: ..." | as above |
| key_features[1] | "HIGH methylation → Angelman: ..." | "LOW methylation → Angelman: ..." | as above |
| interpretation_low | Named Prader-Willi (upd(15)mat, paternal deletion, imprinting-centre defect). Also stated ">99% sensitivity for PW syndrome". Described the same mechanisms as both "loss of the methylated maternal allele" and as upd(15)mat, which is internally contradictory. | Names Angelman syndrome: maternal deletion, upd(15)pat, imprinting-centre defect. States that normal methylation does not exclude Angelman syndrome (UBE3A variants, cases without an identified cause; methylation detects about 80%). | PMID:31235867 |
| interpretation_high | Named Angelman (upd(15)pat, imprinting-centre defects) and said 15q11-q13 deletion AS cases "may show normal or only mildly abnormal methylation". | Names Prader-Willi syndrome: paternal deletion, upd(15)mat, imprinting-centre defect; detected in more than 99%. The statement about deletion cases is removed (see 3.6). | PMID:31235867 |

Evidence. Beygo 2019 (PMID:31235867) states that SNURF, MAGEL2 and NDN
promoter DMRs "are methylated on the maternal chromosome", names them
SNURF:TSS-DMR, MAGEL2:TSS-DMR and NDN:TSS-DMR, and reports "no paternal,
unmethylated band" as the PWS result (paternal deletion, upd(15)mat or imprinting
defect) and "no maternal, methylated band" as the AS result (maternal deletion,
upd(15)pat or imprinting defect). Detection rates: more than 99% of PWS and
about 80% of AS. Proportions: PWS deletion 70-75%, upd(15)mat 25-30%,
imprinting defect about 1%; AS deletion 75%, upd(15)pat 1-2%, imprinting defect
3%, UBE3A variant 5-10%, no identified cause 10-15%.

### 2.2 Silver-Russell syndrome, GRB10 and MEST (`disorders.py`)

| Field | Old | New | Sources |
|---|---|---|---|
| description | "...reflected in GRB10 and MEST methylation changes." | States that GRB10:alt-TSS-DMR and MEST:alt-TSS-DMR are methylated on the maternal allele and that upd(7)mat gives increased methylation. | PMID:12700169, 9192843, 38164354 |
| interpretation_low | "At GRB10/MEST, reduced methylation may reflect upd(7)mat." | Hypomethylation at GRB10/MEST is not the pattern of upd(7)mat and is not an established SRS mechanism; loss of methylation at maternally methylated DMRs is described in MLID. ICR1 hypomethylation is referred to the H19/IGF2:IG-DMR result. | PMID:35296332, 38164354 |
| interpretation_high | "Hypermethylation at these loci is not a recognised mechanism for SRS." | Hypermethylation at GRB10/MEST is the pattern expected with upd(7)mat (5-10% of SRS); confirm by microsatellite or SNP analysis. | PMID:27585961, 38164354 |

Evidence. GRB10: the second CpG island is a DMR with maternal methylation in
human and mouse (Arnaud 2003, PMID:12700169). MEST: the CpG island is
unmethylated on the active paternal allele and methylated on the inactive
maternal allele (Riesewijk 1997, PMID:9192843); in a 183-subject SRS series,
hypermethylated MEST identified maternal UPD7 (Lin 2024, PMID:38164354). The
user's suspicion was correct. The one gap is GRB10 in upd(7)mat: no report with
measured GRB10 values was found, so that row is "probable" (section 3.1).

### 2.3 Kagami-Ogata / Temple (`disorders.py`)

| Field | Old | New | Sources |
|---|---|---|---|
| description | "The paternal allele expresses DLK1 and MEG8; the maternal allele expresses MEG3 (GTL2), RTL1as..." | Paternal allele expresses DLK1 and RTL1; maternal allele expresses MEG3, RTL1as and MEG8. States that IG-DMR and MEG3:TSS-DMR are paternally methylated and MEG8:Int2-DMR is maternally methylated and so moves in the opposite direction. | PMID:26377239, 41276848, 38715103, 24402520 |
| key_features | "Both associated with upd(14) or 14q32 epimutations" | Mechanisms assigned to the correct syndrome (KOS: upd(14)pat, maternal deletion, epimutation; TS: upd(14)mat, paternal deletion, epimutation). Added the MEG8 reversal. | PMID:26583054, 24891339, 26377239 |
| interpretation_low | Hypomethylation at IG-DMR, MEG3:TSS-DMR and/or MEG8:Int2-DMR all read as Temple syndrome. | Hypomethylation at IG-DMR / MEG3:TSS-DMR reads as Temple syndrome; hypomethylation at MEG8:Int2-DMR is the reverse (Kagami-Ogata). | as above |
| interpretation_high | Hyper at all three DMRs read as Kagami-Ogata (upd(14)pat "or paternal duplication of 14q32"). | Hypermethylation at IG-DMR / MEG3:TSS-DMR reads as Kagami-Ogata; at MEG8:Int2-DMR it reads as Temple. "Paternal duplication of 14q32" removed: no source lists it as a Kagami-Ogata mechanism; the consensus mechanisms are upd(14)pat, epimutation and maternal microdeletion. | PMID:26583054, 26377239 |

Evidence. Court 2014 (PMID:24402520) identifies "a maternally methylated DMR
within intron 2 of MEG8". Baena 2024 (PMID:38715103) calls the MEG8 and DLK1
DMRs "somatically acquired DMRs methylated on the maternal allele" and records
Temple syndrome as hypomethylation at IG-DMR and MEG3-DMR with hypermethylation
at MEG8. D'Angelo 2025 (PMID:41276848) states the MEG8:Int2-DMR gains methylation
in Temple syndrome. Urakawa 2025 (PMID:41254794) reports a Kagami-Ogata case
with paternal isodisomy 14 showing MEG3:TSS-DMR hypermethylation and MEG8:Int2-DMR
hypomethylation. In 0.3.0 the MEG8 label and text were both wrong.

### 2.4 PHP / GNAS (`disorders.py`)

| Field | Old | New | Sources |
|---|---|---|---|
| description | "AS1 is expressed from both alleles with paternal predominance"; "PHP type 1b is most commonly caused by methylation defects at GNAS A/B:TSS-DMR and/or GNAS-AS1:TSS-DMR." | AS1 is paternally expressed. Gives the allele of each DMR (A/B, XL, AS1 maternal; NESP paternal). PHP1B: loss of methylation at A/B, alone or with XL/AS1 loss and NESP gain. | PMID:29959430, 20965295 |
| interpretation_low | "This pattern is associated with upd(20)mat or STX16/NESP55 deletions. Reduced NESP55:TSS-DMR methylation supports maternal-allele loss." | Removed both sentences. PHP1B arises with paternal UPD of 20q (not maternal), isolated A/B loss with STX16 deletion, broad loss with NESP/AS deletions. NESP hypomethylation: no established interpretation. Normal methylation does not exclude PHP1A. | PMID:29959430, 26583054, 20965295 |
| interpretation_high | "Hypermethylation at GNAS loci ... may indicate gain of methylation on the paternal allele." | Hypermethylation at NESP (paternally methylated) is gain of methylation on the maternal allele and is the expected companion of PHP1B (sporadic, UPD20q pat). Hypermethylation at A/B, XL or AS1: no established interpretation. | PMID:29959430, 20965295, 39541438 |

Evidence. Mantovani 2018 (PMID:29959430): one paternally methylated DMR
(GNAS-NESP:TSS-DMR) and three maternally methylated DMRs (GNAS-AS1:TSS-DMR,
GNAS-XL:Ex1-DMR, GNAS A/B:TSS-DMR); loss of methylation at A/B in all PHP1B;
familial disease with STX16 deletion is isolated A/B loss; nearly all sporadic
cases lose two or more additional DMRs; about 8-10% of sporadic cases are paternal
UPD(20q). Bastepe 2011 (PMID:20965295): in paternal isodisomy 20 the pattern is
loss of A/B, XL and AS methylation with gain of methylation at NESP55, and
XL, A/B and antisense are paternal transcripts while NESP55 is maternal.
Eggermann 2015 (PMID:26583054) gives paternal, not maternal, UPD of 20q as the
(small) UPD subset of PHP1B, and lists maternal UPD20 as a separate
growth-failure disorder for which "other types of molecular alterations have not
yet been reported".

### 2.5 `regions.py` labels (both genomes; coordinates untouched)

| Region | Old | New | Sources |
|---|---|---|---|
| GRB10_alt_TSS_DMR | paternal | maternal | PMID:12700169 |
| GNAS_AB_TSS_DMR | paternal | maternal | PMID:29959430, 20965295 |
| MEG8_Int2_DMR | paternal | maternal | PMID:24402520, 38715103, 41276848 |

Added a comment on the dataclass field stating the meaning of the label.
`expected_methylation` is still not read by any flagging logic.

## 3. Items that could not be confirmed (clinician decision needed)

3.1 GRB10 in upd(7)mat. Maternal methylation of the GRB10 DMR is confirmed;
that upd(7)mat gives hypermethylation is the expected consequence (two
methylated maternal copies) and matches the confirmed MEST result, but I found
no primary report with measured GRB10 values. Marked probable; the SRS HIGH text
says "consistent with" rather than "diagnostic".

3.2 MEG8:Int2-DMR hypomethylation in Kagami-Ogata. The direction follows from
the confirmed maternal methylation and one published case (PMID:41254794), not
from a series. In addition, isolated MEG8 hypomethylation has been reported with
pathogenic PLAG1 variants in SRS-like patients with normal IG-DMR and MEG3-DMR
(PMID:41276848, 10 of 17 cases; single study). Candidate wording, not added to the
report: "Isolated hypomethylation at MEG8:Int2-DMR with normal IG-DMR and
MEG3:TSS-DMR is not specific for Kagami-Ogata syndrome." Needs a decision.

3.3 PHP patterns with no established meaning. (i) Hypomethylation at
GNAS-NESP:TSS-DMR; (ii) hypermethylation at A/B, XL or AS1. For (ii), upd(20)mat
would be expected to give it by inheritance logic, but no published GNAS
methylation result for upd(20)mat was found (Eggermann 2015 says no other
molecular alteration is reported), so I removed the old upd(20)mat claim rather
than reversing it. Marked unresolved. The text says "no established
interpretation here".

3.4 Copy-number effects on direction. Several copy-number states move the
methylation fraction without being the named disorder: paternal 6q24
duplication (about one third methylated; this is TNDM, and the TNDM LOW text
does not mention it), paternal 11p15 duplication (the consensus statement says methylation is abnormal; ICR1 up and ICR2 down is my inference from copy number), and, by
inference, maternal 15q11-q13 duplication (increased SNURF methylation without
PWS) and maternal 11p15 duplication. The report cannot tell these from
epimutation or UPD by methylation alone. I did not add caveats because I did not
open a source for the maternal duplications. A general sentence recommending
dosage analysis could be added at clinical sign-off.

3.5 MLID description. The sentence "Pathogenic variants in imprinting
maintenance factors (ZFP57, NLRP5, KCNQ1OT1, DPPA3) can cause MLID" is not
direction-related and was left unchanged. KCNQ1OT1 is a cis-acting locus, not a
trans-acting factor, and I found no source for DPPA3. The 2022 compilation of
trans-acting MLID genes (PMID:35296332) lists ZFP57, NLRP2, NLRP5, NLRP7, PADI6,
KHDC3L, OOEP, ZAR1, ARID4A, UHRF1 and ZNF445. The OMIM number given for MLID
(617107) could not be checked (OMIM was not reachable).

3.6 Angelman deletions with normal methylation. The user asked whether some AS
deletions leave methylation unaffected. In the sources opened, the ordinary
(class I / II) deletions and imprinting-centre deletions give the AS (hypo)
pattern, which is why the old "15q11-q13 deletion AS cases may show normal ..."
sentence was removed. No source opened documents a maternal deletion that spares
the SNRPN region and leaves methylation normal; such cases (for example an
UBE3A-only deletion) are not excluded by anything I read. The new text only claims
what is confirmed: UBE3A variants and cases without an identified cause have
normal methylation.

3.7 No hypermethylation mechanism found (marked probable, not confirmed). TNDM
HIGH at PLAGL1 (one MLID family with gain of methylation, PMID:35296332),
BWS HIGH at KCNQ1OT1, MLID HIGH at PEG3, SRS LOW at GRB10/MEST (reported in
MLID). The unchanged HIGH/LOW texts for TNDM, BWS and MLID state these are not
recognised causes, which rests on absence from the reviews consulted.

3.8 Percentages left unchanged. SRS ICR1 loss "about 40%" (consensus statement:
30-60%) and upd(7)mat "about 10%" (consensus: 5-10%) in the SRS description.

## 4. Statements checked and left unchanged

| Entry / region | Result |
|---|---|
| TNDM, PLAGL1:alt-TSS-DMR | Maternal methylation; LOW = TNDM. Consistent. PMID:20803656, 41153775 |
| BWS, KCNQ1OT1:TSS-DMR | Maternal methylation (ICR2); LOW = BWS (about 50%). Consistent. PMID:29377879, 27585961 |
| BWS/SRS, H19/IGF2:IG-DMR | Paternal methylation; LOW = SRS, HIGH = BWS (5-10%). Consistent. PMID:27585961, 29377879 |
| Kagami-Ogata/Temple, DLK1/MEG3:IG-DMR and MEG3:TSS-DMR | Paternal methylation; HIGH = Kagami-Ogata, LOW = Temple. Consistent (apart from items in 2.3). PMID:26377239, 24891339 |
| MLID, PEG3:TSS-DMR | Maternal methylation; LOW = MLID. Consistent. PMID:19397955, 35296332 |
| GNAS-NESP, GNAS-AS1, GNASXL, GNAS A/B region labels | NESP paternal, AS1 and XL maternal: unchanged and correct. |
| SNURF and MAGEL2 labels | Maternal: unchanged and correct. |
| PLAGL1, MEST, KCNQ1OT1, PEG3, DLK1/MEG3 IG-DMR, MEG3 TSS, H19 labels | Unchanged and correct. |

## 5. How the sources were read

Full texts were read from PubMed Central and Europe PMC (Beygo 2019,
Ramsden 2010, Eggermann 2015, Brioude 2018, Mantovani 2018, Court 2014,
Baena 2024, D'Angelo 2025, Urakawa 2025, Eggermann 2022, Bastepe 2011, Huang 2009,
Lin 2024, Pietrusinski 2025, McCullough 2024) and from the publisher page
(Wakeling 2017). Abstracts only: Mackay and Temple 2010, Mackay 2008, Arnaud 2003,
Riesewijk 1997, Ogata and Kagami 2016, Kagami 2015, Ioannides 2014, Li 2024. Not opened: Monk 2018 nomenclature recommendations (PMID:27911167;
only the abstract was available), OMIM, GeneReviews (blocked). Every PMID in the
table was matched to its title by a database lookup before it was entered.

## 6. Other observations

* The package's own suite had 46 tests, not 34; all passed before the change.
* `DISORDER_GROUPS["SRS"]` lists H19/IGF2:IG-DMR, but that region belongs to
  the "BWS/SRS" disorder key, so the SRS card never shows H19 results. The SRS
  text therefore now refers the reader to the H19 result instead of interpreting it.
* The H19/IGF2:IG-DMR region is marked unreliable and never flagged; its rows in
  the table describe biology, not what the tool will report.
