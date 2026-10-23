# MethReport

**Clinical-grade DNA methylation imprinting disorder reporter**

MethReport analyses Oxford Nanopore modbam files and produces an interactive HTML report showing methylation levels across 15 differentially methylated regions (DMRs) associated with known imprinting disorders.

---

## Features

- Reads **modbam files directly** (BAM with MM/ML tags — no pre-conversion needed)
- Supports **haplotype-phased** analysis (HP1 / HP2 tags)
- Covers **15 imprinting DMRs** across TNDM, SRS, BWS, Kagami-Ogata, Temple, Prader-Willi, Angelman, PHP, and MLID
- Two reference genomes: **T2T-CHM13v2.0** and **hg38**
- Outputs:
  - Interactive HTML report (Plotly and webfonts are fetched from a CDN, so
    charts need network access when the report is opened)
  - Per-region summary TSV
  - Per-CpG methylation TSV
  - IGV-ready BED tracks (unphased + HP1 + HP2)
- Control reference ranges supplied at run time with `--controls`. **No control
  file ships with the package** — see Controls below
- Rich terminal output with progress and flagged-region summary

---

## Covered Disorders & Regions

| Disorder | DMRs |
|---|---|
| TNDM | PLAGL1:alt-TSS-DMR |
| SRS | GRB10:alt-TSS-DMR, MEST:alt-TSS-DMR, H19/IGF2:IG-DMR |
| BWS | H19/IGF2:IG-DMR, KCNQ1OT1:TSS-DMR |
| Kagami-Ogata / Temple | DLK1/MEG3:IG-DMR, MEG3:TSS-DMR, MEG8:Int2-DMR |
| Prader-Willi / Angelman | MAGEL2:TSS-DMR, SNURF:TSS-DMR |
| MLID | PEG3:TSS-DMR |
| PHP | GNAS-NESP:TSS-DMR, GNAS-AS1:TSS-DMR, GNASXL:Ex1-DMR, GNAS A/B:TSS-DMR |

---

## Installation

```bash
pip install methreport
```

Or from source:

```bash
git clone https://github.com/j-jamshidi/MethReport.git
cd MethReport
pip install -e ".[dev]"
```

**Requirements:** Python ≥ 3.10, samtools (for BAM indexing)

---

## Quick Start

### 1. Prepare your modbam file

Your BAM must be:
- Produced with ONT's `wf-human-variation` using `--mod` and `--phased` flags
- Sorted and indexed:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

### 2. Run MethReport

```bash
# T2T-CHM13v2.0 reference (recommended)
methreport run --bam sample.sorted.bam --ref t2t --out results/

# hg38 reference
methreport run --bam sample.sorted.bam --ref hg38 --out results/

# With custom sample name
methreport run --bam sample.sorted.bam --ref t2t --sample PATIENT_001 --out results/
```

### Prefer modkit BEDMethyl over the BAM

`--bam` parses MM/ML tags directly, which merges strands less accurately than
modkit does. On a 30x phased HG002 BAM the two paths differ materially: mean
methylation read 5-17 percentage points lower through `--bam`, the mean z-score
against controls sat at -0.93 rather than -0.20, and one region was flagged
abnormal that the BEDMethyl path called normal with both methods agreeing.
Phasing also survives better, with per-haplotype evidence available for 10 of
14 regions against 6 of 14.

Use `--bed` where you can:

```bash
modkit pileup sample.bam unphased.bed --cpg --combine-strands --ref genome.fa
modkit pileup sample.bam hp_dir      --cpg --combine-strands --ref genome.fa \\
    --partition-tag HP --prefix hp

methreport run --bed unphased.bed \\
    --bed-hp1 hp_dir/hp_1.bed --bed-hp2 hp_dir/hp_2.bed \\
    --ref hg38 --out results/
```

### 3. Open the report

```
results/PATIENT_001_methreport.html
```

---

## Options

```
methreport run [OPTIONS]

Input (choose one):
  --bam PATH         Indexed modbam file (.bam)
  --bed PATH         modkit BEDMethyl, unphased (recommended)
  --bed-hp1 PATH     modkit BEDMethyl for haplotype 1 (used with --bed)
  --bed-hp2 PATH     modkit BEDMethyl for haplotype 2 (used with --bed)
  --list PATH        Text file of BAM paths, one per line, for batch analysis

Options:
  --ref  -r TEXT     Reference genome: 't2t' or 'hg38'  [default: t2t]
  --out  -o PATH     Output directory  [default: methreport_output]
  --sample -s TEXT   Sample identifier (default: BAM filename)
  --controls -c PATH User-supplied controls TSV/XLSX
  --replace-controls Discard any controls found in methreport/data/
  --require-controls Fail if no controls are loaded, instead of judging by
                     allelic asymmetry alone
  --min-cov INT      Minimum CpG coverage  [default: 5]
  --call-threshold   MM/ML probability threshold  [default: 0.5]
  --no-bed-tracks    Skip BED track output
  --no-tsv           Skip TSV export
  --verbose -v       Debug logging
```

---

## Controls

Two schemas are accepted and detected by column name:

```
long        region_name  position  sample_id  methylation_pct
aggregated  region_name  position  mean  sd  [n_controls]
```

Use the long form for your own cohort, one row per CpG per control sample:

```
region_name           position     sample_id   methylation_pct
PLAGL1_alt_TSS_DMR    145200350    ctrl_01     48.2
```

Pass either form with `--controls my_controls.tsv`, and add
`--replace-controls` to discard any file found in `methreport/data/` rather than
supplement it. Add `--require-controls` for diagnostic runs, so that a mistyped
path or a missing file is an error instead of a quieter result. Two aggregated summaries of the same position cannot be pooled
without each one's control count, so where both files cover a position the
user's value is taken and a warning is logged.

### Building an in-house cohort

Analyse each control sample with `methreport run`, then pool the per-CpG TSVs it
wrote (unphased rows only) into the long schema, with a provenance record:

```bash
methreport controls build results/*/*_cpg_methylation.tsv \
    --ref hg38 --out controls/controls_hg38.tsv \
    --mod-code 5mC --basecaller-model <model> --pipeline-version <version>
```

The provenance file (`controls_hg38.provenance.json`) lists the samples, the
modification counted, the basecaller and pipeline versions, and the number of
controls at each position. Use the same modification code for every sample
analysed against these controls.

### Obtaining the NanoImprint controls

No control file ships with this package. A cohort built from clinical samples is
patient-derived data and does not belong in a code repository, so keep yours in
your laboratory's reference directory and pass it with `--controls`. Without
controls every region falls back to allelic asymmetry between the haplotypes,
which `confidence` records as `asymmetry only`.

If a cohort ever has to be shared, share the aggregated form. A long-form file is
reduced to per-position mean, SD and count on load, so an aggregated file gives
identical reference ranges without one row per patient. Keep the long form in the
laboratory, since only it carries `n` per position and can be extended sample by
sample.

The published NanoImprint
spreadsheets give per-position mean and SD, with no per-sample values and no
control count, so they convert to the aggregated schema:

```bash
git clone https://github.com/carolinehey/NanoImprint
python scripts/convert_controls.py \
    --t2t  NanoImprint/ctrls_2.xlsx \
    --hg38 NanoImprint/ctrls_hg38.xlsx \
    --out  controls/
```

That yields 2,196 positions across all 14 hg38 regions, and 2,211 for T2T.
Positions are assigned to a DMR by coordinate using this package's own region
definitions, so the two cannot drift apart; a position falling outside every
region is reported and dropped.

**The NanoImprint repository carries no licence**, so its control data is all
rights reserved. Convert it for local use; do not redistribute the converted
file or ship it inside a published image without the authors' permission.

---

## Other Commands

```bash
# Validate a BAM has modbam tags
methreport validate sample.bam

# List all DMR regions
methreport list-regions --ref t2t
```

---

## Output Files

| File | Description |
|---|---|
| `<sample>_methreport.html` | Interactive HTML report |
| `<sample>_summary.tsv` | Per-region mean methylation and flags |
| `<sample>_cpg_methylation.tsv` | Per-CpG values (unphased + HP1 + HP2) |
| `<sample>_unphased.bed` | IGV track — all reads |
| `<sample>_hp1.bed` | IGV track — haplotype 1 |
| `<sample>_hp2.bed` | IGV track — haplotype 2 |

---

## Data Source

Region coordinates and control data derived from
[NanoImprint](https://github.com/carolinehey/NanoImprint) by Caroline Hey.

---

## License

MIT
