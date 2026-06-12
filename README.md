# MethReport

**Clinical-grade DNA methylation imprinting disorder reporter**

MethReport analyses Oxford Nanopore modbam files and produces an interactive HTML report showing methylation levels across 14 differentially methylated regions (DMRs) associated with known imprinting disorders.

---

## Features

- Reads **modbam files directly** (BAM with MM/ML tags — no pre-conversion needed)
- Supports **haplotype-phased** analysis (HP1 / HP2 tags)
- Covers **14 imprinting DMRs** across TNDM, SRS, BWS, Kagami-Ogata, Temple, Prader-Willi, Angelman, PHP, and MLID
- Two reference genomes: **T2T-CHM13v2.0** and **hg38**
- Outputs:
  - Interactive HTML report (Plotly — self-contained, no server needed)
  - Per-region summary TSV
  - Per-CpG methylation TSV
  - IGV-ready BED tracks (unphased + HP1 + HP2)
- Bundled control reference ranges with user-override support
- Rich terminal output with progress and flagged-region summary

---

## Covered Disorders & Regions

| Disorder | DMRs |
|---|---|
| TNDM | PLAGL1:alt-TSS-DMR |
| SRS | GRB10:alt-TSS-DMR, MEST:alt-TSS-DMR, H19/IGF2:IG-DMR |
| BWS | H19/IGF2:IG-DMR, KCNQ1OT1:TSS-DMR |
| Kagami-Ogata / Temple | MEG3:TSS-DMR, MEG8:Int2-DMR |
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
methreport run sample.sorted.bam --ref t2t --out results/

# hg38 reference
methreport run sample.sorted.bam --ref hg38 --out results/

# With custom sample name
methreport run sample.sorted.bam --ref t2t --sample PATIENT_001 --out results/
```

### 3. Open the report

```
results/PATIENT_001_methreport.html
```

---

## Options

```
methreport run [OPTIONS] BAM

  BAM                Indexed modbam file (.bam)

Options:
  --ref  -r TEXT     Reference genome: 't2t' or 'hg38'  [default: t2t]
  --out  -o PATH     Output directory  [default: methreport_output]
  --sample -s TEXT   Sample identifier (default: BAM filename)
  --controls -c PATH User-supplied controls TSV/XLSX
  --replace-controls Replace bundled controls instead of supplementing
  --min-cov INT      Minimum CpG coverage  [default: 5]
  --call-threshold   MM/ML probability threshold  [default: 0.5]
  --no-bed           Skip BED track output
  --no-tsv           Skip TSV export
  --verbose -v       Debug logging
```

---

## Custom Controls

To use your own cohort controls, prepare a TSV with columns:

```
region_name    position    sample_id    methylation_pct
PLAGL1_alt_TSS_DMR    145200350    ctrl_01    48.2
...
```

Then pass it with `--controls my_controls.tsv`.

To convert the original NanoImprint Excel controls:

```bash
python scripts/convert_controls.py \
    --t2t  ctrls_2.xlsx \
    --hg38 ctrls_hg38.xlsx \
    --out  methreport/data/
```

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
