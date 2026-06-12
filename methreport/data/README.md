# Bundled control data

Place control TSV files here with names:
- `controls_t2t.tsv`
- `controls_hg38.tsv`

Expected columns: `region_name`, `position`, `sample_id`, `methylation_pct`

These can be generated from the original NanoImprint Excel control files using
the `scripts/convert_controls.py` helper script.
