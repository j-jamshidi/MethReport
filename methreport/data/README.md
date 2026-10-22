# Bundled control data

- `controls_hg38.tsv` - in-house cohort, long schema (`region_name`, `position`,
  `sample_id`, `methylation_pct`), unphased rows only. Samples are identified by run number.
- `controls_hg38.provenance.json` - samples, modification counted (5mC), basecaller,
  pipeline and modkit versions, and n per position.

No `controls_t2t.tsv` is bundled. Samples tested against these controls must be
counted as pure 5mC with the same modkit settings. Further files can be built
with `methreport controls build`.
