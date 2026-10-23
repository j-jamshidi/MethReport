# Control data

No control file ships with this package, deliberately: a cohort built from clinical
samples is patient-derived data and does not belong in a code repository.

Supply one at run time with `--controls`, in either accepted schema:

    long        region_name  position  sample_id  methylation_pct
    aggregated  region_name  position  mean  sd  [n_controls]

Without controls, `z_score_flag` has nothing to compare against and every region falls
back to the allelic asymmetry between the two haplotypes, which `confidence` records as
`asymmetry only`. See the Controls section of the top-level README for how to build a
cohort and for the licensing position on the published NanoImprint data.
