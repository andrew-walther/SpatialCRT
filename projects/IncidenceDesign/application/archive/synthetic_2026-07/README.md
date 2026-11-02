# Synthetic-incidence application (superseded 2026-10-02)

The July 2026 application runs used synthetic SUD incidence and the pre-revision
design rules. The observed-incidence application (`../../results/real_sud_rev_20261002/`)
and the methods extensions (`../../results/methods_ext_2026-10/`) replace them.
Moved here on 2026-11-02; nothing was deleted.

- `results/{smoke,pilot,full}/`: outputs of `application/code/run_application_profiles.R`.
  Re-running that script writes to `application/results/<profile>/` again.
- `report/proposal_report/`: the synthetic-data proposal report (qmd, html, pdf).

`code/14_manuscript_supplement_figures.R` section (d) still reads
`results/full/application_full_results.rds` from here for its STALE application table.
