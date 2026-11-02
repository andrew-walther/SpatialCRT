# Figures not used by the current CTJ or SI

Moved here on 2026-11-02 so `figures/` holds only what `CTJ_Manuscript.tex`,
`CTJ_Exhibits.tex` and `Supplementary_Information.tex` include (or what
`presentation/tools/refresh_figures.sh` reads). Nothing was deleted.

| Path | What it is | Producer |
|---|---|---|
| `application_maps/` | Synthetic-incidence and k-means maps from the July 2026 synthetic application, plus an older service-area map | `application/code/run_application_profiles.R` (superseded) |
| `si_figures/` | Earlier SI copies of the code/14 statistical-comparison figures; the SI no longer includes them | `code/14` (now writes to `results/six_design_manuscript/<nb>/si_figures/`) |
| `design_samples_6panel.pdf` | Old six-design sample figure, replaced by `design_assignments_panels.pdf` | `code/06` |
| `fig_coverage_tau_6design.pdf`, `fig_mse_by_design_6design.pdf`, `fig_mse_by_design_6design_queen_journal.pdf` | Pre-split (pooled) or journal-variant grid figures | `code/12`, `code/14`, `code/19` |
| `real_sud/*_journal.pdf` | Duplicates of the `real_sud/` figures under their export names | `application/code/real_sud_summary.R` |

Re-running a producer may write a fresh copy into `results/`, not here.
