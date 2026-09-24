# Archive: Pre-Revision (April 2026) Simulation Outputs

**Date archived:** 2026-09-24
**Archived by:** step 0.5, B5 of `docs/plans/simulation-revision-plan.md`

## What is archived here

Everything produced from the April 2026 tau-sweep simulation and its predecessors:
- **Source data:** `sim_data/sim_results_MLE_tau_sweep_*_20260408_191916.rds` (12,800
  scenarios) and `sim_data/sim_results_MLE_*_20260322_151030.rds` (τ = 1 baseline).
- Derived tables, figures and rendered reports built from them.

## Why superseded

The step 0 chapter review found implementation oversights in that simulation. The
2026-09 revision (spec: `docs/plans/simulation-revision-spec.md`) fixes them:
- designs and outcomes now use the same incidence surface (M1);
- Poisson clusters have 100,000 people, not 1,000 (M2);
- tied incidence is broken at random (M3);
- seeding is explicit and key-based (M4);
- Monte Carlo SEs come from surface-level means;
- aliasing is flagged instead of hidden (M6/M7);
- a non-oracle sensitivity estimator is added (M8).

New results: `results/sim_data/sim_results_MLE_tau_sweep_combined_<2026-09 timestamp>.rds`
(oracle) and `sim_results_MLEnonoracle_tau_sweep_combined_<timestamp>.rds`.

**These numbers must not be cited** (e.g. MSE 0.079 / 0.802, Friedman χ² 557.3–800.9).

## Files archived

All moves were made with `git mv`, so the files are tracked here as they were at the old
locations (the plan expected the `.rds` files to be untracked; they weren't).

| File | Original location |
|------|-------------------|
| `sim_data/sim_results_MLE_tau_sweep_{combined,iid,spatial,poisson}_20260408_191916.rds` | `results/sim_data/` |
| `sim_data/sim_results_MLE_{combined,iid,spatial,poisson}_20260322_151030.rds` | `results/sim_data/` |
| `six_design_manuscript/` | `results/` |
| `eight_design_supplementary/` | `results/` |
| `mle_per_config/` | `results/` |
| `00_mathematical_specification.{pdf,html}`, `07_results_summary.{pdf,html}`, `09_MLE_design_recommendation_report.{pdf,html}`, `11_statistical_comparisons_report.{pdf,html}` | `results/` |
| `MLE_*.pdf` (5 files), `tau_sweep_run_*.log` (3 files) | `results/` |
| `figures/design_samples_{8panel,option1_overlays}.{pdf,png}` | `results/figures/` |
| `code_sources/07_results_summary.qmd`, `code_sources/09_MLE_design_recommendation_report.Rmd` | `code/` (stale twins of the live `07_…Rmd` / `09_…qmd`) |

## Not moved

- `results/sim_data/sim_results_DIM_*_20260304_195321.rds`: the DIM baseline wasn't
  re-run (out of scope for the revision). It shares the pre-revision DGP flaws and only
  covers 6 designs.
- `paper/ctj_manuscript/Supplementary_Information.tex` (Section S1) still cites
  `results/six_design_manuscript/…` and `results/eight_design_supplementary/…`. Those
  paths now live under this folder. The CTJ/SI stay untouched until manuscript step 5.
