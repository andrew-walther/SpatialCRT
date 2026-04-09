# Archive: Pre-Tau-Sweep Deliverables

**Date archived:** 2026-04-08
**Archived by:** Post-simulation update workflow

## What is archived here

Reports and figures generated from the MLE_combined baseline:
- **Source data:** `sim_results_MLE_combined_20260322_151030.rds`
- **Scenarios:** 2,560 (τ fixed at 1.0)
- **Date produced:** 2026-03-22

## Why superseded

The tau-sweep simulation completed on 2026-04-08, producing 12,800 scenarios with:
- **True_Tau** swept across {0.8, 1.0, 1.5, 2.0, 3.0}
- **N_Valid_Est** — Monte Carlo SE denominator (count of non-NA estimates per scenario)
- **Power** — P(reject H₀: τ=0) = fraction of CIs with lower bound > 0

New results file: `sim_results_MLE_tau_sweep_combined_20260408_191916.rds`

## Files archived

| File | Original location |
|------|-------------------|
| `00_mathematical_specification.{pdf,html}` | `results/` |
| `07_results_summary.{pdf,html}` | `results/` |
| `09_MLE_design_recommendation_report.{pdf,html}` | `results/` |
| `11_statistical_comparisons_report.{pdf,html}` | `results/` |
| `MLE_combined_design_recommendations.pdf` | `results/` |
| `MLE_combined_incidence_overview.pdf` | `results/` |
| `MLE_statistical_comparisons.pdf` | `results/` |
| `MLE_combined_*.pdf` (5 per-config files) | `results/mle_per_config/` |
| `IncidenceSpatialCRT_Report.{html,pdf}` | `paper/report/` |
