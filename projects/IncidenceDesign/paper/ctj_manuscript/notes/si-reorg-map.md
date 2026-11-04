# SI reorganization map (2026-11-04)

The SI was rebuilt into Part I Supplementary Methods and Part II Supplementary
Results, following the DeSurv SI Appendix layout (Contents, List of Figures,
List of Tables, then the two parts with sections numbered straight through).
Every paragraph, table and figure of the old file was moved as a block; none was
dropped. The `app:*` labels shared with Chapter 3 were kept.

| New section | Old source |
|---|---|
| S1 Simulation setting | S1.1, S1.2 |
| S2 Treatment assignment designs | S1.3; S2.1 = first half of old S4 (excluded-variant definitions, `tab:excluded-designs`) |
| S3 Estimation and identification | S1.4; S3.1 = old S9.1 identification diagnostics |
| S4 Simulation design, resampling and metrics | S1.5; S4.1 = test/reporting paragraphs of S1.5; S4.2 = design paragraph of the allocation-risk study (old S2.9); S4.3 = old S5 simulation mechanics |
| S5 Observed-incidence application methods | S3.1–S3.5; regional criterion, `tab:nc-order`, `fig:nc-regions`, `fig:nc-examples` from old S10; S5.4 = reproducibility/precision paragraphs of old S10; estimation-model and spillover-form definitions from old S11 |
| S6 Estimation error and its dependence on τ | old S9 intro + S9.2 |
| S7 Grid benchmark under queen contiguity | S2.1, S2.3, S2.4, `fig:coverage`; S7.4 = second half of old S4 (`tab:eight-design`) |
| S8 Formal statistical comparisons | S2.2 + old S6 critical-difference diagrams |
| S9 Stability across treatment effect sizes | S2.7 + old S6 τ-coverage figure |
| S10 Rook contiguity | S2.5 + old S7 rook figures |
| S11 Non-oracle sensitivity analysis | S2.6 + old S8 table |
| S12 Robustness, SRS benchmark, allocation-specific risk | S2.8, S2.9 (results paragraph) |
| S13 Application: accuracy, inference and reach | S3.6 + results figures/tables of old S10 (`fig:nc-pop`, `fig:nc-setting-ratios`, `fig:nc-mse-detail`, `tab:nc-annual`) + targeting paragraph of old S11 |
| S14 Application: allocation downside | S3.7 + `fig:nc-risk`, `tab:nc-summary-sens` from old S10 |
| S15 Application: participation, power and sensitivity | S3.8 + results of old S11 (power, calibration, guidance, formal tests, model and spillover-form results) |

Old → new for the `app:*` labels: eight-design S4 → S7.4; mechanics S5 → S4.3;
supp-comparisons S6 → S8; rook-figures S7 → S10; oracle-nonoracle S8 → S11;
technical/technical-notes S9 → S6; nc S10 → S5.4; nc-ext S11 → S15. New label
`app:variants` (S2.1).
