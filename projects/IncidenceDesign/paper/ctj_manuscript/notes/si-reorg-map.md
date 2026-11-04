# SI reorganization map (2026-11-04)

The SI was rebuilt into Part I Supplementary Methods and Part II Supplementary
Results, following the DeSurv SI Appendix layout (Contents, List of Figures,
List of Tables, then the two parts with sections numbered straight through).
Every paragraph, table and figure of the old file was moved as a block; none was
dropped. The `app:*` labels shared with Chapter 3 were kept.

| New section | Old source |
|---|---|
| 1 Simulation setting | S1.1, S1.2 |
| 2 Treatment assignment designs | S1.3; S2.1 = first half of old S4 (excluded-variant definitions, `tab:excluded-designs`) |
| 3 Estimation and identification | S1.4; S3.1 = old S9.1 identification diagnostics |
| 4 Simulation design, resampling and metrics | S1.5; S4.1 = test/reporting paragraphs of S1.5; S4.2 = design paragraph of the allocation-risk study (old S2.9); S4.3 = old S5 simulation mechanics |
| 5 Observed-incidence application methods | S3.1–S3.5; regional criterion, `tab:nc-order`, `fig:nc-regions`, `fig:nc-examples` from old S10; S5.4 = reproducibility/precision paragraphs of old S10; estimation-model and spillover-form definitions from old S11 |
| 6 Estimation error and its dependence on τ | old S9 intro + S9.2 |
| 7 Grid benchmark under queen contiguity | S2.1, S2.3, S2.4, `fig:coverage`; S7.4 = second half of old S4 (`tab:eight-design`) |
| 8 Formal statistical comparisons | S2.2 + old S6 critical-difference diagrams |
| 9 Stability across treatment effect sizes | S2.7 + old S6 τ-coverage figure |
| 10 Rook contiguity | S2.5 + old S7 rook figures |
| 11 Non-oracle sensitivity analysis | S2.6 + old S8 table |
| 12 Robustness, SRS benchmark, allocation-specific risk | S2.8, S2.9 (results paragraph) |
| 13 Application: accuracy, inference and reach | S3.6 + results figures/tables of old S10 (`fig:nc-pop`, `fig:nc-setting-ratios`, `fig:nc-mse-detail`, `tab:nc-annual`) + targeting paragraph of old S11 |
| 14 Application: allocation downside | S3.7 + `fig:nc-risk`, `tab:nc-summary-sens` from old S10 |
| 15 Application: participation, power and sensitivity | S3.8 + results of old S11 (power, calibration, guidance, formal tests, model and spillover-form results) |

Old → new for the `app:*` labels: eight-design S4 → 7.4; mechanics S5 → 4.3;
supp-comparisons S6 → 8; rook-figures S7 → 10; oracle-nonoracle S8 → 11;
technical/technical-notes S9 → 6; nc S10 → 5.4; nc-ext S11 → 15. New label
`app:variants` (2.1).

New sections are numbered plainly (1, 2, …, as in the DeSurv SI); figures,
tables and equations keep the S prefix. "Old" numbers in the right-hand column
are the pre-2026-11-04 S-numbered sections.
