# Step 0 review findings — `Dissertation_Chapter.qmd`, first pass (2026-09-23)

Produced by a fresh reviewer run of `step0-chapter-review-prompt.md` against
the chapter as snapshotted in
`archive/IncidenceDesign_Manuscripts_PreUnification_2026-09-23/` (line numbers
refer to that version). **Status: findings reported, none applied yet** —
several need author decisions (see "Decisions needed" at the bottom).

Totals: 17 ERROR, 31 WARN, 17 STYLE, 13 NEEDS-AUTHOR-CONFIRMATION.

Sources beyond the three ground-truth `.txt` files that the reviewer used for
some findings: `results/sim_data/sim_results_MLE_tau_sweep_combined_20260408_191916.rds`,
`results/six_design_manuscript/six_design_comparison_report.rds`,
`results/eight_design_supplementary/eight_design_summary.txt`, and scratch
re-runs of `03_designs.R` / `02_incidence_generation.R`.

**Independently re-verified by the executing session** (not just the reviewer's
word): #11/#27/#32/#36 (rook vs queen split — sim rds, tau=1, six designs:
Checkerboard queen MSE 1.118 / cov 0.921 / bias +0.047; rook MSE 0.487 / cov
0.133 / bias −0.657; rook bias by gamma = −0.508, −0.607, −0.707, −0.807;
queen MSE flat 1.117–1.118), #15 (incidence regenerated per tau inside
`run_incidence_config()`; X_matrix column k used for outcome replicate k),
#23 (X_matrix and Epsilon_matrix drawn before the per-scenario `set.seed`),
#9 (Chapter 2 draft L1312–1354: "checkerboard paradox", β̂ ≈ β − ψ under
exact collinearity).

## Findings

Categories: 1 citations, 2 numbers, 3 method/fact vs code, 4 internal
consistency, 5 overclaim, 6 standing rules/mannered, 7 voice.

### Citations
1. L38 WARN — "spatial Durbin model formalized by Ord [@ord_estimation_1975]": lit review cites Ord 1975 only for ML computational cost; `lesage_introduction_2009` is in the master bib. Re-attribute.
2. L38 WARN (minor) — Moran's I: lit review cites both `moran_interpretation_1948` and `moran_notes_1950`.
3. L38 WARN — Jarvis 2017 is a review of *analysis* methods; the "little guidance on assignment" point is the author's inference (lit review L314–318), not a Jarvis finding.
4. L40 WARN — covariate-constrained randomization: method is Moulton 2004 (`moulton_covariatebased_2004`); Crisp 2023 extends it to cluster selection and notes it increases spatial dispersion (lit review L374–379).
5. L40 ERROR — McCann 2018 "fried-egg ... ring of untreated units": fried-egg omits a cluster's outer area from analysis; McCann proposed whole-cluster inclusion/exclusion (lit review L414–418).
6. L40 WARN (minor) — `cai_independentset_2023` is a preprint (`@misc`); "published or posted".

### Project 1 / cross-chapter
7. L44 ERROR — Project 1 used a SAR (lag) model with spillover indicator, not a spatial Durbin model (P1 draft L518–549).
8. L44 WARN — Project 1's future-considerations list does not include heterogeneous baseline risk; restrict "explicitly anticipated" to what it lists.
9. L44, L292 WARN — Checkerboard is Chapter 2's block-stratified checkerboard; its failure here under rook is Chapter 2's "checkerboard paradox" (β̂ ≈ β − ψ). The chapter never makes this link and instead blames an intuition "imported from non-spatial CRT practice".
10. L44, L72 WARN — "reported in the supplementary material": the chapter has no supplement (that's CTJ's SI). Point to the chapter appendix (step 3) instead.
62. L44, L66 STYLE — "Project 1 of this dissertation" / "The earlier proof-of-concept study" → "Chapter 2" (lit review and P1 draft convention).

### Method / fact vs code
11. L54, L296 ERROR — "adjacency defined by queen contiguity": code sweeps rook and queen; headline numbers pool both. Checkerboard coverage 0.921 queen vs 0.133 rook.
12. L56–62 WARN — Eq. 1 shows only γWZ; control-only regime uses γ·WZ·(1−Z).
13. L60 ERROR — at ρ=0 the model is not a "spatial lag-of-X model"; only treatment is lagged.
14. L66 WARN — incidence transformations omitted: iid Uniform(0,1); spatial = SAR-filtered normals through `pnorm`, always built with queen W; Poisson rates rank-normalized (rank−0.5)/N, ~0.35 expected deaths/cluster → ~5 distinct values (heavy ties).
15. L68 ERROR — incidence is not held fixed across all replications: 10 columns, outcome replicate k uses column k (designs see column 1 only); regenerated for each τ; neighbor type also varies within a configuration.
16. L72 WARN — Balanced Halves vs Balanced Quartiles: Nemenyi p=0.816 but Wilcoxon-Holm p=1.34e-6, so "statistically indistinguishable" is wrong for that pair; rationale is lower MSE/rank. Give separate rationale per excluded pair.
17. Table 1 L84 ERROR — Incidence-Guided Saturation Quadrants uses literal 5×5 coordinate quadrants (not k-means regions — that's the application's "Saturation Regions"); saturations fixed {0.8,0.6,0.4,0.2} by rank of quadrant mean incidence; exactly 50 treated.
18. Table 1 L87 ERROR — Isolation Buffer % treated: rook 34–44 (mean 38), queen 17–25 (mean 21.6), not "20–30".
19. Table 1 L86; L133 WARN — High Incidence Focus treats incidence > median; with tied Poisson ranks 22–50 treated (mean ~35). The Poisson degradation explanation omits this arm imbalance.
20. L97 WARN — estimator is `lagsarlm` spatial-lag model `Y ~ Z + Spill + X` with the regime-specific true spill covariate, not a "spatial Durbin specification".
21. L101 ERROR — τ=1 count is 5×2×4×4×2×6 = 1,920 scenarios, not 384; 320 blocks = 1,920/6 and every block is complete. Delete the parenthetical.
22. L103, L105 WARN — deterministic designs are copied 25× against the same noise (250 recorded estimates = 25 copies of 10 values); noise drawn once per (nb type, ρ) and reused across designs/γ/regimes (crossed, not independent). Contradicts "no reduced effective sample size".
23. L103 WARN — per-scenario seed is set after X_matrix and Epsilon_matrix are drawn, so "exact reproducibility independent of execution order or parallelization" is overstated.
24. L105 ERROR — 12,800 is scenarios, not iterations (3.2M iterations); Fail_Rate counts errors/NA with warnings suppressed, so "converged" overstates.
25. L105 WARN (low) — MSE computed as mean((τ̂−τ)²); power is lower CI bound > 0 (one-sided). Align definitions.

### Results / numbers / overclaims
26. L111, Tables 2–7 WARN — averages pool 5 incidence configs and both neighbor types, against project CLAUDE.md Invariant 2. Author decision.
27. L111, Fig 2 WARN — 52.7% coverage averages queen (~0.92) and rook (~0.13); Fig 2 has queen/rook panels the caption doesn't mention.
28. L129 ERROR — Incidence-Guided Saturation Quadrants vs Balanced Quartiles Nemenyi p=0.0001 (significant); it is significantly better than every other design.
29. L135 WARN — incidence×ρ interaction explanation for Checkerboard is unsupported; label as conjecture.
30. L139 ERROR — "Five of six flat, < 0.03": High Incidence Focus rises 0.110→0.157 (+43%, monotone). Four of six.
31. L160, L292 ERROR — "maximize physical separation between treated and control ... maximizes treated-control neighbor pairs" is self-contradictory; Checkerboard maximally intersperses arms (under rook WZ = 1−Z exactly). The "separation/buffer" intuition is Isolation Buffer.
32. L160, L292 WARN — γ trend is entirely rook: bias = −γ from exact collinearity, variance ~0.04; queen MSE flat ~1.12. Mechanism is lost identifiability, not leakage.
33. L181 ERROR — Isolation Buffer is identical across regimes (0.1479 vs 0.1476) by construction (no treated cluster has a treated neighbor). "Every design" is wrong.
34. L181 WARN — estimator includes the exact regime-specific spill covariate, so there is no omitted-variable bias; Checkerboard's control-only excess is variance (queen); rook MSE identical across regimes (0.487); rook bias larger under "both".
35. L204 WARN — Friedman rejections don't show ranking stability; mean-MSE order of High Incidence Focus and Isolation Buffer swaps at τ=0.8 and 2.0; average-rank order is constant (report rds); at τ=3 top three not significantly different. Cite average-rank stability explicitly.
36. L204 ERROR — "ρ varied fourfold" is wrong (0 → 0.50); pooled coverage 51.3–54.3% hides rook 12.0–14.5% vs queen 88–96.5%.
37. L243 ERROR — High Incidence Focus is the worst design in 17 of 320 configurations.
38. L243 WARN — iid incidence is not "smooth"; wins 30/64 (iid), 30/64 (spatial 0.20), 28/64 (spatial 0.50): not "decisive".
39. Fig 1 L249 WARN — chapter's figure is stale (alphabetical order); is queen-only, coloured by γ; caption wrong.
40. L266 WARN — "archived with the supplementary results bundle" isn't reader-findable; include clean CD diagram (from `14_manuscript_supplement_figures.R`) or drop.

### Application (mostly moot after step 2 deletes the placeholder run)
41. L272 WARN (UNVERIFIABLE) — Nanavati "validated ... screened against expert physician panel review".
42. L272 WARN (UNVERIFIABLE) — Mirzaei NC base rate / 20–64 / "exceed any individual cancer" (bib title is US-scale); age band differs from dataset's 18–64 (L280).
43. L272 WARN (UNVERIFIABLE, low) — Mounsey "inverse relationship" vs lit review's "varied with household income".
44. L276 WARN — γWZ is contiguity-based; CC coordination network isn't. Soften "exactly as".
45. L286 ERROR — application MSE range: Balanced Quartiles 0.070 (outside 0.075–0.078); 2×2 Blocking 0.078 is in top tier.
46. L286 ERROR — High Incidence Focus does not "ignore incidence"; Isolation Buffer (0.265) also > 0.23.
47. L286, L296 WARN — application ranking differs materially from grid (2×2 5th→3rd; High Incidence Focus 3rd→last; Checkerboard analogue coverage 0.930).
48. L292 WARN (low) — "first systematic ... comparison" novelty claim unverifiable.

### Style / mannered prose / voice
49. L40 — "begun to converge on the exact problem", "signals that the field is actively moving".
50. L48 — padded roadmap sentence.
51. L68 — "Critically,".
52. L111 — "not merely a loss of precision; it represents..." restates.
53. L133 — "Two patterns are notable", "remarkably stable", "more than any single aggregate MSE number, is the clearest evidence".
54. L139 — stacked hedge.
55. L160 — "A clearer and more interpretable pattern emerges".
56. L181 — "runs counter to a naive intuition", "a substantive finding independent of...".
57. L222 — "it is useful to ask", "understates just how badly".
58. L243 — "surfaces a more nuanced picture", "apparent contradiction is resolved"; last sentence restates.
59. L276, L286 — "natural spatial and organizational unit", "direct, real-world analog", "With that caveat stated clearly".
60. L292, L300 — "directly actionable", "actively harmful", "natural culmination".
61. L296 — "at least broadly robust".
63. Whole file — person mixed ("this chapter" ×20, "we" ×9). Chapter 2 draft uses "we" + "this study"; lit review uses neither.
64. YAML — manuscript front matter (4 authors, affiliations, `\linenumbers`, `article`); handled at step 7 port.
65. Results — tense shifts; Table 2 caption writes "tau" in plain text.

No findings: DIM mentions; numbered design shorthand; unresolved @keys; Tables 2–7 values; 9,600/12,800/320; 13 of 15 (both tests); χ² values; win-rate %; table/figure numbering.

## NEEDS-AUTHOR-CONFIRMATION
1. YAML author list / affiliations / corresponding author.
2. L48 — "a planned SUD prevention pilot to be delivered through NC's community college network".
3. L276 — 58 CCs serve all 100 counties; Durham Tech satellite in Orange County.
4. L276 — CCs "communicate and coordinate programmatic activity".
5. L278 — pilot site Lenoir County CC with spillover to "its Kinston and Jones County satellite sites" (note: Lenoir CC's main campus is in Kinston); intervention description.
6. L278 — intervention delivered county-wide.
7. L278 — process measures (medication adherence, clinic visit timeliness).
8. L280 — ~100,000 NC death certificates 2018–2021, ages 18–64, zip-code pooling, SEER denominators.
9. L280 — IRB-protocol access underway.
10. L280 — planned covariates (median income, food desert, physician access, driving distance).
11. L284 — application re-run on 58-cluster geography, dissolved boundaries, Queen contiguity.
12. L284 — Poisson SAR placeholder, 2024 NC OSBM populations.
13. L296, L300 — "substantially heterogeneous" CC-service-area populations.

## Decisions needed before applying (raised with author 2026-09-23)
See session transcript; recorded here once answered.
