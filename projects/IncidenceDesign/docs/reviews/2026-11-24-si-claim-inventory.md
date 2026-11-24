# SI claim inventory (2026-11-24)

Read-only review of `paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex`
(65 pages, built 2026-11-15; source last changed in `d8639d5`). Companion review:
[2026-11-24-si-review.md](2026-11-24-si-review.md).

Conventions. "SI:n" is a line of the SI source. "C:n" is a line of
`manuscript_ctj/CTJ_Manuscript.tex` and "E:n" a line of `CTJ_Exhibits.tex`. Pages come from the
built table of contents (`.build/Supplementary_Information.toc`); exhibit sizes are estimates.
Recommendation codes: **Keep** (abridged SI as is), **Shorten** (abridged SI, cut text),
**Full** (full-length version only), **Merge** (fold into the named item), **Drop** (duplicate
of a main-text exhibit; no version needs it). "Abr. §" is the section of the proposed abridged
outline in the companion review.

## Main-text pointers into the SI

The main text hard-codes SI numbers. Every one changes under any reorganization (see review (c)).

| C/E line | Pointer | Target now | Claim it supports |
|---|---|---|---|
| C:158 (×2) | SI Section 15.2 | Interval calibration | The df factor was chosen by calibration; t + df covered β better than Wald |
| C:165 | "the SI gives full specifications" | Sections 2, 5.2 | Full design rules and NC adaptations |
| C:193 | SI Section 5.5; SI Section 4 | Partial participation; grid metrics | Power design; grid power one-sided at 0.025 |
| C:204 | SI Figure S8 | Bias–variance, queen | Control-only loss is imprecision, not bias |
| C:206 | SI Section 6 | Error decomposition | ρ̂ error scales with β; Isolation Buffer avoids it |
| C:210 | SI Section 12.3 | Grid allocation risk | Outcome noise = median 63% of IGS risk variance |
| C:215 | SI Table S17 | `tab:nc-year` | Absolute yearly MSE |
| C:222 | SI Section 14 | NC allocation downside | Both-arms tail differences are mostly MC noise |
| C:226 | SI Table S5 | `tab:nc-budgets` | Budgets by design and year |
| C:271 | Supplementary material paragraph | whole SI | Lists contents (maps, budgets, orderings, example assignment, sensitivity analyses) |
| E:15 | Table 1 caption "reported in the SI" | Section 2.1, Table S2 | Two variants performed like their counterparts |
| E:89 | SI Table S18 | `tab:nc-annual` | Absolute MSE, bias, coverage |

Main-text claims with an SI home but no pointer: rook grid results (C:197), omitted spillover
term (C:208, SI 11), allocation-risk study design (C:177, SI 4.2), region-selection criterion
(C:187, SI 5.2), `lagsarlm` agreement (C:158, SI 3), incidence-shift and rook sensitivity
(C:228, SI 14 / Table S20).

## Sections and subsections

| SI § | Pages | Contents | Main-text pointer | Claim supported | Recommendation | Duplication notes |
|---|---|---|---|---|---|---|
| Front matter | 1–8 (8) | Title, intro paragraph, Contents, List of Figures (3 p), List of Tables (3 p) | — | — | **Shorten** to ~1.5 p: short captions in the lists (`\caption[short]{long}`) or drop LOF/LOT | Restricted-data sentence also at SI:650 and SI:837 |
| 1.1 Spatial structure and outcome model | 9 (1) | W, outcome model (S1), regimes (S2), parameter values | C:151–156 (no pointer) | Model definition | **Shorten** to the parts the main text lacks (σ, ρ = 0 case) | Repeats main Eq. 1 and regime text almost word for word |
| 1.2 Heterogeneous baseline incidence | 9–10 (1) | Three incidence mechanisms, base-rate justification, matched surfaces, target estimand | none (C:173 has the formulas) | Incidence generation | **Keep**; becomes the home of trim T1/T11 | C:173 already has the same formulas; after T1 the main text needs a pointer here, not "SI Section 2" as the trim plan says |
| 2 Treatment assignment designs | 10–12 (2.5) | Six design rules (Table S1, Figure S1), Checkerboard history, randomization details | C:165 | Full specifications | **Shorten**: one table for all eight designs; Checkerboard history to Full | Table S1 + Figure S1 duplicate main Table 1 (same rules, same illustrative surface) |
| 2.1 The two related variants | 12–13 (1) | Saturation Quadrants and Balanced Halves, test results, why excluded | E:15 | Variants perform like counterparts | **Shorten** to one paragraph | Rationale repeated in 7.4; "incidence guidance matches..." sentence repeated at SI:1817 and SI:1904; **stale claim** SI:227 (see review (c)) |
| 3 Estimation and identification | 13–14 (1) | Oracle/non-oracle models, own ML routine vs `lagsarlm`, t + df interval, Checkerboard identification | C:158 (no pointer) | Estimator validity; interval form | **Keep**; single home of the interval definition | Interval defined again at SI:567 and SI:817; identification repeated in 3.1 |
| 3.1 Identification diagnostics | 14 (0.75) | Exact WZ = 1 − Z, queen boundary information, Isolation Buffer aliases | C:160 | β − ψ under rook Checkerboard | **Merge** with 3 (cut the paragraph-3 preview) | Rook Checkerboard non-identification restated in ~15 captions and paragraphs |
| 4 Simulation design, resampling, metrics | 14–15 (1.5) | Scenario count, blocks, surfaces, CRN, metrics, MC SEs | C:193 (grid power) | Precision and one-sided grid power | **Shorten** | Seed/CRN text repeated in 4.3 |
| 4.1 Statistical comparisons and reporting | 15–16 (0.5) | Friedman, Nemenyi, Holm–Wilcoxon, 50-unit tests, reporting rules | none | Test methods | **Keep** | Same tests restated for NC at SI:862 |
| 4.2 Allocation-risk study | 16 (0.3) | Grid allocation-risk design | none (C:177) | m(a) definition, 100×100 design | **Keep** | m(a) defined again at SI:830 and C:162 |
| 4.3 Simulation mechanics | 16–18 (2) | Parameter grid (S3), seed keys (S4), checkpoints, metric and SE formulas | none | Reproducibility | **Full** except Table S3 and one summary sentence | Overlaps Section 4 prose; "checkpoints/code" text is pipeline detail |
| 5 (heading) | 18 | — | — | — | — | — |
| 5.1 Purpose and observed inputs | 18–19 (1.5) | Purpose, 58 clusters, county conventions, numerator correction, rate formulas, Figures S2–S3 | C:180 | Data provenance | **Shorten**; home of trim T5 | Purpose duplicates C:147; data statement repeats SI:22 |
| 5.2 Geographic adaptations and budgets | 19–21 (2.5) | Adjacency check, each design's NC rule, region criterion C, blocking groups, Table S5, Table S6, Figures S4–S5 | C:165, C:226 | NC rules, budgets, region choice | **Keep** (rules, criterion, S5, S6); S4/S5 figures to Full | Region map = main Figure 1B; rounding rule wording differs from C:187 |
| 5.3 Outcome model and MC evaluation | 21–23 (2) | Education model (S3), estimation, interval, settings counts, replication schedule, m(a), tail refinement | C:185, C:193 | Precision targets; 1,000 outcomes; refinement | **Merge** with 5.4 and **Shorten**; home of trim T2 | Interval, one/two-sided power, year-sharing caveat all repeated |
| 5.4 Reproducibility and MC precision | 23–24 (1) | Refinement schedule, duplicate-allocation weights, raw/corrected risk variance, split-half check, Jaccard | C:193 | Split-half check | **Merge** into 5.3 | Data statement third copy (SI:837); split-half also at SI:2062 and C:193 |
| 5.5 Partial participation and sample size | 24–25 (1) | Subsets, volunteering models, W_P, exposure with non-participants, κ, MC sizes, ICC translation | C:193 | Power analysis design | **Keep**, split the 400-word paragraph | — |
| 5.6 Sensitivity analyses (methods) | 25–26 (1.5) | NC formal tests, CAR/Bayes models and priors, MCMC settings, spillover forms | none | Model-sensitivity design | **Shorten**: priors and MCMC settings to Full (or a compact table) | Test method duplicates 4.1 |
| 6 Estimation error and dependence on β | 26–27 (1.5) | Error decomposition (S4), g and h terms, IGS vs Isolation Buffer, rook Checkerboard coverage rise | C:206 | ρ̂ error scales with β | **Shorten** to ~0.75 p and **Merge** into results with Section 9 (it reports results, yet sits in Part I) | IGS control-only MSE 0.131 → 0.182 also at SI:1495; IB regime equivalence also at SI:1230 |
| 7.1 Design performance (queen) | 28–30 (2) | Pooled queen results, Table S7, Figures S6–S7, block wins | none | Pooled ranking | **Merge** with 7.4 (one nine-strategy table) | Table S12 repeats Table S7's MSE/coverage |
| 7.2 Rankings by incidence configuration | 30 (0.75) | Table S8, configuration stability | C:197 (no pointer) | "Configuration changed MSE < 10% except HIF" | **Shorten** to two sentences; Table S8 to Full | Main Figure 2 already shows configurations |
| 7.3 ρ, ψ and regime | 30–33 (2.5) | Tables S9–S11, regime reversal, Figure S8, IB regime equivalence | C:204 (Figure S8) | Regime drives ranking; loss is variance | **Shorten**: keep Figure S8 and the rook regime column; S9/S10 to Full | Queen regime column of S11 = main Table 2 |
| 7.4 Eight-design comparison | 33 (0.75) | Table S12 | E:15 | Variants adjacent to counterparts | **Merge** with 7.1 | Points back to 2.1; **stale** "four and three decimal places" (SI:1310) |
| 8 Formal statistical comparisons | 33–35 (2) | Friedman/Nemenyi/Wilcoxon (queen), rank vs mean reversal, 50-unit adjacent-pair tests, Figure S9 | none | Rankings are statistically separated | **Keep**, shorten the figure-description paragraph (SI:1409–1420) | Figure S9 caption repeats it; adjacent-pair tests omit SRS (review (c)) |
| 9 Stability across β | 35–38 (3) | Order across β, Friedman by β, margins, IB vs IGS, worst cases, coverage, Figures S10–S11 | C:206 (indirectly) | IGS advantage narrows as β grows | **Shorten** and **Merge** with 6 | "β levels share draws" caveat stated twice (SI:1464, SI:1536); SI:1528–1537 is a caption in prose |
| 10 Rook contiguity | 38–42 (4) | Rook ranking and tests, Checkerboard bias by ψ, Table S13, Figures S12–S14 | C:197 (no pointer) | Rook conclusions match queen; Checkerboard biased | **Shorten** to ~0.75 p; Figures S12–S14 and Table S13 to Full | Three "Figure X corresponds to Figure Y" paragraphs (SI:1642, 1658, 1673) |
| 11 Non-oracle sensitivity | 42–43 (1) | Omitted-spillover bias, coverage, lower MSE paradox, Table S14 | C:208 (no pointer) | Bias −0.121 to −0.400, coverage 54–91% | **Keep**; cut SI:1719–1721 | "these results do not change the recommendation" conflicts with the no-single-design rule |
| 12.1 Robustness and recommended design | 43 (0.75) | Table S15 block-level distribution, recommendation paragraph | none | IGS lowest worst-case block MSE | **Full** (Table S15); cut "recommended design" framing | Framing sentence third copy (SI:1817) |
| 12.2 Comparison with SRS | 43–45 (1.5) | Pooled vs SRS, Figure S15, Table S16, rook vs SRS, preview of NC tails | C:197 (main Table 2) | Ratios 0.63 / 1.08 | **Keep** Figure S15 and rook sentence; **Drop** Table S16 | Table S16 (label `tab:grid-srs`) duplicates main Table 2; SI:1839 previews Section 14 |
| 12.3 Allocation-specific risk (grid) | 45 (0.5) | Grid tail numbers, 63%/92% noise share | C:210 | Noise share of risk variance | **Keep** | Mean/worst numbers repeat C:210 and Table S19 |
| 13 Case study: accuracy, inference, reach | 45–53 (8) | Run totals, mean MSE by regime, coverage, budgets, targeting/rate share, Figures S16–S20, Tables S17–S18 | C:215, E:89 | Yearly MSE; coverage; reach; incidence share | **Shorten** to ~2.25 p | Four exhibits for one quantity (S16, S17, S18, S19); framing sentence second copy (SI:1904) |
| 14 Case study: allocation downside | 53–56 (3) | Tail reductions, Table S19, Figure S21, Figure S22, tail uncertainty, incidence-shift/rook/summary sensitivity, Table S20 | C:222, C:228 | Worst-decile reductions; MC noise share; sensitivity ratios | **Keep** text, **Drop** Figure S21 | Figure S21 is main Figure 3 (same file) |
| 15 (intro) | 56 (0.2) | Run provenance | — | — | **Shorten** to one sentence (pipeline jargon) | — |
| 15.1 Colleges needed and power | 56–58 (2.5) | Colleges for 80% power, MDE, ICC translation, unequal n_c, Wald size, Figures S23–S24, Tables S21–S22 | C:224 | 35–37 colleges; no design reaches 80% under control-only | **Keep** S21 and one figure; S22 and the second figure to Full | Calibration numbers repeat 15.2 |
| 15.2 Interval calibration | 58–59 (1) | Methods tried, coverage rule, Tables S23–S24 | C:158 | t + df chosen by calibration | **Keep** | Wald/t ranges also at SI:1915 (nine strategies) and C:215 (six + SRS), which differ |
| 15.3 Guidance under partial participation | 59–60 (1) | Table S25, saturation advantage vs N, incidence-weighted volunteering, region-size check | C:224, C:234 | 22% → 5% advantage | **Keep**, shorten the region-size check | — |
| 15.4 Sensitivity analyses (results) | 60–63 (3) | CAR/Bayes results, Tables S26–S28, spillover forms, NC formal tests | none | Comparison robust to model and spillover form | **Shorten** (split 450-word paragraph); Table S27 to Full | — |
| 15.5 Implications | 63 (0.3) | Planning advice | — | — | **Drop** (main Discussion covers it) | Repeats C:232 |

## Figures

| Fig | Label | Content | Main pointer | Approx. size | Recommendation | Notes |
|---|---|---|---|---|---|---|
| S1 | `fig:design-samples` | Example assignment per design, grid | none | 0.75 p | **Full** | Same surface and designs as main Table 1 thumbnails |
| S2 | `fig:service-areas` | County → college map with names | C:271 ("cluster maps") | 0.75 p | **Keep** | Only place the 100 → 58 mapping is shown |
| S3 | `fig:nc-observed` | Annual rates, four maps | none | 0.6 p | **Drop** | Main Figure 1C shows the same maps |
| S4 | `fig:nc-regions` | Four fixed regions | none | 0.4 p | **Drop** | Main Figure 1B |
| S5 | `fig:nc-examples` | Example NC assignment per strategy | C:271 ("an example assignment") | 0.6 p | **Full** (or Keep if C:271 keeps the phrase) | Panel titles numbered 1–9 (old design IDs) |
| S6 | `fig:mse` | Queen MSE by design × ρ × regime × ψ | none | 0.8 p | **Full** | Main Figure 2 covers configurations; text covers ρ/ψ |
| S7 | `fig:coverage` | Queen coverage boxes | none | 0.5 p | **Full** | One sentence suffices; possible in-figure "Queen Contiguity" title |
| S8 | `fig:biasvar` | Bias² + variance, queen, Poisson ρ_X = 0.2 | C:204 | 0.5 p | **Keep** | — |
| S9 | `fig:cd` | Nemenyi CD diagrams, queen + rook | none | 0.8 p | **Keep** | Caption and SI:1409–1420 overlap |
| S10 | `fig:tau-mse` | MSE across β, queen + rook | C:206 (indirect) | 0.5 p | **Keep** | Axis label carries a note ("band: ± ...") and minus-sign hyphens |
| S11 | `fig:tau-coverage` | Coverage across β | none | 0.8 p | **Full** | Text states the ranges |
| S12 | `fig:mse-rook` | Rook version of S6 | none | 0.8 p | **Full** | — |
| S13 | `fig:coverage-rook` | Rook version of S7 | none | 0.5 p | **Full** | — |
| S14 | `fig:biasvar-rook` | Rook version of S8 | none | 0.5 p | **Full** | — |
| S15 | `fig:vs-srs` | Ratio to SRS with t49 CIs, by regime, eight designs | (C:197 via Table 2) | 0.5 p | **Keep** | Adds the two variants to main Table 2's ratio column |
| S16 | `fig:nc-mse` | NC yearly MSE ratio to SRS, eight designs | none (same label as main Figure 4) | 0.5 p | **Drop** (or Full) | Main Figure 4 shows six; Table S18 has the numbers |
| S17 | `fig:nc-pop` | Population share by design/year | none | 0.6 p | **Full** | Table S5 has the same means |
| S18 | `fig:nc-setting-ratios` | Setting-level MSE ratios | none | 0.7 p | **Full** | In-figure "gamma=", "rho=" labels; teal/red palette differs from S16 |
| S19 | `fig:nc-mse-detail` | NC absolute MSE by design/year | none | 0.6 p | **Full** | Same numbers as Table S18 |
| S20 | `fig:spread-gain` | Guidance gain vs regional spread | none | 0.5 p | **Keep** | Supports the AGENTS targeting framing; natural partner for the IGS sensitivity study |
| S21 | `fig:alloc-dist` | Allocation-specific MSE distribution | C:210, C:222 (as main Figure 3) | 0.6 p | **Drop** | Identical file to main Figure 3, same label |
| S22 | `fig:nc-risk` | NC refined risk at ρ = 0.5 | none | 0.7 p | **Keep** | In-figure "Estimated q90" |
| S23 | `fig:colleges-needed` | Colleges for 80% power | none | 0.5 p | **Keep** | — |
| S24 | `fig:power-n` | Power vs N | none | 0.6 p | **Full** (keep one of S23/S24) | — |

## Tables

| Tab | Label | Content | Main pointer | Approx. size | Recommendation | Notes |
|---|---|---|---|---|---|---|
| S1 | `tab:designs` | Six design rules | C:165 | 0.5 p | **Merge** with S2 (one eight-design table) | Overlaps main Table 1 |
| S2 | `tab:excluded-designs` | Two variant rules | E:15 | 0.25 p | **Merge** into S1 | — |
| S3 | `tab:parameter-grid` | Parameter grid | none | 0.5 p | **Keep** | Caption still says "original eight-design simulation"; benchmark is nine with SRS |
| S4 | `tab:seed-keys` | Seed keys | none | 0.3 p | **Full** | — |
| S5 | `tab:nc-budgets` | NC budgets / population share | C:226 | 0.4 p | **Keep** | Counts printed as 27.00, 29.00 (integers need no decimals) |
| S6 | `tab:nc-order` | Regional ordering under three summaries | none (C:228 "reordered some regions") | 0.25 p | **Keep** | — |
| S7 | `tab:performance` | Six designs, queen + rook, MSE/SE/bias/coverage/power/rank | none | 0.6 p | **Merge** with S12 into one table with SRS and variant rows | Average rank among 7 here, among 9 in S12 |
| S8 | `tab:by-config` | MSE by configuration | none | 0.5 p | **Full** | Main Figure 2 |
| S9 | `tab:rho` | MSE by ρ | none | 0.3 p | **Full** | One sentence in text suffices |
| S10 | `tab:gamma` | MSE by ψ | none | 0.4 p | **Full** | Label still `tab:gamma` (harmless; rename for hygiene) |
| S11 | `tab:regime` | MSE by regime, queen + rook | none | 0.3 p | **Shorten** to the rook columns, or fold into merged S7 | Queen columns = main Table 2 |
| S12 | `tab:eight-design` | Eight designs + SRS, MSE/coverage/rank | E:15 | 0.5 p | **Merge** with S7 | Text claims four decimals; table shows three |
| S13 | `tab:checkerboard-coverage` | Checkerboard coverage by ρ | none | 0.2 p | **Full** | Text gives the range |
| S14 | `tab:oracle-nonoracle` | Oracle vs non-oracle | none (C:208) | 0.5 p | **Keep** | Add SRS row (rule: SRS in design comparisons) |
| S15 | `tab:robustness` | Block-level MSE distribution | none | 0.3 p | **Full** | Scenario-level worst cases, easily confused with allocation risk |
| S16 | `tab:grid-srs` | Nine-strategy regime MSE/coverage | none (same label as main Table 2) | 0.4 p | **Drop** | Main Table 2 + S15 + merged S7 cover it |
| S17 | `tab:nc-year` | NC yearly MSE, SRS + saturation | C:215 | 0.3 p | **Merge** into S18 and repoint C:215 | Subset of S18 |
| S18 | `tab:nc-annual` | NC MSE/bias/coverage by design, year, regime | E:89 | 1 p | **Keep** | Bias printed with hyphen-minus "-0.007" (use $-$) |
| S19 | `tab:downside` | Mean and worst-decile ratios, grid + NC | (C:210, C:222 via text) | 0.5 p | **Keep** | Ratios at 2 decimals (see review (c) on the MSE-precision rule) |
| S20 | `tab:nc-summary-sens` | IGS sensitivity ratios | C:228 (no pointer) | 0.3 p | **Keep** | Row label "Baseline / adjusted" is internal shorthand |
| S21 | `tab:nc-power` | Colleges needed, MDE | C:224 | 0.4 p | **Keep** | — |
| S22 | `tab:nc-translation` | Power by n_c, ICC, d | none | 0.5 p | **Full** (keep the worked example in text) | — |
| S23 | `tab:nc-calibration` | Size and coverage vs N | C:158 | 0.4 p | **Keep** | — |
| S24 | `tab:interval-comparison` | Wald vs t + df, grid + NC | C:158 (via 15.2) | 0.4 p | **Keep** (move next to S18) | AGENTS points to it by this label |
| S25 | `tab:nc-guidance` | MSE ratio vs N | C:224 | 0.4 p | **Keep** | — |
| S26 | `tab:nc-car` | SAR/CAR data × fit | none | 0.4 p | **Keep** | — |
| S27 | `tab:nc-bayes` | ML vs Bayes SAR vs Bayes CAR | none | 0.4 p | **Full** | Text carries the key numbers |
| S28 | `tab:nc-forms` | Spillover forms | none | 0.4 p | **Keep** | — |

## Cross-section duplications (consolidated)

1. Restricted-data statement: SI:22, SI:650, SI:837 (and main Data availability).
2. "Incidence guidance matches Saturation Quadrants' accuracy while directing saturation to high-incidence regions": SI:240–244, SI:1817, SI:1904; targeting numbers again at SI:1925.
3. Rook Checkerboard non-identification ("kept and flagged", "estimates β − ψ"): SI:314–322, SI:328–337, captions of Tables S7, S8, S10–S14 and Figures S10, S12–S14, prose at SI:1419, SI:1535, SI:1596, SI:1644.
4. t + df interval definition: SI:312, SI:567, SI:817, plus C:158.
5. One-sided grid power vs two-sided NC power: SI:407, SI:817, C:193.
6. Isolation Buffer "two regimes define the same model": SI:949–953, SI:1230–1235, C:204.
7. "Designs ignoring incidence share simulations across years": SI:827, SI:862, SI:1938, SI:1960, SI:1981, SI:1989, SI:2425 and Figure S16 caption.
8. Split-half worst-decile check: SI:848, SI:2062, C:193, C:222.
9. Worst-decile uncertainty (10 draws, Jaccard, noise share): SI:1887–1889, SI:2062, SI:2112–2114.
10. NC yearly MSE: Figure S16, Table S17, Figure S19, Table S18, main Figure 4.
11. Budgets and reach: Table S5, Figure S17, SI:1917–1923, C:226.
12. IGS vs Isolation Buffer control-only MSE across β (0.131 → 0.182): SI:955–957, SI:1494–1500, C:206.
13. "β values share assignments and noise": SI:1464–1467, SI:1536–1537, SI:397.
14. Main-text exhibits reproduced in the SI: Table 2 ≈ Table S16 (same label), Figure 3 = Figure S21 (same file and label), Figure 4 ≈ Figure S16 (same label), Figure 1B/C ≈ Figures S4/S3, Table 1 ≈ Table S1 + Figure S1.
