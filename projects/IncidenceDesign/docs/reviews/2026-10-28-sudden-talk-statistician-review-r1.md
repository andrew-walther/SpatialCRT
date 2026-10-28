# SUDDEN talk (2026-10-29): independent statistician review, round 1

**Date:** 2026-10-28
**Deck reviewed:** `presentation/sudden_talk_2026-10-29.qmd` (slides and `{.notes}` speaker script), with the figures in `presentation/figures/`
**Sources checked against:** `paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex`, `paper/thesis_chapter/Thesis_Chapter.qmd` (Methods, Results, Application, Appendices A8 and A9), `results/srs_benchmark/manuscript_regime_means.csv`, `results/srs_benchmark/manuscript_paired_ratios.csv`, `results/grid_tdf_2026-10/interval_comparison.csv`, `results/allocation_risk/summary/surfaces10/downside_by_design.csv`, `results/methods_ext_figures/tables/*.tex` and `tables_log.txt`, `application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv`, `observed_cluster_inputs.csv`, `sensitivity_setting_matched_comparison.csv`

**Verdict:** Blocking issues: 3 (plus 5 major, 14 minor)

The numbers hold up well. Of 52 spot checks, 52 match their sources (table below). The three blocking problems are all about definitions and framing, and each takes one line to fix. One symbol is defined incompletely (τ), one is defined wrongly (X in NC), and the take-home slide credits random-order saturation with a property only the incidence-guided version has.

---

## BLOCKING

### B1. τ is not stated as the effect *before* spatial propagation
- **Slide:** "One outcome model links treatment, spillover and neighbors" (glossary). The same reading carries into "A one-SD college effect is a modest effect for each trained professional" and "About 36 colleges detect a one-SD effect...".
- **Text:** "$\tau$ — Direct effect of a cluster's own treatment, in standard deviations of the noise $\epsilon$ (the quantity we estimate)". Notes: "Tau times Z is the direct effect of the cluster's own treatment."
- **Problem:** Both sources define τ as the coefficient inside the brackets. It is the effect *before* outcomes spread through $(I-\rho W)^{-1}$, not a total effect. CTJ: "the effect of a cluster's own treatment before outcomes spread between neighbors through $(I-\rho W)^{-1}$, rather than a total effect". Chapter, line 506: "τ is the coefficient before spatial propagation, not an average total impact". The power slides use ρ = 0.5, where the realized own-cluster effect is larger than τ, and the effect-size slide reads "a college effect of τ SDs" as a plain shift in the college score. Someone in the audience could fairly ask whether "1 SD" means the effect on the score or a model parameter. For the prelim audience this is the first question.
- **Fix:** Glossary entry: "Direct effect of a cluster's own treatment *before* outcomes spread between neighbors through $(I-\rho W)^{-1}$ (not a total effect), in SDs of the noise $\epsilon$". Notes: add one sentence: "Tau is the push the program gives a college's own score before that score spreads to neighbors, so it is not the total effect." On the effect-size slide, add a fig-note: "College effect = τ, the pre-propagation coefficient, in noise-SD units."

### B2. The X definition is wrong for North Carolina
- **Slide:** Outcome model glossary.
- **Text:** "$X$ — Cluster incidence on a 0 to 1 scale (simulated on the grid, observed SUD rates in NC)".
- **Problem:** Observed SUD rates run from about 45 to 195 per 100,000. They are not on a 0–1 scale. In NC, X is in the model only in the β = 1 sensitivity, and there it is the rank-scaled rate $X = (\text{rank}-0.5)/58$ (chapter line 508). The grid Poisson configurations are also rank-normalized. The β = 1 backup already says "incidence rank", which contradicts the glossary.
- **Fix:** "$X$ — Cluster incidence scaled to 0–1 (simulated on the grid; in NC, the rank of the observed SUD rate, (rank − 0.5)/58, used only in the β = 1 sensitivity)".

### B3. Take-home slide credits random-order saturation with burden targeting
- **Slide:** "Choose the design for the spillover you expect", left card, and its notes.
- **Text:** "**Saturation** (Incidence-Guided or random-order): lowest or near-lowest error, a lower risk of an unlucky allocation, and the largest treated share in the highest-burden region". Notes: "saturation designs ... still gave the highest-burden region the largest share of treatment."
- **Problem:** Saturation Quadrants assigns 80/60/40/20% in random order, so it does not give the highest-burden region the largest share. Only Incidence-Guided Saturation does (CTJ conclusions credit this property to Incidence-Guided Saturation alone). The deck's own targeting slide shows 50.4% for random order.
- **Fix:** "**Saturation** (Incidence-Guided or random-order): lowest or near-lowest error and a lower risk of an unlucky allocation; the incidence-guided version also gives the highest-burden region the largest treated share". Make the same change in the notes.

---

## MAJOR

### M1. "Fine in one regime" contradicts the table and the take-home
- **Slide:** "On the grid, the spillover regime decides which design wins", amber card. Notes: "They're tolerable when spillover reaches both arms."
- **Text:** "Checkerboard and 2×2 Blocking are fine in one regime and poor in the other"
- **Problem:** Under both arms both designs have 1.47 times SRS's MSE (CI 1.44–1.51), shown orange in the same table. The take-home says they "did not improve on SRS in either regime", and the chapter Results agree. (The chapter Discussion's "close to the best designs" compares them with the other designs, not with SRS.)
- **Fix:** "Checkerboard and 2×2 Blocking were about 1.5 times SRS's error under both arms and far worse under controls only (2.7 and 9.8 times)."

### M2. Rook backup overstates agreement for "every other design"
- **Slide:** "Backup · Under rook neighbors, Checkerboard cannot estimate the effect"
- **Text:** "Every other design: same broad conclusions under rook; Incidence-Guided Saturation and Balanced Quartiles lead (MSE 0.072, 0.084)"
- **Problem:** Under rook, Isolation Buffer goes from 26% *lower* control-only MSE than SRS (queen) to 5.8% *higher* (chapter line 372). The saturation advantage also shrinks from 37% to 23–24%. AGENTS.md (Invariant 2) asks for the split to be shown wherever a conclusion changes.
- **Fix:** "Under rook the saturation advantage under controls only is smaller (23–24% lower MSE than SRS, against 37% under queen), and Isolation Buffer loses its advantage (5.8% higher than SRS)."

### M3. Limitations fig-note and notes overstate robustness to other models
- **Slide:** "What these results do not tell us yet", fig-note. Notes: "The results held with other spatial models, including Bayesian ones."
- **Text:** "The design ordering held under conditional autoregressive and Bayesian models and with rook neighbors (backup slides)."
- **Problem:** The deck's own backup has an exception. A Bayesian CAR (Leroux) fit to strongly dependent spatial lag data was biased by 0.14–0.22, covered 88–90%, and lost the saturation advantage. Rook also changes Isolation Buffer (M2).
- **Fix:** "The regime-specific guidance held under CAR models, a Bayesian spatial lag model and rook neighbors, except a misspecified Bayesian CAR fit (backup slides)." Notes: "...including a Bayesian spatial lag model; a Bayesian CAR model that assumes a different dependence structure was the exception."

### M4. Sample-size claims lose their conditions on the title and take-home
- **Slides:** "About 36 colleges detect a one-SD effect if spillover reaches every college" (title). Take-home: "about **36 colleges** detect a one-SD college effect".
- **Problem:** The 35–37 figure holds only for 80% power, the corrected t test, random volunteering, ρ = 0.5 and γ = 0.8τ, with τ in noise-SD units before propagation (B1). The caption says this, but the title and the take-home do not. The audience will remember "36 colleges". The glossary also defines γ as the fixed range 0.5–0.8, while these slides use γ = 0.8τ, which is called only "spillover proportional to the effect".
- **Fix:** Title: "About 36 colleges give 80% power for a one-SD effect if spillover reaches every college". Caption: "spillover γ = 0.8τ". Take-home: "about 36 colleges give 80% power for a one-SD college effect (ρ = 0.5, random volunteering)".

### M5. The allocation-risk title is unconditional
- **Slide title:** "A trial gets one allocation, and saturation also lowers the risk of an unlucky one"
- **Problem:** The claim holds only under controls-only spillover. Under both arms, Incidence-Guided Saturation's grid worst-decile mean is 1.07 times SRS's (0.0583 vs 0.0544), and the NC differences are mostly noise. The bullets carry this qualification but the title does not, and the title is what the audience takes away.
- **Fix:** "A trial gets one allocation; under controls-only spillover, saturation also lowers the risk of an unlucky one"

---

## MINOR

1. **Style: colons in prose** (AGENTS/memory rule). Slide lines 193 ("is classic contamination: untreated..."), 222 ("Queen contiguity: two clusters..."), 224, 227 (glossary "Spatial dependence: ...", "Spillover magnitude: ..."), 623 (take-home card), 781–782 (rook backup), 837, 964, and notes line 672 ("together: define..."). Replace each with a comma, a dash or a sentence break. Bold labels such as "**Controls only:**" in bullets are acceptable.
2. **Figure label:** `fig_colleges_needed.png` prints "2x2 Blocking" with a letter x. Every other exhibit uses "2×2 Blocking". Regenerate with the × sign.
3. **"SD" is never defined.** It appears in titles from the "Colleges needed" slide onward. Write "standard deviation (SD)" at first use, on the "Colleges needed" fig-note.
4. **ICC is defined only in the notes.** The effect-size slide says "intraclass correlation ICC" but not what it measures. Add "(the share of score variance between colleges)".
5. **The "Strict alternation (Checkerboard, 2×2 Blocking)" label** on the take-home is inaccurate for 2×2 Blocking, which is blocking, not alternation (and becomes spatial blocking in NC). Use "Checkerboard and 2×2 Blocking did not improve on SRS in either regime."
6. **"With no loss of accuracy"** (targeting slide) and "costs nothing in accuracy" (notes). The author's framing rule (AGENTS.md, 2026-10-21) is that guidance *matches Saturation Quadrants' accuracy*. Under both arms, Incidence-Guided Saturation is 8.3% above SRS and Saturation Quadrants 4.3% above. Use "matching the accuracy of random-order saturation".
7. **Background slide:** "the treated-versus-control gap understates the effect" holds for spillover to controls. Under both-arms spillover the direction depends on the allocation. Add "when spillover reaches untreated clusters" or leave it as is for a clinical audience.
8. **Design slide notes:** "Balanced Quartiles ... so the arms have the same incidence". Say "similar incidence distributions".
9. **Grid-study notes:** "The overall ordering of designs was the same at every effect size we tried, with one exception at large effects." The pooled order had no exception (chapter line 315). The exception is within the controls-only regime. Say "Pooled over regimes the ordering was the same at every effect size; under controls only, Isolation Buffer caught up at large effects."
10. **Shared simulations across years:** the NC accuracy slide says saturation was better "in every year". SRS, Checkerboard, Isolation Buffer and Saturation Quadrants reuse identical simulations across years (identical cells in the heat map), so the four years are not independent replications. Add a fig-note: "Years share one map and some designs share simulations across years; not independent replications."
11. **22% vs 21–30%:** the partial-participation table gives "58 → 22% lower" while the NC accuracy slide gives "21–30% lower" at 58 colleges. These come from different settings (ρ ∈ {0, 0.5}, γ = κτ vs the six primary settings). Label the table "power-analysis settings (ρ = 0, 0.5; γ = 0.5τ, 0.8τ)" or expect the question.
12. **Formal-tests backup:** "Incidence-Guided Saturation ranked best, significantly ahead of every other design and SRS" is pooled over both regimes. Add "pooled over spillover regimes" so it does not read as a single recommended design (author rule 2026-10-14).
13. **Rook backup wording:** "treatment and exposure are the same variable". Precisely, exposure is exactly 1 − Z, so it is collinear with treatment and the intercept. Use "exposure equals 1 − Z, so it cannot be separated from treatment".
14. **Timing:** the notes' running total reaches 20.5 min before questions for a 20-minute slot. Trim about 0.5 min (for example, the design-gap slide or the budgets slide notes).

Checked and clean: no DIM mention, no Claude/Codex mention, no sentence starts with "Because", design names correct throughout the slide text (the Saturation Quadrants name is used, not the CSV's "Saturation Regions"), β values correct (grid β = 1; NC β = 0 primary, β = 1 hypothetical), both spillover forms correct and matching the chapter/CTJ equations, the interval formula matches chapter line 137 and CTJ, the m(a) definition matches CTJ, rook Checkerboard identification (WZ = 1 − Z, estimate centers on τ − γ) correct, no causal SUD-reduction claim (explicitly disclaimed), no single recommended design.

---

## Number spot-checks (52)

| # | Number | Slide | Source | Match |
|--:|---|---|---|:-:|
| 1 | IGS/SRS controls only 0.63 | Grid results | manuscript_paired_ratios.csv 0.631 | yes |
| 2 | IGS/SRS both 1.08 | Grid results | 1.083 | yes |
| 3 | Isolation Buffer 0.74 / 3.80 | Grid results | 0.743 / 3.796 | yes |
| 4 | Balanced Quartiles 1.02 / 0.99 | Grid results | 1.023 / 0.995 | yes |
| 5 | High Incidence Focus 1.53 / 3.82 | Grid results | 1.534 / 3.818 | yes |
| 6 | 2×2 Blocking 2.68 / 1.47 | Grid results | 2.685 / 1.471 | yes |
| 7 | Checkerboard 9.84 / 1.47 | Grid results | 9.839 / 1.473 | yes |
| 8 | SRS MSE 0.215 / 0.042 | Grid results | regime_means 0.2146 / 0.0424 | yes |
| 9 | 37% lower, 0.136 vs 0.215 | Grid results | 0.1355 vs 0.2146; chapter 36.9% | yes |
| 10 | BQ grey (CI includes 1) | Grid results | 0.995–1.050 | yes |
| 11 | Coverage 94–95% | Grid results | regime means 0.940–0.952 | yes |
| 12 | Backup Table 2 MSE/coverage, all 14 cells | Grid benchmark backup | regime_means.csv / chapter tab:grid-srs | yes |
| 13 | Ratio CIs e.g. 9.84 (9.65, 10.03), 0.63 (0.61, 0.65) | Grid benchmark backup | paired_ratios.csv | yes |
| 14 | Pooled 0.0907 vs 0.1285 | Grid backup notes | chapter line 370 | yes |
| 15 | 23,523 deaths, 25.6M person-years | NC incidence | observed_cluster_inputs.csv 23,523 / 25,594,321 | yes |
| 16 | Rates about 45 to 195 | NC incidence notes | pooled min 45.3 (Wake Tech), max 194.9 (Edgecombe) | yes |
| 17 | 29 of 58; Checkerboard 27; blocking 27–31; Isolation Buffer about 16 | NC allocations | yearly means: 27; 28/31/27/29; 15.89 | yes |
| 18 | IGS reach 53/45/48/52% | Budgets | 53.2/45.1/47.9/51.8 | yes |
| 19 | HIF reach 22/24/23/22% | Budgets | 22.2/23.5/23.3/22.4 | yes |
| 20 | Isolation Buffer reach 27% | Budgets | 26.7 | yes |
| 21 | BQ 49/50/51/50; Checkerboard 50/50/50/51 | Budgets | 48.5/49.7/50.6/49.8; 50.4/50.4/50.4/50.7 | yes |
| 22 | IGS controls only 21–30% lower | NC accuracy | ratios 0.789/0.701/0.716/0.771; chapter 21.1–29.9% | yes |
| 23 | HIF 0.70 to 1.40 | NC accuracy | 0.697–1.398 | yes |
| 24 | Both arms SRS/BQ/IGS within ~3% | NC accuracy | BQ 0.976–0.998, IGS 0.972–1.006 | yes |
| 25 | Checkerboard ~3×, Isolation Buffer 3.7× | NC accuracy | 2.83/2.99; 3.74 | yes |
| 26 | Worst 10% 21–43% lower (NC) | Allocation risk | chapter/CTJ 20.7–43.4% | yes |
| 27 | Grid worst decile 0.19 vs 0.35 | Allocation risk | downside csv 0.1912 vs 0.3462 | yes |
| 28 | Allocation backup table, all 64 cells | Allocation risk backup | downside_by_design.csv / chapter tab:downside | yes |
| 29 | 35–37 colleges both arms | Colleges needed | tab_nc_power.tex 35–37 | yes |
| 30 | 28 vs 32 at 1.5 SD controls only | Colleges needed | tab_nc_power 28 / 32 | yes |
| 31 | 22%, 9–10%, 4–5%, ~1% | Partial participation | tab_nc_guidance (0.78; 0.90–0.91; 0.95–0.96; 0.99); chapter line 656 | yes |
| 32 | ≤2% with incidence-weighted volunteering | Partial participation | tables_log IGS 1.00/0.99 | yes |
| 33 | σ_c ≈ 0.30 (m = 25, ICC = 0.05) | Effect size | √(0.05 + 0.95/25) = 0.297 | yes |
| 34 | Power 0.51/0.82, 0.71/0.96, 0.80/0.98 | Effect size | tab_nc_translation.tex | yes |
| 35 | Controls only 0.49/0.53 SRS, 0.51/0.64 IGS | Effect size | chapter line 652 | yes |
| 36 | 62.2/35.4, 51.1/49.2, 50.4/49.6 | Targeting | chapter line 538 | yes |
| 37 | Mean ranks 24–34; grid gain 2–3 pts | Targeting | chapter 24–34; 2.1–3.1 | yes |
| 38 | Grid spread correlation 0.84–0.99 | Spread backup notes | chapter line 538 | yes |
| 39 | Wald/t coverage, all 42 cells | Interval backup | interval_comparison.csv | yes |
| 40 | Test size 11.1–12.8 / 4.4–5.4 (15), 8.3–8.8 / 5.0–5.5 (29), 6.3–7.0 / 5.0–5.5 (58) | Interval backup | tab_nc_calibration.tex | yes |
| 41 | Corrected coverage 94.0–95.6% | Interval backup | chapter line 654 | yes |
| 42 | Translation table, all 81 cells | Power translation backup | tab_nc_translation.tex | yes |
| 43 | CAR 0.76–0.82 | CAR/Bayes backup | tab_nc_car.tex controls-only columns | yes |
| 44 | Bayes r = 0.999, coverage 94.6–95.3; Leroux bias 0.14–0.22, coverage 88–90% | CAR/Bayes backup | chapter line 658 | yes |
| 45 | Spillover-forms table, all 36 cells; 4.93; 90–92%; 0.70–0.71 | Spillover forms backup | tab_nc_forms.tex; chapter line 1630 | yes |
| 46 | Friedman χ² 560.3, df 6, 160 blocks; Wilcoxon p ≤ 4.9e-9 | Formal tests backup | chapter line 214 | yes |
| 47 | NC Friedman 135–180, df 8, 24 blocks | Formal tests backup | tables_log 135.3–179.8 | yes |
| 48 | HIF 3.7–4.4× and 1.4–2.5× (β = 1); randomized designs little affected | β = 1 backup | sensitivity csv 3.68–4.44, 1.42–2.52; others 0.95–1.08 | yes |
| 49 | Omitted spillover bias −0.12 to −0.40, coverage 54–91% | Omitted spillover backup | chapter line 311 | yes |
| 50 | τ = 2: 0.160 / 0.160; τ = 3: 0.160 vs 0.182 | Larger effects backup | chapter line 333 (0.1596/0.1604; 0.160/0.182) | yes |
| 51 | Rook IGS 0.072, BQ 0.084; Checkerboard coverage 0.16 | Rook backup | chapter line 305; interval_comparison 0.1615 | yes |
| 52 | Isolation Buffer about 21 (queen grid); SRS population 31–69% | Design rules / population backups | chapter 21; yearly means 0.313–0.687 | yes |

The six figures I opened (`nc_annual_accuracy_six.png`, `allocation_distribution_by_design.png`, `nc_design_examples.png`, `fig_colleges_needed.png`, `fig_power_vs_n.png`, `nc_planning_map.png`) agree with the slide text and the speaker-note descriptions of panel positions.
