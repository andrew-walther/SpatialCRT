# SUDDEN talk (2026-10-29): independent statistician review, round 2

**Date:** 2026-10-28
**Deck reviewed:** `presentation/sudden_talk_2026-10-29.qmd` (slides and `{.notes}` script, 928 lines), figures in `presentation/figures/` (opened `fig_colleges_needed.png`, `fig_power_vs_n.png`, `design_assignments_panels_wide.png`, `nc_design_examples.png`, `allocation_distribution_by_design.png`, `nc_annual_accuracy_six.png`)
**Sources checked against:** `paper/thesis_chapter/Thesis_Chapter.qmd` (lines 214, 305, 315, 372, 480–492, 508, 518, 575–587, 630–658, 1421–1490), `paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex` (151, 167, 187, 198, 229, 235, 272), `results/methods_ext_figures/tables/tab_nc_power.tex`, `tab_nc_car.tex`, `tables_log.txt`, `results/allocation_risk/summary/surfaces10/downside_by_design.csv`, `results/grid_tdf_2026-10/interval_comparison.csv`, `application/results/real_sud_rev_20261002/exhibits/sensitivity_setting_matched_comparison.csv` (NC rook ratios recomputed), `DECISIONS.md` (IRB)

**Verdict:** Blocking issues: 0 (plus 1 major, 11 minor)

The revision fixed all three round-1 blocking findings and four of the five major ones. The remaining major issue is that the colleges-needed title now drops the effect size and the power target altogether ("suffice" for what?), and the take-home still says "detect" without "80% power". Every number I checked in the rewritten wording traces to a source. The minor items are wording precision and one figure label.

---

## Round-1 findings: status

| Finding | Status | Evidence in the revised deck |
|---|---|---|
| B1 τ not stated as pre-propagation | **Fixed** | Glossary line 213 ("before outcomes spread between neighbors through $(I-\rho W)^{-1}$ (not a total effect)"); notes line 229. The suggested fig-note on the effect-size slide was not added (see m4). |
| B2 X wrong for NC | **Fixed** | Line 219: "in NC, the rank of the observed SUD rate scaled to 0 to 1, used only when β = 1". Matches chapter line 508. |
| B3 take-home credits random-order saturation with targeting | **Fixed** | Line 579 credits the highest-rate region only to Incidence-Guided Saturation; notes line 596 likewise. |
| M1 "fine in one regime" | **Fixed** | Line 333: "worse than SRS in both regimes, and far worse under controls only"; notes line 346. |
| M2 rook backup overstates agreement | **Fixed** | Line 738 gives 23–24% vs 37% and Isolation Buffer 5.8% higher (chapter line 372). |
| M3 robustness overstated | **Fixed** | Fig-note line 621 names the misspecified Bayesian CAR exception; the notes no longer make the claim. |
| M4 sample-size claims lose their conditions | **Partly fixed** | Caption (line 467) and γ glossary now carry the conditions, but the title lost even more (now "About 36 participating colleges suffice...") and the take-home still omits 80% power. See MAJOR M1 below. |
| M5 allocation-risk title unconditional | **Fixed** | Title line 440: "Under spillover to controls only, saturation also lowers the risk...". |

Round-1 minors: fixed 1 (mostly; one prose colon remains, m9), 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14. **Not fixed:** 2 (`fig_colleges_needed.png` still prints "2x2 Blocking"). Item 7 left as is, which is acceptable for this audience.

---

## BLOCKING

None.

---

## MAJOR

### M1. The colleges-needed title and the take-home drop the power target (and the title drops the effect size)
- **Slides:** "Colleges needed" (title, line 463) and "Take-home" (navy card, line 584).
- **Text:** Title: "About 36 participating colleges suffice if spillover reaches every college". Take-home: "about **36 participating colleges** detect a one-SD college effect (about 0.3 SD per professional with 25 trained per college), under our assumptions".
- **Problem:** "Suffice" has no object. The title no longer says for what effect or what power, which is what an audience remembers and what a clinician would quote. The number 35–37 holds only for 80% power at τ = 1 college-level SD under both-arms spillover, ρ = 0.5, γ = 0.8τ and random volunteering (chapter line 634, `tab_nc_power.tex`). "Detect" on the take-home reads as certainty, not an 80% chance. At 1.5 SD the same designs need about 19, so the effect size changes the answer by half.
- **Fix:** Title: "About 36 participating colleges give 80% power for a one-SD effect if spillover reaches every college". Take-home: "about **36 participating colleges** give 80% power for a one-SD college effect (about 0.3 SD per professional with 25 trained per college), under our assumptions". Notes line 596: "...about 36 participating colleges give an 80 percent chance of detecting a one-SD college effect under our assumptions."

---

## MINOR

1. **Figure label (round-1 minor 2, not fixed).** `presentation/figures/fig_colleges_needed.png` prints "2x2 Blocking" with a letter x; every other exhibit and the deck text use "2×2 Blocking". Regenerate the figure (`code/23_methods_extension_figures.R`) with the × sign and refresh the presentation copy.

2. **"Power stalls past about 40 colleges"** (Partial participation slide, line 500). True for SRS (and Checkerboard falls), but in `fig_power_vs_n.png` both saturation designs keep rising, from about 0.54 at 35 colleges to 0.63 at 58 (1 SD). The chapter says power "rose slowly beyond about 40 colleges" (line 634). Fix: "Under controls only, power rises slowly past about 40 colleges, and SRS stalls, since nearly every control then has treated neighbors."

3. **High Incidence Focus under both arms in NC is omitted from the accuracy bullets.** The heat map the audience sees shows High Incidence Focus at 0.92–0.95 of SRS under both arms in every year, the lowest row on the left (and first in the NC Friedman ranks, `tables_log.txt`). The bullets (lines 421–426) and the take-home ("High Incidence Focus varied by year", line 588) describe only its controls-only swing. Someone will point at the 0.94s. Fix: add to the left bullet "High Incidence Focus was 5–8% below SRS, but that rests on incidence not affecting scores (β = 1 backup)"; take-home "High Incidence Focus varied by year under controls only and is not randomized".

4. **The effect-size slide equates the college-score SD with the noise SD** (line 519–521: "its standard deviation is σ_c ... a college effect of τ SDs"). In the model τ is in SDs of ε, before propagation (glossary, B1). At ρ = 0.5 the observed college score has a larger SD than ε because of $(I-\rho W)^{-1}$. The chapter (line 518) maps σ_c to the noise, which is fine, but the slide's wording invites a prelim question. Fix: "its noise standard deviation (the SD of $\epsilon$ in the model) is σ_c" and "a direct college effect of τ noise SDs is an individual effect of $d = \tau\sigma_c$".

5. **Conditional conversion stated as general** (Colleges needed, line 469). "one college-level SD ≈ 0.3 SD per professional" holds only for m = 25, ICC = 0.05. Fix: "one college-level SD ≈ 0.3 SD per professional with 25 trained per college and ICC 0.05".

6. **NC accuracy notes understate the β = 1 increase** (line 437). "its error rose up to fourfold": the source range is 3.68–4.44 under both arms (chapter line 630; CTJ line 235 "roughly fourfold"). Fix: "its error rose roughly fourfold".

7. **Allocation-risk bullet is ambiguous about which designs** (line 450). "Both arms: differences in the worst allocations are small and mostly simulation noise" is true for saturation, Balanced Quartiles and SRS (0.93–1.07), but the figure shows Isolation Buffer's worst decile at 3.97 (grid) and 4.75 (NC) times SRS's. Fix: "Both arms: saturation's worst allocations are close to SRS's (0.93–1.07 of SRS), and the differences are mostly simulation noise".

8. **Formal-tests backup misstates the unit-level result** (line 900). "Adjacent designs in the MSE order differed in all 50 independent configuration-by-surface units" is wrong as written. Chapter line 214 says every adjacent pair differed *significantly* and the better design had lower MSE in all 50 units *except* High Incidence Focus vs 2×2 Blocking (48 of 50). Fix: "Each pair of designs adjacent in the MSE order differed significantly over the 50 independent configuration-by-surface units (Wilcoxon $p \le 4.9\times10^{-9}$); the better design was lower in all 50 units for every pair but one (48 of 50)".

9. **Style: one prose colon and a double "and"** (Rook backup, line 736). "Under rook, Checkerboard gives $WZ = 1 - Z$: exposure equals $1 - Z$, so it cannot be separated from treatment and $\tau$ is **not identified** and the estimate centers on $\tau - \gamma$". Fix: "Under rook, Checkerboard gives $WZ = 1 - Z$, so exposure cannot be separated from treatment. $\tau$ is **not identified**, and the estimate centers on $\tau - \gamma$ (coverage 0.16)".

10. **β = 1 backup grammar** (line 781). "adds incidence rank to the outcome ($\beta = 1$) and adjusted for in the analysis" mixes tenses. Fix: "adds incidence rank to the outcome ($\beta = 1$) and adjusts for it in the analysis". Its notes (line 788) say High Incidence Focus "is not recommended"; that is consistent with the no-single-recommended-design rule but reads better as "is a risky choice".

11. **Take-home bullet phrasing** (line 588). "Checkerboard never improved on SRS, and 2×2 Blocking was no better than SRS" says the same thing twice in different words. Fix: "Checkerboard and 2×2 Blocking never improved on SRS; High Incidence Focus varied by year under controls only and is not randomized" (combines with m3).

---

## Accuracy, notation and framing: checked and clean

- **Notation.** Every symbol in the outcome equation is defined (Y, Z, W, τ, ρ, S(Z), $(I-\rho W)^{-1}$, γ, β, X, ε, ∘). Both spillover forms match CTJ/chapter. m(a) is defined as an allocation's expected squared error (matches CTJ). σ_c formula and ICC definition match chapter line 518. SD is defined at first use (line 467). The interval formula matches chapter line 1470 and the SI.
- **β values.** Grid β = 1; NC β = 0 primary, β = 1 hypothetical only (glossary line 218, notes line 231, backup lines 779–788). Correct.
- **γ.** Glossary now distinguishes the grid (0.5–0.8) from the power analysis (0.5τ or 0.8τ). The NC accuracy slide's γ ∈ {0.5, 0.8} lies within the grid range stated.
- **Overclaiming.** "No single design is best; the spillover regime decides" (line 573). No causal SUD claim; line 610 and notes line 626 disclaim it. Shared simulations across years and blocks are flagged (lines 417, 902). Formal tests are labeled "pooled over spillover regimes" (line 899).
- **Style.** No DIM, Claude or Codex mention; no sentence begins with "Because"; design names match AGENTS.md everywhere in the deck text (the only "saturation regions" is the alt text for the map panel describing the four regions, which is fine). Prose colons are gone except line 736 (m9); label colons in bullets and fig-notes are acceptable.
- **IRB 14-2036** on the thank-you slide is the author-supplied protocol (`DECISIONS.md` line 146; CTJ line 272), not invented.
- **Timing.** Main body runs to 18.5 min before questions.

---

## Number spot-checks this round (34)

| # | Number (deck) | Slide | Source | Match |
|--:|---|---|---|:-:|
| 1 | 35–37 colleges, both arms 1 SD (HIF 35, 2×2 37, rest 36) | Colleges needed | tab_nc_power.tex col 1 | yes |
| 2 | Checkerboard and Isolation Buffer >58 at 1 SD both arms | Colleges needed notes | tab_nc_power.tex; figure | yes |
| 3 | No design 80% at 1 SD controls only, even with 58 | Colleges needed | MDE at 58 under controls 1.24–>2 (tab_nc_power); chapter 634 | yes |
| 4 | Saturation 28 vs SRS 32 at 1.5 SD | Colleges needed | tab_nc_power 28/28/32; figure | yes |
| 5 | Grey = not reached in at least one year | Colleges needed fig-note | figure legend; chapter 1421 | yes |
| 6 | 22% / 9–10% / 4–5% / ~1% lower | Partial participation | chapter 656 | yes |
| 7 | ≤2% at 20–29 with incidence-weighted volunteering | Partial participation | tables_log IGS 1.00/0.99, SQ 1.01/0.98; chapter 656 | yes |
| 8 | SRS and Checkerboard lose power near 58 (controls only) | Partial participation notes | fig_power_vs_n.png (SRS 0.54→0.52; Checkerboard 0.30→0.17) | yes |
| 9 | Power stalls past ~40 | Partial participation | figure: SRS yes, saturation rises (m2) | partly |
| 10 | σ_c ≈ 0.30 at m = 25, ICC = 0.05 | Effect size | √(0.05+0.95/25) = 0.297 | yes |
| 11 | Power 0.51/0.82, 0.71/0.96, 0.80/0.98 | Effect size | backup table rows m = 10/25/50, ICC 0.05; chapter 652 | yes |
| 12 | All nine σ_c values in the translation backup | Power translation backup | recomputed from formula (0.330, 0.381, 0.436, 0.223, 0.297, 0.369, 0.173, 0.263, 0.344) | yes |
| 13 | "Lowered power by about 0.1" (unequal m) | Power translation backup | tables_log A4 −0.097 | yes |
| 14 | Rook saturation 23–24% lower, Isolation Buffer 5.8% higher | Rook backup | chapter 372 | yes |
| 15 | Rook pooled IGS 0.072, BQ 0.084; coverage 0.16 | Rook backup | chapter 305; interval_comparison 0.1615 | yes |
| 16 | NC rook gives the same broad conclusions | Rook backup notes | recomputed from sensitivity_setting_matched_comparison.csv: controls only SQ 0.73, IGS 0.74, Checkerboard 3.76; both arms Checkerboard 3.27, Isolation Buffer 3.90, others 0.98–1.03 | yes |
| 17 | Worst 10% 21–43% lower in NC | Allocation risk | chapter 575 (20.7–43.4%); CTJ 229 | yes |
| 18 | Grid worst decile 0.19 vs 0.35 | Allocation risk | downside csv 0.1912 / 0.3462 | yes |
| 19 | Allocation backup cells (spot: Checkerboard grid 1.41/1.11/9.94/6.22; Isolation Buffer NC 3.69/4.75/0.93/0.91; IGS NC 0.98/0.93/0.75/0.68) | Allocation backup | downside_by_design.csv | yes |
| 20 | Grid units "50 ... at ρ ∈ {0, 0.5}, γ ∈ {0.5, 0.8}" | Allocation backup fig-note | chapter caption, 200 units = 50 × 4 | yes |
| 21 | Orange bar = SRS; saturation bar left of 1 | Allocation risk notes | allocation_distribution_by_design.png | yes |
| 22 | IGS 21–30% lower, bottom row; HIF 0.70–1.40; Checkerboard ~3×; Isolation Buffer 3.74 | NC accuracy | nc_annual_accuracy_six.png | yes |
| 23 | Both arms SRS/BQ/IGS within ~3% | NC accuracy | heat map 0.97–1.01 | yes |
| 24 | β = 1 HIF "up to fourfold" | NC accuracy notes | chapter 630, 3.68–4.44 (m6) | partly |
| 25 | 2×2 Blocking reach 50/55/48/50, counts 28/31/27/29 | Budgets | chapter 487 | yes |
| 26 | SRS reach 50% | Budgets | chapter 49.5% (rounds to 50) | yes |
| 27 | Isolation Buffer 14 in 2018 example, ~16 average; HIF top right | NC allocations | nc_design_examples.png; chapter 487 (15.89) | yes |
| 28 | Percentages only on IGS panel; quartile numbers on BQ panel; SRS bottom right | Designs | design_assignments_panels_wide.png | yes |
| 29 | Saturation 0.76–0.82 under CAR; Checkerboard and Isolation Buffer worst both arms | CAR backup | tab_nc_car.tex (IGS/SQ 0.76–0.82; Checkerboard 2.83–2.93, Isolation Buffer 4.11–4.45) | yes |
| 30 | Grid Friedman 560.3, IGS significantly best, SRS ≈ BQ | Formal tests backup | chapter 214 | yes |
| 31 | "Differed in all 50 units" | Formal tests backup | chapter 214 (48 of 50 for one pair) (m8) | no |
| 32 | NC Friedman 135–180, df 8; saturation two best ranks controls only | Formal tests backup | tables_log (SQ 1.67/1.29, IGS 1.79/1.71) | yes |
| 33 | Bootstrap and profile intervals did not meet the 94–96% rule | Interval backup notes | chapter 1470 | yes |
| 34 | Wald / t coverage, IGS and Checkerboard rows | Interval backup | interval_comparison.csv (0.9366/0.9469 etc.) | yes |
