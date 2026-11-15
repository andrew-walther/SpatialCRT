# Chapter 3 (incidence design) review and first CTJ revision: session plan

Written 2026-11-15 for the session started from
`docs/prompts/2026-11-14-chapter3-review-session-prompt.md`. Author: Andrew Walther.
Modeled on `multiomicsGEP/docs/plans/2026-11-14-chapter4-review-session-plan.md`. Section 4
(content map) and Section 5 (feedback triage) are proposals until Andrew confirms them.
Section 6 lists the decisions he still has to make.

## 0. Where edits go

- **Chapter 3.** The only editable source is `paper/thesis_chapter/Thesis_Chapter.qmd` (and the
  figure scripts for exhibits built by code). Committing it runs the post-commit hook, which
  re-renders `bios-dissertation/prelim/project-proposals/project2-incidence/draft/` and commits
  it there. Never edit or commit the draft folder by hand (the bios-dissertation pre-commit
  guard refuses it). Andrew's later edits in the draft `.tex` are carried back with the
  `prelim-edit-sync` skill.
- **CTJ.** Edit `paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`;
  rebuild `CTJ_Manuscript.pdf` and `CTJ_Reading.pdf` after each edit (`latexmk` in that folder,
  TinyTeX on PATH). The supplement is `paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex`.
  The hook copies the committed manuscript to `bios-dissertation/manuscripts/project2-incidence/`.
- **Markup copies.** When Andrew approves a CTJ version as final for this round, copy
  `CTJ_Reading.pdf` to `revision_markup/CTJ_Reading_Markup_V6.pdf` (V1 to V5 exist), commit and
  push it. Never edit an existing copy.
- Render with RStudio's Quarto 1.10. After every render, compare page counts with the baseline
  below. After any edit to Chapter 3 or CTJ, say whether the other needs it.

## 1. Baseline (Phase 1, 2026-11-15)

Deadline: Andrew confirmed a full working day is available (2026-11-15); the current chapter is
the fallback if needed.

| Item | Commit | Pages |
|---|---|---|
| SpatialCRT `main` (clean, matches origin) | `74f8d48` | |
| Chapter 3 source, last change | `9e64dca` (2026-11-13) | |
| Chapter 3 standalone draft PDF (synced, bios-dissertation `a7e5127`) | | 95 |
| Chapter 3 in the full prelim (`prelim.pdf`, bios `204d342`, 245 pp) | | body pp. 75–126 (52 pp); Appendix B pp. 176–212 (37 pp) |
| Chapter 3 exhibits | | 27 figures, 28 tables |
| CTJ_Manuscript.pdf (source `a685b4a`, PDF current) | | 22 |
| CTJ_Reading.pdf | | 9 |
| Supplementary_Information.pdf (`9f3f392`) | | 64 |

Both repositories are clean apart from two untracked files in SpatialCRT and two modified
`.gitignore` files in bios-dissertation that belong to other work. The Chapter 3 draft PDF is
dated after the last source commit. The CTJ PDFs were built from the current `.tex`.

The prompt's figures (51 body and 34 appendix pages, as of 2026-11-05) are out of date: the
chapter has grown to 52 and 37.

## 2. Goals

1. Make the targeted CTJ revisions from Dr. Lin's meeting and written annotations, then carry
   them into Chapter 3. CTJ gets its first revision round on the way. The CTJ 3,500-word limit
   is not a goal now.
2. Chapter 3 ready for the prelim: a finished study, with a body that carries the argument and
   an appendix that holds detail a reader can consult. Curate, don't dump; hold-outs go to
   `paper/thesis_chapter/notes/trimmed-material.md`, never deleted.
3. Notation and math correct for the code (Phase 3b), consistent across Chapter 3, CTJ and the
   supplement (Phase 3), with the β/τ/γ reparameterization applied if approved (Decision A).
4. Every FCL annotation and meeting item has a recorded disposition (Section 5).

## 3. Settled decisions and rules (from the prompt and the meeting notes)

- D-1 to D-5 in `docs/meetings/2026-11-14-ctj-feedback-meeting-notes.md` are settled:
  β on **Z** (intervention effect), τ on **X**, spillover γ S(**Z**) with S(**Z**) = **WZ**;
  the incidence coefficient is pre-specified, its value stated in the simulation settings;
  "Application" becomes "Case Study"; t and normal intervals cited; estimating the incidence
  coefficient after a trial is future work. Q-01 is settled (binary spillover only where it
  already exists, the case-study forms comparison).
- Edit order: CTJ and the supplement first (sign-off list, Phase 4a), then the Chapter 3
  carry-over list (Phase 4c) with the content-map restructuring.
- Standing author rules in `AGENTS.md`: one design-name set; t + df intervals with Wald beside
  them; SRS in every design test; Checkerboard in the main comparison; "Chapter 2" in the chapter
  only, never in CTJ; no repository links; incidence guidance framed as matching Saturation
  Quadrants' accuracy while targeting high-burden regions; three-decimal rule for MSE values.
- Writing: human-writer skill; no mannered prose; avoid semicolons and dashes as punctuation;
  define abbreviations at first use; Title Case headings; American English; "fitted" only
  before a noun ("is fit", "we fit", "a line fit to").
- Tense: completed work past, future work future.
- Appendix: one numbered list (B.1, B.2, ... in the prelim), methods sections first, then
  results, no "Part" headings.

## 4. Content map (confirmed 2026-11-15 with Andrew's edits)

Word counts are from the current source. "Keep" means keep in the body; "shorten" means keep
with less text; "move" means move to the appendix; "hold out" means move to
`trimmed-material.md`.

### 4.1 Introduction (1,350 words)

| Subsection | Proposal |
|---|---|
| Background & Motivation (832) | Keep. Add citations for "how much contamination occurs can depend on the spatial arrangement" (I-01). Keep the existing Chapter 2 / Checkerboard definition (I-02 done 2026-11-11). |
| Aims & Contributions (516) | Keep. "Application" → "case study" (D-3); third aim and second contribution reworded as planning a SUD-prevention education pilot trial (I-06, FCL-012, FCL-015). Outline paragraph updated to the confirmed headings. |

### 4.2 Methods (2,980 words): five subsections, renamed

| Proposed heading | Content | Now | Proposal |
|---|---|---|---|
| Spatial Weights and Outcome Model | **W** defined with its elements (M-03); one-line model with ε, σ and ρ described in prose (M-02); β **Z**, τ **X**, γ S(**Z**) (D-1); "cluster-level incidence" (M-04); the incidence term as an outcome shift, not a baseline (M-05); the coefficient's value moved out (M-06); queen/rook choice moved out (M-07); spillover regimes; one sentence that Chapter 2 used a binary spillover indicator and where the binary form is compared (Q-01) | Spatial Structure, Outcome Model and Estimand (580) | Keep, revise |
| Incidence Surfaces | The three generating models written out in full, with where ρ_X enters (S-02); why incidence is on the (0, 1) scale (S-04); why five configurations (S-03); grid only, the SUDDEN data enter only the case study (S-01); "surface" and "configuration" defined (FCL-035, FCL-042) | Heterogeneous Baseline Incidence (514) | Keep, revise; heading drops "Baseline" |
| Treatment Assignment Designs | `tab:designs`, `fig:design-samples`. First paragraph's eight-design comparison result shortened to one sentence pointing to the appendix (it is a result) | (397) | Keep, shorten paragraph 1 |
| Estimation and Performance Criteria | Oracle and non-oracle fits; MSE defined in words (M-12); t + df interval and Wald with citations, general quantile notation (M-10, D-4); identifiable / not identifiable (M-08); m(a) defined in words and formula (M-11); drop the "only one possible assignment" sentence (M-09) | (628) | Keep, revise |
| Simulation Design | Parameter values including the incidence coefficient's value and why it is pre-specified (M-06, D-2, Q-02); queen primary, rook sensitivity (M-07); scenario counts, the 250-fit structure, MC standard errors, design tests, allocation-risk sub-study. Seeding paragraph shortened (detail already in appendix) | Simulation Design, Resampling, and Metrics (862) | Keep, shorten |

### 4.3 Results (5,000 words, 9 subsections → 6)

| Proposed heading | Exhibits | Now | Proposal |
|---|---|---|---|
| Design Performance Under Queen Contiguity | `tab:performance`, `fig:mse` | (436) | Keep. Move the Nemenyi / Wilcoxon detail of paragraph 3 to the appendix, keep one sentence. The rank-vs-mean reversal paragraph shortened |
| The Spillover Regime Changes the Ranking | `tab:grid-srs` (absorbs `tab:regime`) | Spillover Regime (451) | Keep. Merge `tab:regime` into `tab:grid-srs`, which already has regime-specific MSE for all designs and SRS; one table instead of two |
| Sampling Designs Compared to Simple Random Sampling | `fig:vs-srs`, `fig:grid-srs` | Comparison With Simple Random Sampling (291) | Keep |
| Robustness to Incidence Pattern, Spatial Dependence and Spillover Magnitude | `tab:robustness` | Incidence Configuration, Spatial Dependence and Spillover Magnitude (230) + Variation Across Parameter Settings (153) | Merge. `tab:by-config` moves to the appendix (`fig:grid-srs` already shows the configurations). States that the configurations did not change the best design (S-03) |
| Stability Across Treatment Effect Sizes | `fig:tau-mse` | (537) | Shorten to about 300 words; block-level counts move to the appendix |
| Allocation-Specific Risk on the Grid | (points to `tab:downside`, `fig:alloc-dist`) | Allocation-Specific Risk (254) | Keep |
| (moved to Appendix B) | | Rook Contiguity (371) + Non-Oracle Sensitivity Analysis (215) | Andrew: not a primary result. Both move to the appendix (rook with B.9, non-oracle with B.10); the Results keep one sentence each pointing there. Results now has 6 subsections |

### 4.4 Case Study (5,300 words, 8 subsections → 7)

Heading (D-3): **Case Study: Simulation Based on Sudden Unexpected Death in North Carolina.**

| Proposed heading | Exhibits | Now | Proposal |
|---|---|---|---|
| Purpose, Observed Inputs and Simulated Outcomes | `fig:service-areas`, `fig:nc-observed` | Purpose and Observed Inputs (372) | Keep; says plainly what is observed and what is simulated (A-02) |
| Design Adaptations to the North Carolina Map | `tab:nc-budgets` | Geographic Adaptations and Treatment Budgets (268) | Keep; how the blocks and regions were built, from the code (A-03) |
| Outcome Model and Evaluation | | Outcome Model, Evaluation and Sensitivity Analyses (807) | Shorten to about 450 words; the CAR, Bayesian and spillover-form model detail moves to the appendix with one summary sentence left |
| Partial Participation and Sample Size | | (549) | Shorten; the ICC translation paragraph moves to the appendix with one sentence left |
| Accuracy, Coverage and Reach | `fig:nc-mse` | (554) | Keep. `tab:nc-year` moves to the appendix (beside `tab:nc-annual`; the figure carries the comparison) |
| Allocation-Specific Risk in North Carolina | `tab:downside`, `fig:alloc-dist` | (469) + Sensitivity Analyses (357) | Keep; the sensitivity paragraph becomes a short closing paragraph here or its own subsection (Andrew's choice) |
| Power Analysis for Partial Participation | `fig:colleges-needed` | Colleges Needed Under Partial Participation (602) | Keep. `fig:power-n` moves to the appendix |

### 4.5 Discussion (1,850 words)

| Subsection | Proposal |
|---|---|
| Summary of Findings (602) | Shorten to about 400 words: it restates many Results numbers |
| Design Guidance by Spillover Regime (380) | Keep. Add where Incidence-Guided Saturation is weakest (FCL-059): both arms, larger τ, few participants, volunteers favoring high-incidence colleges |
| Connection to Prior Spatial CRT Work (88) | Keep as its own subsection, retitled; gains the Chapter 2 symbol mapping if not placed in Methods |
| Limitations (681) | Keep; trim items that repeat Results |
| Future Considerations (101) | Expand to about 250 words (Decision E): estimating the incidence coefficient after a trial (D-5); setting the coefficient for other outcomes (Q-02 alternatives); binary spillover on the grid; noisy planning incidence; the count outcome |

### 4.6 Appendix B (one numbered list, methods first)

| Proposed section | From |
|---|---|
| B.1 Simulation Mechanics: Parameter Grid, Seeds and Metric Formulas | A2 |
| B.2 Identification Under Checkerboard and Isolation Buffer | A7, first part |
| B.3 Estimation Error Depends on the Treatment Effect | A7, second part |
| B.4 Design Rules and Regions on the North Carolina Map | A8 method text, `fig:nc-regions`, `fig:nc-examples`, `tab:nc-order` |
| B.5 Case Study Methods for Participation, Power and Model Sensitivity | A9 method text plus the model detail moved from the body (CAR, Bayesian, spillover forms, ICC translation) |
| B.6 Comparison of All Candidate Sampling Designs | A1 |
| B.7 Pairwise Tests and Critical-Difference Diagrams | A3 plus the test detail moved from the body |
| B.8 Results by Incidence Configuration, Spatial Dependence and Spillover Magnitude | A4 plus `tab:by-config` |
| B.9 Coverage, Bias and Variance, and Rook Contiguity Figures | A5 plus the rook paragraph moved from Results |
| B.10 Sensitivity Analysis Without the Spillover Term (Andrew, 2026-11-15) | A6 plus the non-oracle paragraph moved from Results |
| B.11 Additional Case Study Results | A8 results (`fig:nc-setting-ratios`, `fig:nc-risk`, `tab:nc-annual`, `tab:nc-summary-sens`) plus `tab:nc-year`, `fig:power-n` |
| B.12 Case Study Participation, Power and Sensitivity Results | A9 tables (`tab:nc-power` ... `tab:nc-forms`, `fig:spread-gain`), including the binary/count/proportional spillover comparison (Q-01) |

Page budget, as a guide: body about 45 pages, Appendix B about 38. Andrew (2026-11-15): 95 pages
is too long; the Phase 4c list also names appendix content that can be held out (moved to
`trimmed-material.md`, never deleted) to shorten the chapter and appendix together.

## 5. Feedback triage (proposed)

Status: **to do**, **addressed** (already handled; listed for confirmation), **decline**
(with reason), **after** (after the prelim). "Where": C = CTJ, S = supplement, 3 = Chapter 3.
Chapter 3 locations are source lines in `Thesis_Chapter.qmd` as of `9e64dca`; CTJ lines are in
`CTJ_Manuscript.tex`. Rows from the consistency and math checks are added under 5.9 when those
files are in.

### 5.1 Introduction and abstract

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-01 | FCL-013 | "Spatial interference challenges the assumption" | C 141 | Reword ("Spatial interference violates the assumption ...") | to do |
| T-02 | FCL-008, I-01 | Cite "Two allocations with the same number of treated clusters ..." and add references to the first paragraph | C 141; 3 l. 49 (last two sentences) | Two or more verified citations (candidates: the BMC paper, Baird et al. 2018, Leung 2022, Halloran and Hudgens 2016); checked on Crossref before use | to do |
| T-03 | FCL-001, FCL-009 | Citation placement; restructure the design-stage methods sentence | C 143 | Adopt Dr. Lin's wording with the citation at the end | to do |
| T-04 | FCL-002, -003, -004, -010, I-03 | "In a strict checkerboard with edge neighbors ... its own treatment" is unclear | C 143 | Rewrite in plain terms: under rook contiguity every neighbor of a treated cluster is a control, so a cluster's share of treated neighbors is set by its own arm and the two effects cannot be separated; cite the BMC paper. Chapter 3 has no such sentence | to do (C only) |
| T-05 | FCL-011, I-02, TODO | Checkerboard = the earlier study's block stratified sampling; define once, then "Checkerboard" only | C 145; S ~158, ~197; 3 l. 55, 61 | State once in CTJ and SI. Chapter 3 done 2026-11-11 | to do (C, S); addressed (3) |
| T-06 | FCL-005, I-04 | "That work weighed ..." unclear; "that work" | C 145 | Name the prior study ("In that study, we judged each design by ...") keeping CTJ self-contained | to do (C); addressed (3, l. 55 already clear) |
| T-07 | FCL-006, FCL-007, I-05 | Move the Checkerboard definition earlier, drop "meaning complete randomization ...", reorder sentences, make clear the comparison is this paper's and covers more than two designs | C 145 | Rewrite the paragraph: prior work and its two designs first, then "In this paper, we compare six designs and SRS ..." as the closing sentence | to do (C); check 3 l. 61 |
| T-08 | FCL-012, I-06 | "Our application supports planning" → "In this paper, we plan ..." | C 147; 3 l. 449, 59 | Reword with "case study" | to do |
| T-09 | FCL-014 | Abstract: "error under an unfavorable realized assignment (allocation risk)" unclear | C abstract | Say it in words: how large the error can be for the particular assignment a trial happens to draw | to do (C) |
| T-10 | FCL-015 | "pilot" → "pilot trial" | C abstract; 3 l. 59, 445 | "education pilot trial" at first use | to do |

### 5.2 Model and notation

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-11 | M-01, D-1, FCL-030 | Reparameterize: β **Z**, τ **X**, γ S(**Z**); update derived statements (Checkerboard estimates β − γ) | C, S, 3; figure and table labels; code comments | One item per passage in the sign-off lists, driven by the symbol list in `consistency-check.md` | to do (Decision A on timing) |
| T-12 | M-02, FCL-016, FCL-018, -019, -020 | One-line equation; ε in prose; "noise; τ is in units ..." without parentheses | C 152–163 (the reading layout splits the equation); 3 l. 72–76 | One-line equation in both layouts; prose: "ε is independent standard normal noise, the intervention effect is in standard-deviation units, and ρ is ..." | to do |
| T-13 | M-03 | Define **W**; avoid "the ith element of **WZ**" | C 151; 3 l. 69 | Define w_ij (1/number of neighbors if j is a neighbor of i, else 0), then (**WZ**)_i in words | to do |
| T-14 | M-04 | "cluster incidence" → "cluster-level incidence" | C 163; 3 l. 76 | As stated | to do |
| T-15 | M-05, FCL-021 | "baseline outcome" | C 163; 3 l. 76, 644, 702 | "shifts the outcome" / outcome shift; "baseline" kept only for historical incidence where it means incidence before the trial (to confirm) | to do |
| T-16 | M-06, D-2, FCL-022, FCL-047, FCL-053, Q-02 | Value of the incidence coefficient: move to the simulation settings, explain why it is pre-specified, why 1 on the grid and 0 / 1 in the case study | C 163, 192; 3 l. 76, 515, 522 | Text from Decision D; scale meaning from math check question (a) | to do |
| T-17 | M-07, FCL-017 | Queen/rook sentence "came from nowhere" | C 151; 3 l. 69 | Move to Simulation Design / Simulation study | to do |
| T-18 | M-13, FCL-023 | Split "We estimated τ, giving τ̂, by ML, fitting ..." | C 165 | Two or three sentences | to do (C); addressed (3, l. 138) |
| T-19 | M-10, D-4, FCL-024 | Cite the t interval, the df-corrected SE √(N/(N−p)) and the Wald interval; general quantile notation | C 165; 3 l. 140, 520; S | Write θ̂ ± q·SE with q = t_{0.975, N−p} or z_{0.975}; verified citations for the small-sample variance correction (to be found and checked); the "covered more reliably" claim points to the interval comparison table | to do |
| T-20 | M-12, FCL-025 | Define MSE at first use (Methods) | C 165; 3 l. 59 (defined in Aims) | Spell out at the Methods definition in CTJ; Chapter 3 already defines it at first use | to do (C); addressed (3) |
| T-21 | FCL-026, FCL-027 | "confidence intervals"; "coverage was measured" | C 165; 3 l. 156 | As suggested | to do |
| T-22 | M-08, FCL-028 | "identified" → "identifiable"; rook: not identifiable; queen: identifiable but imprecise; drop "(cannot be estimated)" | C 167, 204; 3 l. 142, 392, 437, 546, 692; S | Every use | to do |
| T-23 | M-09, FCL-029 | Delete "A design with only one possible assignment has no spread ..." | C 169, 229; 3 l. 639 | Delete in C 169; shorten the parallel sentences | to do |
| T-24 | M-11, FCL-014 | Define m(a) clearly at first use | C 169; 3 l. 144 | One sentence in words before the formula | to do |

### 5.3 Simulation design

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-25 | FCL-031 | Missing comma after "regimes" | C 180 | As suggested | to do |
| T-26 | S-02, FCL-032, FCL-033, FCL-034, FCL-036 | Write the model for **X**, including ρ_X; "One normal and one Poisson?"; "How?" | C 180; 3 l. 87 (formulas present); S | CTJ gets the formulas (short form); Chapter 3 clarifies that the spatial mechanism is a normal SAR field mapped by Φ and the Poisson mechanism is counts from a log-normal SAR rate, each at ρ_X = 0.2 and 0.5 | to do (C, S); check (3) |
| T-27 | S-03 | Justify or reduce the five configurations | C 180; 3 l. 87, 229 | Keep all five (already run): state what each tests and that none changed the ranking; detail table to the appendix (content map) | to do |
| T-28 | S-04 | Why incidence is on (0, 1) | C 180; 3 l. 87 | State the reason (a common scale, so one coefficient value means the same in every configuration; rank-scaling for counts) | to do |
| T-29 | S-01 | SUDDEN data used only in the case study | C 180; 3 l. 69, 87 | One sentence | to do |
| T-30 | FCL-035, FCL-042, W-04 | "surface k", "configuration" undefined | C 180, 182; 3 l. 87–89 | Define "incidence surface" (one draw of **X** over the grid) and "incidence configuration" (generating mechanism and ρ_X) at first use | to do (C); check (3) |
| T-31 | FCL-041, W-03 | Remove "so" | C 182; 3 l. 156 | "... and therefore carries ..." | to do |
| T-32 | FCL-043 | Is the grid allocation-risk paragraph a result? | C 184 | Decline the move: it describes the design of a separate study. Reword to say so ("We also ran a separate allocation-risk study ...") | decline move; to do (wording) |

### 5.4 Designs and Table 1

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-33 | FCL-037, FCL-038 | Table 1's "Incidence surface" row is not a design; "cluster" → "region"? | C Exhibits Table 1 | Relabel the row ("Example incidence surface") so it reads as context, not a design. Decline "region": the cluster is the unit; "region" means the four saturation groups | to do; decline (FCL-038) |
| T-34 | FCL-039 | NC Checkerboard rule unclear | C Exhibits; 3 l. 483 and Appendix A8 | Plain wording: visits clusters from most to fewest neighbors and gives each the arm held by fewer of its already-assigned neighbors | to do |
| T-35 | FCL-040, FCL-054, A-03 | How the blocks and regions were built (k-means) | C 194 and Exhibits; 3 l. 483 and A8; S | Code (math check (e)): the 14 blocks come from `kmeans` on four standardized features (projected coordinates, that year's rank-scaled incidence, log person-years), 25 starts, rerun each year (`application_designs.R:134-139`, `real_sud_setup.R:109-112`). It is **not** k-means on pairwise distances, as the meeting notes say. The four saturation regions are a separate k-means on coordinates only: 250 candidate partitions, scored by the criterion C. CTJ 194 already describes the blocks correctly; Chapter 3 l. 483 says only "nearby clusters". Wording: explain both in plain terms, in the body of both documents | to do |
| T-36 | FCL-044 | Isolation Buffer is not incidence-guided yet performs well: fair? | C Table 1, Results | Yes, under control-only spillover (no treated cluster has a treated neighbor), but not under both arms, and it treats fewer clusters. Already stated in Chapter 3 l. 698 and CTJ 239; reply only | addressed |

### 5.5 Case study

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-37 | D-3, A-01, FCL-048 | Rename "Application" → "Case Study"; should it follow the simulation results? | C 186, 219 and all "application" uses; 3 §Application and all uses; S | Rename everywhere. CTJ keeps its design in Methods and results in Results (meeting decision). Chapter 3 placement is Decision F | to do |
| T-38 | A-02, FCL-045 | Observed vs simulated; drop "rather than simulated" | C 187; 3 l. 453 | "Observed incidence and geography were fixed inputs. Outcomes were simulated." | to do |
| T-39 | FCL-046 | "Equation 1" → "Equation (1)" | C 192 and every equation reference | Check SAGE style, then apply | to do (C) |
| T-40 | FCL-060 | Results heading "Incidence-Guided Saturation is strong" → "the best" / "preferred"? | C 219 | Neither word fits, since Saturation Quadrants did as well: "The saturation designs had the lowest error under control-only spillover" | to do (claim group) |

### 5.6 Results and exhibits

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-41 | FCL-049 | "Incidence configuration changed each design's MSE by less than 10% ..." mentioned too early | C 204 | Move after the regime results | to do |
| T-42 | FCL-050, FCL-051 | Pooled 0.091 is not in Table 2; "separating the regimes was essential" | C 204; 3 l. 313 | Three-decimal part addressed. Drop the pooled figure from the CTJ results paragraph and lead with the regimes ("The pooled average hides opposite results under the two regimes"). Chapter 3 reports the pooled value in `tab:performance`, so it stays there | partly addressed; to do |
| T-43 | FCL-052 | "(SI)" is unclear | C 211, 213, 217, 222, 229 | Cite the SI section and table each time | to do |
| T-44 | FCL-055 | "outcome-dependence values" → "spatial dependence coefficient" | C Exhibits Table 2, Figure 2 captions; 3 captions | As suggested | to do |
| T-45 | FCL-056 | "the true spillover term" → "true γ?" | C Exhibits Table 2 caption; 3 l. 138 | After D-1: the fitted model includes the true exposure S(**Z**) and estimates γ. The code fits the column γ**WZ** (coefficient 1). Fitting **WZ** with coefficient γ gives the same τ̂ and SE, because rescaling a covariate leaves the other ML estimates unchanged, so the D-1 description is exact | to do (with T-11) |
| T-46 | FCL-057 | "N/A" → "Reference" for SRS | C Exhibits Table 2; 3 `tab:grid-srs` | As suggested | to do |
| T-47 | FCL-058 | Under both arms, no design beats SRS: fair to say? | C 204, 239 | Already stated in both; the answer recorded in the annotation agrees | addressed |
| T-48 | FCL-059 | How sensitive is Incidence-Guided Saturation; where can it fail? | C Discussion; 3 l. 698 | One short paragraph from existing results: 8% above SRS under both arms, control-only MSE rising with τ (Isolation Buffer lower at τ = 3), advantage shrinking with fewer colleges and gone when high-incidence colleges volunteer more | to do (no new analysis) |

### 5.7 Writing throughout

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-49 | W-01 | Unneeded parentheses | all | Sweep in the sign-off lists; the Phase 5 writing review checks again | to do |
| T-50 | W-02 | Long sentences with many commas | all | Same | to do |
| T-51 | W-03 | "so" as a connective | all | Same | to do |
| T-52 | W-04 | Undefined terms ("tractable" in the meeting; none found in the current CTJ text) | all | Same | to do |

### 5.8 Open questions, future work and pending TODO items

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-53 | Q-01 (settled) | Binary vs proportional spillover | 3 Methods and Appendix A9 (`tab:nc-forms`); C (previous work, self-contained); S | One sentence in Methods: the earlier study used a binary indicator; the case study compares binary, count and proportional forms (appendix); no new analysis | to do |
| T-54 | Q-02 | What the incidence coefficient's value represents | 3, C | Decision D | decision |
| T-55 | D-5 | Estimating the incidence coefficient after a trial | 3 Future Considerations; C Continued work (optional) | One or two sentences | to do |
| T-56 | FCL-050 (rest), FCL-058 | see T-42, T-47 | | | |
| T-57 | CTJ word limit | 3,500 words | C | Not now | after |
| T-58 | TODO "Full SI prose review" | | S | Only passages touched by this round now; full review after the prelim | after |

### 5.9 Consistency and math checks

Both checks were completed on 2026-11-15
(`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/`). The math check
found no wrong or undefined math: 107 rows, 71 verified, 15 cited, 15 uncited and 6
inconsistent. The consistency check found 11 mismatches. Every row below becomes a sign-off item.

| Row | IDs | Issue | Where | Proposed change | Status |
|---|---|---|---|---|---|
| T-59 | CC-01 | CTJ case-study coverage 94.2–94.9% (t) and 92.7–93.6% (Wald) omits Saturation Quadrants; the correct values are 94.1–94.9% and 92.6–93.6% | C 53, 222 | Use the chapter's values (`a5_stage1_coverage_by_design_regime.csv`) | to do (number) |
| T-60 | CC-02 | CTJ "median 63%" noise share is cited to the SI, which does not state it | C 217; S ~1890; 3 l. 431 | Add the grid medians (63% control only, 92% both arms; `conditional_summary.csv`) to the SI and Chapter 3 | to do |
| T-61 | CC-03 | CTJ says the SI has per-year cluster allocations; it has budgets, orderings and one example | C 276 | Reword | to do |
| T-62 | CC-04, MC-55 | n and N both used for the number of colleges | 3 l. 520 vs 534; S 818 vs 855 | N throughout | to do |
| T-63 | CC-05, MC-54, SI-31 | Reused symbols: **S** (spillover and CAR matrix), λ, m, d, n, P, and ρ written for the CAR parameter | 3, S | Rename the secondary uses (CAR matrix, professionals per college, generic metric) | to do |
| T-64 | CC-06 | σ² **I** (chapter, SI) vs **I** (CTJ) | 3 l. 73; S 59 | Settles with T-12: equation without ε's distribution, σ = 1 stated in prose | to do |
| T-65 | CC-07 | Bias and SD at two decimals | C 215; 3 l. 437, 441, 646, 1600; S | Three decimals per the MSE precision rule | to do (number format) |
| T-66 | CC-08 | Ratios to SRS: two decimals in tables, three in 3 l. 313 and S 1832 | C Table 2; `tab:downside`; 3; S | Decision H | decision |
| T-67 | CC-09 | Grid power one-sided (0.025), case-study power two-sided (0.05); the case-study sidedness is never stated | 3 l. 520; S 818; C 200 | State both sidednesses; CTJ 200 names the t interval | to do |
| T-68 | CC-10 | Coverage and power captions do not name the interval | 3 `tab:performance`, `tab:by-config`, `tab:regime`; S; C 215 | Add "t-based" | to do |
| T-69 | CC-11 | "BSS" reused at 3 l. 686 | 3 l. 686 | "SRS and Checkerboard" | to do |
| T-70 | MC-32, SI-22 | `rate` in the rank formula prints as italic | 3 l. 522; S 823 | `\text{rate}` | to do |
| T-71 | MC-56 | "**WZ** = 1/2" sets a vector equal to a scalar | 3 l. 1233; S 342 | (**WZ**)_i = 1/2 | to do |
| T-72 | math check, uncited | t + df interval, γ**WZ** spillover form, MSE decomposition, MC SE formulas (morris_using_2019 is in the bib), SAR equation in the SI | 3, C, S | Verified citations, staged in `round5/ch3/new-citations.bib`; merges with T-19 | to do |
| T-73 | math check, question | "The application uses the same interval" follows N = 100, p = 5 | 3 l. 140 | Say the case study uses N = 58, p = 4 (df 54) | to do |
| T-74 | math check, question | 4 aliased rook Isolation Buffer fits keep p = 5 in the df | 3 Appendix A7 | Footnote or leave (a difference of one degree of freedom) | to do (minor; propose leave plus one clause) |

Facts the checks settled, for the sign-off wording:
- **Scale of X (Q-02, S-04).** X lies in (0, 1) in every configuration. With σ = 1, an
  incidence coefficient of 1 means the outcome rises by one noise standard deviation from the
  lowest to the highest incidence, before spatial spread. That is the same size as the treatment
  effect at its primary value. The incidence term carries about 8% of the variance of the
  incidence term plus noise. The spatial configuration is only roughly uniform (SD 0.29–0.30).
- **Model for X (S-02).** The spatial configuration is **X** = Φ{(**I** − ρ_X **W**_Q)⁻¹**u**}.
  The Poisson configuration takes log λ = log(35/100,000) + (**I** − ρ_X **W**_Q)⁻¹**u**, draws
  counts C_i ~ Poisson(100,000 λ_i), then rank-scales the counts. The independent configuration
  has no ρ_X. Queen **W** is always used.
- **Checkerboard (M-01, M-08).** Under rook contiguity **WZ** = **1** − **Z** exactly, so the
  spillover column is dropped and the coefficient estimates τ − γ (β − γ after D-1, a plain
  symbol swap). Under queen contiguity the effect is identifiable in both regimes: imprecise
  under control-only spillover (the spillover column has R² 0.981 on [1, **Z**]) and precise
  under both arms (R² 0.347).
- **Partial participation.** The spillover term is κτ e(**Z**), with e divided by each
  college's full neighbor count, so it is not **W**_P **Z**. The text states this correctly.

## 6. Decisions for Andrew (Phase 2d)

- **A. Reparameterization before the prelim?** It changes symbols in prose, equations, table
  headers and figure axis labels (figures rebuilt by `code/28_thesis_chapter_figures.R` and the
  CTJ figure scripts). Recommendation: yes, for Chapter 3, CTJ and the supplement together, so
  the three never disagree. Size, from `consistency-check.md`: τ appears on about 200 lines of
  the chapter, 17 of CTJ, 3 of the CTJ exhibits and about 250 of the SI. β appears on about 25,
  5 and 23 lines. 19 grid scripts contain τ tokens, and the application scripts were not
  scanned. Only figures that print τ or β in a label need rebuilding (at least `fig:tau-mse`,
  `fig:tau-coverage`, `fig:power-n`, `fig:colleges-needed`; to be listed exactly before
  editing). The swap is mechanical but large. It must be a true swap (τ→β and β→τ in one pass),
  leaving the earlier study's β alone (SI 336–338). Recommendation still yes, as its own
  commit, checked by grep and by the verifier.
- **B. Chapter 2 notation.** Chapter 2 (from the published paper) writes
  y = α + β x + ψ z + ρ **W**y + ε, with x the treatment indicator, z the spillover exposure and
  ψ the spillover coefficient. After D-1, β means the intervention effect in both chapters, but
  Chapter 2's x, z and ψ still differ from Chapter 3's **Z**, S(**Z**) and γ, and Chapter 3
  uses **X** for incidence. Recommendation: leave Chapter 2 as published and add one sentence in
  Chapter 3 mapping the symbols. The alternative is changing Chapter 2 to **Z**, S(**Z**) and γ.
- **C. What waits until after the prelim.** Proposed: the CTJ word limit, the full SI prose
  review beyond touched passages, any new analysis (none is proposed).
- **D. Q-02, the incidence coefficient.** Proposed framing: the coefficient is set by the
  investigator according to whether incidence is believed to shift the outcome. The grid uses 1
  (incidence shifts the outcome by one noise standard deviation across the full incidence range,
  the harder case for incidence-based designs); the case study uses 0 for the education outcome
  and 1 as a hypothetical sensitivity analysis. Confirmed by the math check: with incidence in
  (0, 1) and unit noise variance, a value of 1 moves the outcome by one noise standard deviation
  from the lowest to the highest incidence, the same size as the primary treatment effect.
- **H. Precision of ratios to SRS (CC-08).** Does the three-decimal rule cover ratios to SRS?
  Table 2 and `tab:downside` use two decimals, while Chapter 3 l. 313 and the SI use three.
  Recommendation: ratios at two decimals everywhere, a stated exception to the rule, because
  they are already relative and their intervals are wide. The alternative is three decimals
  in every table.
- **I. A-03 wording.** The meeting notes say k-means on pairwise distances. The code clusters
  standardized coordinates, incidence and log population. The text will follow the code.
- **E. Future Considerations scope.** Proposed: about 250 words covering D-5, alternative ways
  to set the coefficient, binary spillover on the grid, noisy planning incidence, other
  geometries and the count outcome.
- **F. Case study placement in Chapter 3.** Keep it as its own top-level section after Results
  (as Chapter 2 does with its Application), renamed; or split it, with its design in Methods and
  its results in Results (as CTJ does). Recommendation: keep it as its own section.
- **G. Content map.** Confirm or edit Section 4, including the merge of `tab:regime` into
  `tab:grid-srs` and the moves of `tab:by-config`, `tab:nc-year` and `fig:power-n`.

### Answers (Andrew, 2026-11-15)

- **A. Yes.** Reparameterize in the text so β is the intervention effect in Chapters 2 and 3.
  τ becomes the incidence coefficient. The spillover coefficient γ becomes ψ, as in Chapter 2
  (checked: γ means only the spillover magnitude in the chapter, CTJ, exhibits and SI, and ψ
  appears only in the two sentences quoting Chapter 2's notation, so the swap is mechanical).
  Code internals (`MLE_tau_sweep`, column names) stay as they are; only plot label strings
  change, and the affected figures are rebuilt from stored results. Anything heavier is deferred.
- **B.** Chapter 2 is not changed. Chapter 3 keeps **Z** (treatment) and **X** (incidence),
  since renaming them to Chapter 2's x and z would collide with **X**, and adds one sentence
  mapping the symbols (β and ψ now match; Chapter 2's x and z are Chapter 3's **Z** and
  S(**Z**)).
- **C. Yes.** CTJ word limit and the full SI prose review wait until after the proposal. The SI
  still gets the reparameterization and the items from the checks, so it never disagrees with
  CTJ.
- **D. Yes.** The coefficient is the investigator's choice; future work can estimate it. Explain
  the meaning of 1 clearly (the current implementation) and say it may need adjusting in practice.
- **E. Yes.** About 250 words.
- **F. Yes.** Its own section, framed as a case study using simulation with real incidence data.
- **G.** Confirmed with the renames above and the sensitivity subsection moved to the appendix.
- **H.** Leave ratios as they are now; no note. Chapter 3 l. 313 and SI 1832 (three decimals)
  are left alone unless feedback asks.
- **I. Yes.** Follow the code and explain how each design is adapted.

Answers on the check rows (T-59 to T-74):
- **T-72 citations.** t + df interval: cite (verified on Crossref) MacKinnon and White 1985,
  J Econometrics 29:305–325, doi:10.1016/0304-4076(85)90158-7, for the N/(N−p) variance
  correction, and Leyrat et al. 2018, Int J Epidemiol 47(1):321–331, doi:10.1093/ije/dyx169,
  for t intervals with reduced degrees of freedom in trials with few clusters. MC SE formulas:
  `morris_using_2019` (in the bib). SI SAR equation: `ord_estimation_1975` and
  `lesage_introduction_2009` (in the bib). γ**WZ** form: no citation (follows from the SAR
  setup). MSE decomposition: no citation (standard).
- **T-63.** **S** for the CAR matrix stays, with a clear definition where it is used. The other
  reused symbols (m, n) are renamed.
- **T-62.** N for the number of colleges (clusters analyzed) everywhere; n kept for samples
  within a cluster.
- **T-59.** Andrew asked whether the gap is rounding. It is not: the CTJ range is the
  seven-design range (Saturation Quadrants and Balanced Halves left out). Saturation Quadrants'
  control-only coverage is 0.94115 (t) and 0.92625 (Wald), Incidence-Guided Saturation's
  0.94205 and 0.92651. Use 94.1–94.9% and 92.6–93.6%.
- **T-60, T-65, T-67.** Yes: add the 63% / 92% medians to the SI and chapter; bias and SD to
  three decimals; state the sidedness of both power analyses.

## 7. Phases

1. Set-up: done (Section 1).
2. Plan, content map and triage: this file. Decisions received 2026-11-15.
3. Consistency check (subagent): `bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/consistency-check.md`.
3b. Math check (subagent): `.../round5/ch3/math-check.md`.
4a. CTJ sign-off list: `paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/proposed-edits.md`. Written 2026-11-15; decisions received 2026-11-15 (all accepted, six edited).

Model use (prompt, added 2026-11-15): from 4b on, approved edits are applied by Sonnet subagents item by item, and I check each diff before committing. Render, page-count, cross-reference and citation checks and the writing review also go to Sonnet. Restructuring decisions, new prose, and the committee review stay on Opus.
4b. Apply accepted CTJ and supplement items; rebuild; commit.
   Notation swap committed 2026-11-15 (`7e51487`). Figure labels: the SI's 6-design figures
   (ranked MSE, coverage, β sensitivity, `app_tau_coverage`) are THESIS_BUILT, drawn by
   `code/28_thesis_chapter_figures.R` and copied to the CTJ folder by `paper/tools/copy_figures.sh`.
   Their "True τ" labels change only when script 28 is edited and rerun, which also redraws the
   chapter's figures. Deferred to Phase 4d (chapter notation), then rerun 28 and the copy script
   once for both. Script 28 also has a "rho=…, gamma=…" label (l. 466); `real_sud_summary.R`
   has ASCII "rho/gamma" in captions.
4c. Chapter 3 carry-over list with the restructuring: `.../round5/ch3/carry-over-edits.md`. **Stop.**
4d. Apply section by section; **stop after each section**; commit after each approved group.
   **Status 2026-11-15:** 4b done (CTJ/SI: 7e51487, 43b6f99, 09803e1, 324cd80; trim plan in the
   revisions folder, deferred). 4c/4d done in SpatialCRT: globals fcf0c70, Intro+Methods 1eaf4da
   (+ bib fix 44eb69a, bios 1ba1609), Results a1558af, Case Study 227e3db, Discussion 67d9b39,
   Appendix 6b8da5a. Chapter 3 in the bios draft: 88 pp (was 95). Cut text is in
   `paper/thesis_chapter/notes/trimmed-material.md`. Still open from Phase 4: rebuild CTJ PDFs
   (bib title changes from 44eb69a). Next: Phase 5.
5. Writing review and committee review (parallel subagents, without this session's
   conclusions): `round5/ch3/proposed-edits.md`, `committee-proposals.md`, `new-citations.bib`.
   **Stop.**
6. Apply accepted items; re-render everything; check numbering, cross-references and page
   counts; V6 markup copy; update the triage table, ROADMAP.md, docs/TODO.md, DECISIONS.md here
   and bios-dissertation/ROADMAP.md.

## 8. Done when

- The content map is confirmed and followed.
- Every FCL and meeting item has a recorded disposition.
- The reparameterization (if approved) is consistent across all documents.
- The consistency check and the math check have no open rows.
- The CTJ and Chapter 3 sign-off lists have a decision on every item and accepted items are applied.
- Both review files have decisions and accepted items are applied.
- Chapter 3, CTJ and the supplement are in American English, with "fitted" only before a noun.
- Everything renders, with page counts recorded.
