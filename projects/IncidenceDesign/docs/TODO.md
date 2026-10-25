# IncidenceDesign to-dos — completed study and drafts

Updated 2026-10-13. The author superseded earlier deadline dates and requested ASAP
completion with CTJ remaining a required deliverable. Historical interview and
execution to-dos are preserved in Git history; settled choices are not reopened.

## Completed

- [x] SRS benchmark added and verified without altering original eight-design results.
- [x] Average MSE and allocation-specific downside interpreted together; grid pilot
  distinguished from the full study and finite tail uncertainty acknowledged.
- [x] Checkerboard retained in main comparisons; technical diagnostics in appendix.
- [x] Approved conceptual Chapter 2 link; CTJ standalone preceding-work citation.
  Project 1 is final/read-only; historical numerical discrepancy remains unresolved.
- [x] Author-approved education application: observed 2018–2021 inputs, 58 clusters,
  β=0 primary, constant τ, no synthetic incidence, all eight designs and SRS,
  queen primary/both spillover regimes/yearly reporting, geographic adaptations,
  unequal budgets and focused sensitivities.
- [x] Independent alignment review and concrete implementation approval.
- [x] Implementation, behavioral/ML checks, manifests, named cached weights,
  reproducible seeds and explicit warnings/alias/failure accounting.
- [x] Full production: 1,248 rows/936 computational sources/7,598,000 independent
  outcome fits; all complete and mean/coverage MC gates passed; no invalid fits.
- [x] Separate same-allocation tail extension: 144 rows/96 sources, 2,280,000
  additional outcomes; tail-membership uncertainty remains explicitly reported.
- [x] Offline HTML companion: development, model, maps, mean/tail performance,
  budgets, yearly filters and links to completed drafts; DOM behavior tested.
- [x] Full thesis chapter and appendix with number/exhibit provenance, resolved
  former application comments and fresh PDF review; approved Chapter 3 sync.
- [x] Full CTJ manuscript/supplement derived from the chapter; old results replaced,
  live journal requirements checked, confirmed declarations and explicit ethics
  placeholder, figure selection, legends and cover-letter draft.
- [x] Independent scientific source review, aggregate-to-table verification,
  fresh document builds and visual review.
- [x] Requested comprehensive Claude review prompt saved.
- [x] Author-requested paired CTJ layouts: shared-source two-column reading preview
  and double-spaced submission draft, rendered and numerically checked.

## Remaining author/submission work

- [x] Assemble full-length prelim in bios-dissertation, retaining the full
  Chapter 3 appendix; fix CTJ reading bibliography flow/last-page columns.
  See [assembly/reference review](findings/2026-10-02-prelim-assembly-ctj-reference-review.md).
- [x] Add author-requested List of Tables and List of Figures to the full prelim;
  verify all 63 current exhibit page references and preserve main-matter numbering.
- [x] Remove Chapter 4's fixed-date timeline in canonical multiomicsGEP source,
  sync and rebuild the prelim; frame continuing work as methodology and expected evidence.
- [ ] Author review of full prelim and initial front abstract. Chapter 4's later
  substantive revision will follow its October 2 plan and verified new results.
  Rebuild the assembly after later chapter syncs.

- [x] Implement approved title restoration and six-exhibit revision in both CTJ
  layouts; carry applicable clarity changes into the chapter and derived supplement.
  See [exhibit review](findings/2026-10-02-ctj-exhibit-revision-findings.md).
- [x] Chapter 3 length review (2026-10-12/13): concise version independently reviewed
  and adopted (prelim body 45.5 to 41 pp; Appendix B 27 pp; standalone chapter 68 pp;
  SI 45 pp with merged annual table); synced, full prelim rebuilt (192 pp), pushed.
  See [plan/status](plans/2026-10-12-chapter3-length-review-plan.md); cut text in
  `paper/thesis_chapter/notes/trimmed-material.md`.
- [x] Reference-list length check: strict APA 7 with DOIs and full author lists kept
  (a six-author et-al. variant saved no pages; reverted).
- [ ] Author full review of the revised (concise) chapter, both CTJ layouts and supplement.
- [x] Comprehensive review (four tracks: numbers, journal compliance/clarity,
  chapter consistency, code) completed 2026-10-12; approved fixes applied. See
  [findings](reviews/2026-10-12-ctj-comprehensive-review-findings.md).
- [x] Application study report for the advisor meeting:
  `application/report/IncidenceDesign_Application_Study_Report.{qmd,pdf}`.
- [ ] **Ethics/IRB statement** (CTJ requires one in Methods): protocol filled in
  2026-10-21 (UNC IRB 14-2036, "Sudden Unexpected Death in North Carolina", PI Ross
  Simpson, Department of Medicine) in CTJ Methods and the Chapter 3 application
  section. Still pending from Simpson (expected in the coming days): the determination
  type and consent wording (highlighted CTJ placeholder; the chapter states the
  protocol only), any required disclosures, the title-page statement, and permission
  to publish 58-cluster aggregates and maps.
- [ ] **Generative-AI disclosure** (Sage policy requires one for AI-assisted code,
  code checking and visualization): draft in CTJ Acknowledgements (2026-10-14);
  author to confirm tool list and wording.
- [ ] Coauthor ORCIDs if the submission system requests them (Walther's is on the
  title page). Lin's telephone (919) 966-9431 added from his UNC profile.
- [ ] Confirm contribution wording; obtain all coauthors' final approval and
  submission attestations; confirm the pilot description with Simpson.
- [ ] **Coverage below nominal (~93%; added 2026-10-14; plan A5):** Wald intervals
  τ̂ ± 1.96·SE are slightly too narrow in every design (likely ML SE ignoring
  uncertainty in ρ̂; untested). Quick fix tried 2026-10-14 (application only): a
  t critical value (df = 58 − 4 = 54, t ≈ 2.005; 53 for β = 1) recomputed from the
  stored production fits (stored Wald coverage reproduced exactly, 936/936 sources)
  raised primary coverage only from ~93.1% to ~93.7% (per design/regime
  92.6–93.6% → 93.2–94.1%; β = 1: 92.9% → 93.5%); rankings unchanged, ~1.3 points
  of the gap remain. Grid not checked. Closing the gap would need a heavier option:
  parametric bootstrap; randomization inference (not possible for the
  single-allocation designs). Ask Lin whether a stated limitation (optionally with
  one sentence on the t-value check) suffices.
- [ ] Submit CTJ with the author-completed
  [checklist](../paper/ctj_manuscript/submission_checklist.md)/cover letter.
- [x] Completion work is pushed: `main` matches `origin/main` (verified 2026-10-13,
  SpatialCRT `f67ab4b`, bios-dissertation `7d10596`). Future pushes still require explicit author authorization.

No additional analysis was required by the final independent source review.
Actual pilot endpoint/enrollment/budget/network specification and interval
calibration belong to planning the implemented trial; this simulation evaluation
is complete and does not claim observed intervention effectiveness.

Findings: [application](findings/2026-10-02-nc-sud-application-findings.md),
[manuscripts and teach-back](findings/2026-10-02-manuscript-completion-findings.md).
Commands and safeguards: project/application README and AGENTS.

## Advisor meeting 2026-10-14 follow-ups

Context, Lin's feedback, unanswered questions and fact-checks:
[meetings/2026-10-14-advisor-meeting.md](meetings/2026-10-14-advisor-meeting.md).
**Approved plan for the new analyses (2026-10-20):
[plans/2026-10-20-methods-extensions-plan.md](plans/2026-10-20-methods-extensions-plan.md)**
(decisions in DECISIONS.md, 2026-10-20). Analyses A0–A7 ran 2026-10-20/21 (commits
23f1d7a–40efbf8, plus the `a3v2` rerun); see `findings/2026-10-20-methods-extensions-findings.md`.
Items overlapping the submission list above (IRB, coauthor approval, pilot
description, ORCIDs) are not repeated here.

**Decide first (it shapes the rest)**
- [x] Outcome–incidence relationship, decided 2026-10-20: Y is the education
  endpoint, β = 0 primary with β = 1 beside it, constant τ, no count model (plan D1).
  Confirm with Lin/Simpson at the upcoming meeting. Original options: (Lin's "big question"; needs careful
  consideration, depends on how we frame incidence vs outcome): (a) keep β = 0 and
  state the endpoint-tracks-SUD assumption, citing grid β = 1 and the existing NC
  β = 1 runs; (b) promote the existing NC β = 1 runs (all 432 queen settings, no new
  simulation; reverses the 2026-10-02 rule, numbers change across documents);
  (c) model a future-incidence endpoint (new model/runs); (d) report both β = 0 and
  β = 1 side by side (author's leaning, 2026-10-14; no new simulation; decide which
  leads); (e) incidence-dependent treatment effect (Z·X; new runs). First settle what Y
  is: an education measure (β = 0) or future SUD incidence (β > 0); see the meeting
  record's "Follow-up reasoning". Record the choice in DECISIONS.md.
- [ ] Get Lin's sign-off on the open prep questions (design choices, regime-based
  framing, six designs + SRS, deterministic designs as comparators, code availability,
  target date/backup journal): list in the meeting record.

**Coauthors and data (critical path)**
- [ ] Meet with Habib and Simpson: review the application study, share the CTJ draft
  and results, ask for approval.
- [ ] From Habib/Simpson: IRB protocol holder and determination, permission to publish
  58-cluster aggregates/maps, required disclosures for derived data, and funding
  (does "no specific grant" hold?).
- [ ] Replace "personal communication" for Habib's corrected counts with a citation
  Lin finds acceptable (e.g. a data-preparation description or acknowledgement).
- [ ] **Simpson feedback on the application (awaited, 2026-10-21):** whether our
  implementation of the SUD application matches his pilot idea; revise the
  application framing in CTJ, SI and Chapter 3 accordingly.
- [ ] Frame the application as a statewide expansion of Simpson's small-scale (Lenoir
  County) proposal; fix the "pilot" wording flagged before the meeting.

**New analyses** (sample size/power and reduced participation must be addressed in the
thesis chapter and the prelim talk, main or backup slides; author, 2026-10-14)
- [x] Partial participation (plan A3; run 2026-10-20, rerun as `a3v2` 2026-10-21; findings Section 3): sample subsets of colleges (about half; Zoom said ~20),
  rebuild W on the possibly non-contiguous subset, and check whether the regime
  guidance (incidence-guided saturation under control-only spillover) still holds.
- [x] Sample size and power (plan A0, A3, A4; run 2026-10-20; findings Sections 2–3): answer "how many clusters are needed?". Power at τ = 1 for
  the current NC runs can be computed from stored fits (τ̂, SE) without new
  simulation; subset runs would need it too. Candidate for the main CTJ text, not
  only the SI.
- [x] Design comparisons with sample size (plan A2, A3; author, 2026-10-19; run 2026-10-20, findings Sections 2–3): show with formal tests
  which designs are better or worse, and the number of participating colleges needed to
  make that claim (power curves by design). Note that more simulation replicates alone
  make any difference "significant". Planned-extensions reasoning is in the meeting record.
  The application report now lists these analyses (Section 8, "Planned extensions").
- [x] Estimation-method sensitivity (plan A6; run 2026-10-20, long-chain Bayesian CAR 2026-10-21; findings Section 4): design guidance under CAR and/or Bayesian fits
  versus SAR (existing related evidence: grid non-oracle and omitted-spillover fits).
  Likely committee question; scope with Lin.
- [x] Binary and additive-count spillover sensitivities (plan A7; run 2026-10-20; findings Section 4): refit with Project 1-style on/off
  spillover (any treated neighbor). Expected without running it: Isolation Buffer and
  Checkerboard not identified, SRS/balanced worse, saturation relatively better. Could
  be a sensitivity analysis or a prepared committee answer.
- [x] Equal-weight pooling check (2026-10-14): Σd/ΣP versus `mean_of_yearly_rates`
  across 58 clusters: rank correlation 0.9999, 2 one-place rank swaps, no quartile or
  top-29 changes. No meaningful impact; pooled rates do not drive yearly allocation.

**Manuscript and exhibits (edit the chapter first, then re-derive CTJ/SI; run verifiers)**
- [ ] Group the designs into categories in Table 1 (e.g. deterministic, random,
  incidence-driven, saturation-driven; some designs span two). Consider one unified
  table + illustration exhibit.
- [ ] State explicitly that the endpoint is assumed correlated with the guiding
  incidence (wording depends on the outcome decision).
- [x] Update keywords (done 2026-10-19: dropped "simulation" and "mean squared error";
  added interference, design-stage treatment allocation, incidence-guided design).
- [ ] Optional: one sentence contrasting this project's proportional spillover
  (γ × share of treated neighbors) with Chapter 2's on/off indicator (keep brief).
- [ ] Keep CTJ word counts current after every edit (texcount, text + headers, math
  excluded; title page and checklist). 2026-10-14: body 3,362; 3,628 with
  declarations; abstract 393 (394 counting "10×10").
- [ ] Decide whether statistical tests of design differences belong in the CTJ/SI
  (the SI has Nemenyi/Holm-Wilcoxon for variant consolidation; script 12 has six-design
  Friedman/Nemenyi; none for the application).
- [ ] Re-run the saved read-only review prompt
  (`prompts/2026-10-12-claude-comprehensive-review-prompt.md`) after these edits.

- [x] Burden share treated by design (plan A1; run 2026-10-20, grid check 2026-10-21; findings Section 2).

**Prelim and beyond**
- [ ] Prelim slides (bios-dissertation/slides/): structure Project 1 / 2 / 3, brief
  literature review, depth on completed projects; prepare for general-science
  questions and the SAR vs CAR vs Bayesian question.
- [ ] Ask Lin which project to submit as an ENAR 2027 abstract.
- [ ] Final defense: focus on progress since the prelim (mainly Project 3).

### Methods-extension decisions for the author (2026-10-20 run)

Decided 2026-10-21 (`DECISIONS.md`). Details: `findings/2026-10-20-methods-extensions-findings.md`.
- [x] A5 coverage: report the t + df-corrected interval (Stage 1 main, Stage 2 in SI).
- [x] A1/D2 targeting: reworded (small gain on NC; larger where burden is regional);
  rate share primary, death share beside it.
- [x] V6: ≥ 2-college region rule, in the A3 rerun.
- [x] Power headline: size-adjusted, plus t + df-corrected power (A3 rerun); Wald in SI.
- [x] V8: no rerun; report dips with subset-bootstrap intervals.
- [x] Bayesian CAR: rerun ρ = 0.5 with 4× longer chains.
- [x] Overnight 2026-10-21: A3 rerun (`a3v2`), grid burden share, Bayesian CAR long
  chains, summaries, V9 (QUEUEDONE 10:41; one check-spec fix, logged).
- [x] Author: headline test = t + df; design name = "Saturation Regions"; pushed
  (2026-10-21).
- [x] Rename to "Saturation Regions" at the display layer (2026-10-25): chapter, SI and
  table labels, new figures, verifier labels. Older figures: see the item below.
- [x] Author reviewed Section 4 of the report (2026-10-25 feedback): SRS in every design
  test, model specs and citations, narrative results, no colons, ICC spelled out.
- [x] Step 7a (2026-10-25): figures first. `code/23_methods_extension_figures.R` (power
  vs N, colleges needed, spread vs gain), `code/15` vs-SRS forest plot, SRS in the grid
  tests (`code/12`, `code/14`).
- [x] Step 7b first pass (2026-10-25): Chapter 3, SI and CTJ (both layouts) drafted;
  verifiers pass; `CTJ_Reading_Markup_V3.pdf` made for markup.
- [ ] Author: review the step 7b draft (chapter Application sections and Appendix A9,
  SI S10 subsection, CTJ +~95 words; CTJ body now 3,483 of 3,500).
- [ ] Older manuscript figures still print "Saturation Quadrants" / "Incidence-Guided
  Saturation Quadrants" (design samples, tau sensitivity, coverage, ranked MSE, NC
  real_sud and compact figures). Their current script outputs also differ in content
  (SRS reference lines), so swapping them needs a review pass; the text notes the old
  names for now.
- [x] Step 8 (2026-10-25): draft talk outline in `presentation/OUTLINE.md` and
  `presentation/outline.html` (figures inline). The prelim Part II outline in
  bios-dissertation is unchanged; adapt it from the talk deck later.
- [ ] 20-minute project talk for clinicians and epidemiologists, Thursday 2026-10-29:
  full overview, discussion weighted toward the application; built in `presentation/`
  and later adapted into prelim Part II. Kickoff prompt:
  `prompts/2026-10-25-talk-kickoff-prompt.md`.
- [ ] Push the unpushed commits on main (ask the author).

## Known minor code issues (low priority; moved from AGENTS.md 2026-10-14)

- `02_incidence_generation.R`: `generate_incidence_poisson()` still defaults `pop_per_cluster = 1000` (and its roxygen says so). 05 always passes 100,000, so results aren't affected. It was left unchanged after the full run so the run manifest's code hash still matches.
- `DESIGN_FULL_NAMES` in `00_design_names.R` (used only by the two design-sample figures) labels Design 1 "Block Stratified Sampling (Checkerboard)" (user decision 2026-09-24: Chapter 2's name plus this study's). The manuscripts, 03 and 12–14 say "Checkerboard", and `DESIGN_SHORT_NAMES` still says "Block Stratified". This is not a bug to fix: the author uses "Block Stratified Sampling" and "Checkerboard" interchangeably, and the chapter says so once (Phase D).
- `load_latest_results()` comment (line ~87 of `06_visualizations.R`): clarify `_combined_` preference applies per estimation-mode, not globally
- `07_results_summary.Rmd` compare-table caption: should note that DIM only ran 6 designs (D7/D8 NAs are expected)
- `IncidenceSpatialCRT_Report.qmd` caption/text alignment: MC SEs table uses tau=1.0 slice (2,560 rows), not full 12,800 — prose now correctly clarifies this distinction
