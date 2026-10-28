# DECISIONS.md — IncidenceDesign

Author decisions with their reasons, newest first. `AGENTS.md` keeps only the
one-line rule each decision implies.

---

## 2026-10-27 — CTJ V5 revision: full-surface allocation risk, paired SRS column, Methods order

**Decision.** From the author's V5 review of the CTJ reading draft:
- The grid allocation-risk study is rerun on all ten surfaces per incidence
  configuration (same queen ρ ∈ {0, 0.5} × γ ∈ {0.5, 0.8} settings, nine designs, J = 100,
  R = 100 plus the R = 400 precision step) into `results/allocation_risk/surfaces10/`.
  Its numbers replace the two-surface pilot in the CTJ, SI and Chapter 3, and the word
  "pilot" leaves the manuscripts (except "education pilot"). This reverses the
  2026-10-02 choice not to extend the grid allocation-risk study.
- Table 2 keeps its MSE and coverage cells and adds each design's MSE ratio to SRS with
  a 95% paired interval per regime (t₄₉ over the 50 surface units, from `code/15`).
  Figure 2 gains per-point 95% Monte Carlo intervals (t₉ over the ten surfaces) and new
  shapes. Figure 1 may gain yearly maps (panel C) if legible.
- Methods order: outcome model, estimation and allocation risk; designs; simulation
  study; application. The grid's β = 1 makes incidence prognostic, the harder case for
  designs that use it; the application's education outcome uses β = 0 because it is not
  expected to depend on SUD incidence.
- The education outcome stays generic ("a continuous cluster-level education outcome").
  The supplement is "supplementary information (SI)". Funding uses Sage's standard
  no-funding sentence pending coauthor confirmation. The AI disclosure names Claude Code
  and Codex only.
- Main-text Results keep key results; run totals and replication counts move to the SI.
  Yearly budgets and population reach stay in the main text.

**Why.** "Pilot" made a finished study read as preliminary, and the full-surface run costs
about 1.5 hours. Readers need to see which designs differ significantly from SRS, and
designs share surfaces and noise, so only paired intervals answer that. Stating the
model before the simulation lets readers follow the setup (author review, 2026-10-27).

**Source.** `revision_markup/CTJ_Reading_Markup_V5.pdf`; plan
`~/.claude/plans/context-i-completed-a-wild-crab.md` (session 2026-10-27).

---

## 2026-10-25 (afternoon) — "Saturation Quadrants" restored; figure conventions; whole colleges

**Decision.** The plain saturation design is called "Saturation Quadrants" everywhere, in
text, tables and figures, reversing the 2026-10-21 "Saturation Regions". On the NC map the
chapter says that both saturation designs use four fixed geographic regions in place of the
quadrants. Manuscript figures carry no in-figure title, subtitle or note; the caption holds
them. Colleges needed for 80% power are reported as whole colleges, rounded up.

**Why.** Consistency between text and figures matters more than the name, and the NC
figures already said "Saturation Quadrants", so restoring it needed the fewest figure
regenerations. Whole colleges read more intuitively than interpolated means (author
feedback, 2026-10-25).

---

## 2026-10-25 — SRS in every design hypothesis test; figures before text; outline-only slides

**Decision.** Simple Random Sampling is a member of every hypothesis-test comparison
of designs and of its figures. The grid Friedman, Nemenyi and Holm-adjusted Wilcoxon
tests now cover the six designs plus SRS (`code/12`; the critical-difference diagrams
in `code/14`) and the eight designs plus SRS (`code/14`), and `code/15` adds paired
surface-level tests of each design against SRS within each spillover regime, drawn
as a forest plot. SRS stays a benchmark, not a proposed design, and remains a gray
reference line in the descriptive figures. Methods-extension step 7 builds the new
figures first, for review, before any text; step 8 drafts the slide outline (text plus
an HTML page with figures) in `presentation/` and builds no slides. The 2026-10-29
talk deck lives in `presentation/` and is adapted into prelim Part II later.

**Why.** SRS is the standard of care, so readers need to see whether each design
differs significantly from it, not only from the other designs (author feedback on
the methods-extension report, 2026-10-25). This reverses the 2026-09-26 rule that SRS
never enters the rank-based tests.

---

## 2026-10-26 — Application tables also use the t + df interval; both intervals reported

**Decision.** The application's coverage (the yearly table, Results text, CTJ abstract
and Results) now reports the t + df interval (94.1–94.9% by design and regime),
reversing the 2026-10-21 choice of Wald in the main tables. A new appendix/SI table,
`tab:interval-comparison`, gives Wald and t + df coverage side by side for the grid
(queen and rook) and the application (both regimes). Wald values stay stored
(production summaries; `*_Wald` columns; the comparison CSV).

**Why.** Author rule (2026-10-26): report the interval whose coverage is closest to
nominal, and keep both so either can be reported later. The application values come
from the same stored fits; `code/27_interval_comparison.R` reproduces the production
Wald yearly means exactly before writing the t + df values.

---

## 2026-10-26 — Grid intervals use the t + df correction

**Decision.** Grid coverage and power are reported with the same interval as the
application, τ̂ ± t(0.975, N − p) × SE × √(N/(N − p)), N = 100, p = coefficients + ρ
(oracle 5; non-oracle and rook Checkerboard, whose aliased spillover column is dropped,
4). Wald stays as a stated comparison (about one point lower). Scripts 12, 13, 14 and 19
call `apply_tdf_intervals()`; MSE, bias and every ranking are unchanged.

**Why.** Author: use the interval whose coverage is closest to nominal (2026-10-26).
The grid run had stored only summaries, so `code/26_grid_tdf_intervals.R` re-fit all
28,800 scenario-estimator rows from the same seed keys; every stored Wald metric was
reproduced exactly. Identified-design coverage moved from 0.936–0.942 to 0.946–0.952.

---

## 2026-10-26 — Habib's counts: "corrected", no personal communication

**Decision.** The CTJ, SI and chapter describe Habib's numerator as corrected counts
that include heart-failure deaths formerly excluded as presumed non-sudden, without a
"personal communication" citation. CTJ body now 3,477 words.

**Why.** Author instruction (2026-10-26); the correction is a property of the data
Habib prepared, already credited through his report and the author contributions.

---

## 2026-10-25 — Hybrid Table 1 with example maps; new design-sample figure

**Decision.** CTJ Table 1 gains an "Example (grid)" column with one example
assignment per design and an incidence-surface row; its rule text is unchanged and it
stays "Table 1" (the maps are part of the table), so the six-exhibit count holds. The
SI and Chapter 3 replace the old six-panel figure with a 2×4 panel figure (incidence
surface, six designs, SRS). Treated clusters are solid blue and controls white;
overlays show the 2×2 blocks, incidence quartiles and saturation quadrants with their
treated shares. All maps use one illustrative spatially autocorrelated surface
(ρ_X = 0.5, fixed seed), labeled as not one of the simulated surfaces. One script,
`code/25_design_assignment_figures.R`, draws every version, including a wide slide copy.

**Why.** The old figure overlaid incidence fill and ring/X markers, hid each design's
structure and used the retired "Block Stratified Sampling (Checkerboard)" label. The
independent-uniform surface showed little regional burden, so High Incidence Focus and
Incidence-Guided Saturation looked arbitrary. Visual examples beside the prose rules
explain the designs faster (author, 2026-10-25).

---

## 2026-10-21 — Habib affiliation, SUD case algorithm citation, IRB protocol

**Decision.** Habib's affiliation is the Department of Epidemiology, UNC Chapel Hill
(CTJ and Chapter 3 title blocks; authors now Walther¹, Habib², Simpson³, Lin¹). His
master's paper entry (`habib_temporal_2026`) keeps "Department of Biostatistics",
because it was written for the BIOS master's paper course. Cite Gan et al. (2019,
`gan_factors_2019`) as the source of the SUD case-adjudication algorithm Habib used,
in the CTJ, SI and Chapter 3 data paragraphs. State UNC IRB protocol 14-2036, "Sudden
Unexpected Death in North Carolina" (PI Ross J. Simpson, Jr., Department of Medicine)
in CTJ Methods and in the Chapter 3 application section (reversing "Chapter 3 has
neither statement").

**Why.** Author-supplied facts (2026-10-21). The CTJ no longer asserts a
not-human-subjects determination: only the protocol is confirmed, so the determination
and consent wording stay a highlighted placeholder until Simpson supplies the
disclosures and permissions. A dissertation chapter using restricted human data
should state its protocol.

---

## 2026-10-21 (afternoon) — Headline test and the plain saturation design's name

**Decision.**
- **Headline power uses the t + df-corrected test.** This supersedes the morning's
  "size-adjusted headline". Size-adjusted power stays as a check beside it, and Wald
  power goes to the SI.
- **The plain saturation design is renamed "Saturation Regions".** It replaces both
  "Plain Saturation" (AGENTS.md) and "Saturation Quadrants" (chapter, SI and older
  code) in every manuscript, slide and new figure. "Incidence-Guided Saturation" is
  unchanged.
  - The hashed production code (`code/03_designs.R`, `application/code/*`) and the
    stored CSVs keep their old labels. Rename at the display layer only: manuscript
    text, figure and table labels.

**Why.**
- On the A3 rerun the t + df test held 4.4–5.6% size at every N from 15 to 58, and it
  is the test a trial analyst can run. Size adjustment needs a simulated critical
  value.
- "Regions" describes both the grid quadrants and the NC k-means regions, and it
  removes the two-name conflict.

**Source.** Author review session, 2026-10-21.

---

## 2026-10-21 — Methods-extension results: decisions after review

**Decision.** Review of `docs/findings/2026-10-20-methods-extensions-findings.md`, one
decision at a time.
- **Targeting wording (D2 revised).** Incidence-Guided Saturation matches Plain
  Saturation's accuracy and sends more treatment to high-burden regions when burden
  differs between regions. On the NC map the four regions have similar mean incidence
  ranks (about 24–34 of 58 in every year), so the targeting gain there is small (rate
  share 0.511 vs 0.504). The per-college rate share is the primary targeting measure,
  with the death share beside it and one sentence on why they differ (high-rate
  colleges are small). A burden-share check on the stored grid surfaces tests whether
  the gain is larger when incidence clusters spatially.
- **Power headline.** Size-adjusted power, described as the power of a test with
  correct size. Wald power goes in an SI table with one sentence that Wald is
  anti-conservative below about 40 colleges. Also store and report power for the
  t + df-corrected test (an A3 rerun), so the trial-usable test sits beside it.
- **Story points adopted.** Under control-only spillover, power stops rising as
  nearly every college joins, and τ = 1 is not detectable with 58 colleges. The
  saturation designs' control-only advantage shrinks with fewer participants (22%
  lower MSE than SRS at N = 58, 4–5% at N = 29).
- **Coverage (A5).** Report the t + df-corrected interval: τ̂ ± t(0.975, n − p) ×
  SE × √(n/(n − p)). Stage 1 (all production fits) is the main evidence: 94.1–95.0%
  for every design × regime. The Stage 2 comparison (bootstrap, profile likelihood)
  goes in the SI. A reasoned exception to the 94–96% rule: three smaller Stage 2
  cells fall at 92.8–93.9%. No bootstrap rerun.
- **V6.** Rebuilt saturation regions must have at least 2 colleges, through a wrapper
  that replays the same 250 k-means candidates and skips smaller partitions;
  included in the A3 rerun.
- **V8.** No rerun. Report the dips along N with the subset-bootstrap intervals.
- **Bayesian CAR.** Rerun the ρ = 0.5 cells with chains four times longer.
- **Steps 7–8.** No unattended manuscript or slide edits. A report with proposed
  appendix/SI/main-text drafts and Part II slide ideas comes first, plus a prelim
  slide outline with an ideas bank.

**Why.** The NC regions balance incidence, so targeting depends on the map, while
the guidance is meant for maps where burden is regional. Stage 1 is the largest
evidence and the expensive intervals do no better. Size adjustment removes the small-N
Wald distortion, but a trial analyst needs a test they can compute.

**Source.** Author review session, 2026-10-21.

---

## 2026-10-20 — Methods extensions after the 2026-10-14 advisor meeting

**Decision.** Plan: `docs/plans/2026-10-20-methods-extensions-plan.md` (D1–D9).
- **Outcome (reaffirms the 2026-10-02 rule, with reasoning).** The trial's outcome Y is
  the measured education endpoint (an assessment of the trained first responders and
  healthcare professionals, averaged per college). It is an intermediate endpoint,
  hoped but not assumed to lower future SUD incidence; it is not future incidence.
  β = 0 stays primary. β = 1 stays a hypothetical sensitivity analysis, reported side
  by side: the endpoint tracks prior burden in both arms, which is not a claim that the
  education works better in high-incidence areas. τ is constant (no Z·X term). No
  future-incidence (count / rate-ratio) endpoint model; the manuscripts call it future
  work. Revisit only if Dr. Lin or Dr. Simpson ask.
- **Targeting.** Incidence guidance is justified by sending treatment to high-burden
  areas at plain saturation's accuracy (ethical rationale). Report the share of SUD
  burden in treated areas by design.
- **Partial participation (NC only, not the grid).** Colleges volunteer after the
  pilot. Primary: random volunteering. Sensitivity: the chance of joining rises with
  incidence rank (top about 3× bottom). Spillover exposure = treated neighbors / all
  geographic queen neighbors (non-participants count as untreated). The spatial lag
  links participants only; islands get zero rows. Designs are rebuilt on the
  participants. Full 58 stays primary.
- **Power.** Standardized effects. Fixed trained professionals per college (m) as the
  main case, translated to colleges × m × ICC tables; m proportional to population as
  a sensitivity. Fast `fit_sar_lag()` fitter, run overnight locally.
- **Tests.** Grid tests stay in the appendix/SI with one main-text sentence. NC tests
  are added and labeled a consistency summary, since the four years are related
  observed surfaces. Colleges needed for 80% power is the headline design comparison.
- **Sensitivities.** Coverage: try the t value, df-corrected variance, parametric
  bootstrap (percentile and studentized) and profile likelihood; adopt one only if
  coverage is 94–96% for every design and regime with rankings unchanged, otherwise
  state the limitation and what was tried. Estimation: CAR and Bayesian CAR
  (`CARBayes` install approved). Spillover forms: binary and additive count besides
  proportional; no second-order spillover.
- **Placement.** CTJ: power in main-text prose, with the figure in the SI (six main
  exhibits unchanged). Anything a committee member could ask about gets a main or
  backup prelim slide.

**Why.** Dr. Lin's questions (outcome versus incidence, partial participation, sample
size, tests, estimation method). The project recommends a design for a hypothetical
statewide trial with no future-incidence data, so the measurable endpoint is the
education assessment. Monte Carlo replicates make any difference "significant", so the
trial-relevant comparison is measured in colleges.

**Source.** Author planning interview, 2026-10-19/20.

---

## 2026-10-19 — Revisions apply to every manuscript version; versioned markups

**Decision.** Every revision to the CTJ paper is applied to both the submission layout
(`CTJ_Manuscript`) and the reading layout (`CTJ_Reading`), and both PDFs are rebuilt.
Whenever either the CTJ paper or Chapter 3 is revised, check whether the other version
needs the same edit. The author marks up a copy of the reading PDF saved as
`manuscript_ctj/revision_markup/CTJ_Reading_Markup_V<n>.pdf` (gitignored), and the clean
`CTJ_Reading.pdf` stays committed and current for sharing with collaborators.

**Why.** Content edits often matter to both the paper and the chapter, while edits tied to
the journal's structure may not; checking each time keeps the versions consistent. Markups
are personal working copies referenced in revision requests, so a new numbered copy is
made for each round of review.

Same day: the ethics placeholder now reads "[Placeholder for the UNC IRB determination:
protocol number and date]", and the AI acknowledgement says "testing code" and
"reviewing manuscript text" (author wording). Chapter 3 has neither statement.

---

## 2026-10-14 — CTJ ethics and coauthor-approval placeholders

**Decision.** Add a visible ethics placeholder at the end of the application
Methods, assuming a UNC IRB not-human-subjects determination (deidentified data on
deceased individuals), with the IRB number and date left as `[IRB #XX-XXXX, date]`.
Add "All authors critically revised the manuscript and approved the final version."
to Author contributions, followed by "[Placeholder: pending coauthor approval.]".

Add an Acknowledgements section with the Sage-required generative-AI disclosure
(draft wording: Claude Code and Codex assisted with code, code checking, figures and
text editing; the authors reviewed and verified all output).

**Why.** The author judged the 2026-10-12 removal of the ethics placeholder a mistake:
the draft shared with Lin should show the required statements. The ethics and approval
statements remain placeholders until Habib/Simpson confirm the determination and all
coauthors approve. Still never invent a protocol number. Both are highlighted yellow
(`\hl`) in the PDFs until finalized. The 3,500-word limit is taken to apply to the body
(declarations excluded; author, 2026-10-14).

**Replaces.** The 2026-10-12 removal of the in-text ethics placeholder (see the
2026-10-02 CTJ author declarations entry).

---

## 2026-10-14 — Best designs by scenario, not one recommended design

**Decision.** The chapter and manuscripts report the best designs by scenario,
mainly by spillover regime, instead of naming one design for all settings.

**Why.** The rankings change with the spillover regime. Checkerboard and 2×2
Blocking are near the best designs under both-arms spillover and the worst two
under control-only spillover, while Incidence-Guided Saturation stays first or
second under both. A trial team may not know the regime in advance, so an
answer given for each scenario is more robust than a single recommendation.

**Replaces.** The "one recommended design" in
`docs/plans/2026-09-24-simulation-revision-plan.md` (left as a dated record) and in the
Phase D ROADMAP item.

**Source.** `paper/thesis_chapter/Thesis_Chapter.qmd` (Discussion, guidance by
spillover scenario, added 2026-10-11).

---

## 2026-10-12 — Chapter 3 length-review targets

**Decision.** A body of 36–38 pp is an acceptable target, with Chapter 2's
length or shorter as the guide. Keep any result with real insight; keep the
error-vs-τ appendix unchanged; merge the annual NC tables (every year, every
design); keep the NC maps at full size; leave the SI unchanged except the
merged table. The concise copy was adopted 2026-10-13.

**Why.** To make Chapter 3 shorter and easier to follow without losing results
the author needs; cut text is kept in a side file.

**Source.** `docs/plans/2026-10-12-chapter3-length-review-plan.md` (status and
"Goal"); ROADMAP.md "Concise Chapter 3 adopted (2026-10-13)".

---

## 2026-10-12 — CTJ submission declarations and manuscript conventions

**Decision.** Data availability: simulation data and code on request; NC
death-certificate data restricted; no repository link. The
generative-assistance acknowledgement was removed from the draft; an AI
disclosure (required by Sage policy) will be written before submission. The
Habib citation for the corrected counts was left unchanged. The manuscript
conventions (abbreviations, one design-name set, incidence-guidance framing,
no repository hyperlink, no Codex/Claude mentions) date from the same review.

**Why.** Not recorded in the source, except that the death-certificate data are
restricted and, for the Habib citation, that Habib produced the corrected counts
after amending the adjudication method.

**Source.** `docs/reviews/2026-10-12-ctj-comprehensive-review-findings.md`
("Author decisions and resolution"); `paper/ctj_manuscript/submission_checklist.md`;
`paper/thesis_chapter/manuscript-declarations.md` ("Generative assistance");
AGENTS.md "Standing Author Rules" (manuscript conventions).

---

## 2026-10-02 — Application outcome: no SUD baseline term in the primary model

**Decision.** Primary continuous
education outcomes have no SUD-incidence baseline term (β = 0); observed rates
inform allocation. An explicitly hypothetical β = 1 baseline sensitivity is
accepted with rank-scaled X and an X-adjusted fit; primary ML fits treatment and
spillover without X. τ is constant within
a scenario in both cases. This supersedes older matched-primary β = 1 proposals;
the completed grid study remains unchanged.

**Source.** Author clarification, 2026-10-02 (moved from AGENTS.md 2026-10-14).

---

## 2026-10-02 — Application implementation approved

**Decision.** Implementation approved October 2. Continue the
approved plan without re-interviewing settled choices. The companion is
`application/report/real_sud_companion.html`; computation lives in
`application/code/real_sud_{setup,simulation}.R` and `run_real_sud.R`.
Follow `application/README.md` and the implementation findings for run status;
smoke/pilot outputs are verification stages, not full-study conclusions.
The full chapter/appendix, CTJ/SI and comprehensive review (2026-10-12) are done.
Old deadline dates no longer govern work.

Application design authority: start with
`docs/plans/2026-10-02-nc-sud-application-plan.md`; its status records settled
versus proposed choices. The integrated prompt is
`docs/prompts/2026-10-02-continuation-prompt.md`. Do not treat proposals as approved
model/design changes.

**Source.** Moved from AGENTS.md "Standing Author Rules" 2026-10-14.

---

## 2026-10-02 — Checkerboard stays in the main comparison

**Decision.** Retain Checkerboard in
Project 2's main design comparison, with detailed identification diagnostics in
the appendix. A brief Chapter 2 link motivates the extension; no direct numerical
cross-project comparison is required. CTJ must cite the verified Project 1
publication to establish the extension while remaining self-contained.

**Source.** Author-confirmed 2026-10-02 (moved from AGENTS.md 2026-10-14).

---

## 2026-10-02 — CTJ author declarations

**Decision.** Retain Walther, Habib, Simpson and Lin in that order,
with Lin corresponding. Author confirmed no current funding or conflicts.
The in-text ethics placeholder was removed (author, 2026-10-12); the IRB statement
is a TODO before submission. Never invent a protocol number or claim approval/exemption.

**Source.** Moved from AGENTS.md 2026-10-14.

---

## 2026-10-02 — Application uses observed incidence and a continuous SAR outcome; completion endpoint

**Decision.** The NC application uses each year's observed 2018–2021 SUD
incidence as fixed input (no synthetic, resampled or latent incidence), with the
four yearly surfaces reported separately. Outcomes come from a continuous
cluster-level SAR model with τ = 1 for the primary comparison. Also approved:
the regional-summary rule for Incidence-Guided Saturation (mean of average rate
ranks per region), treating 14 and 15 in the two Balanced Halves strata, and
keeping the grid allocation-risk pilot as supporting evidence with the NC
application as the next confirmation. Endpoint: complete the application, put
its results in full chapter/appendix drafts, then derive the CTJ manuscript/SI;
a proposal-only application section is not the final deliverable.

**Why.** There are no observed intervention outcomes, so outcomes are simulated
under known parameters to assess estimation error. The Balanced Halves change
keeps exactly half of the 58 clusters treated. Other reasons are not recorded.

**Source.** `docs/plans/2026-10-02-nc-sud-application-plan.md` (author
clarification, author-confirmed endpoint, "Author decisions to settle",
"Recommended primary analysis"); `docs/prompts/2026-10-02-continuation-prompt.md`
(author-confirmed completion goal).

---

## 2026-09-27 — How CTJ and the chapter refer to Project 1

**Decision.** Dissertation prose may refer directly
to "Chapter 2". CTJ must use self-contained wording (e.g., "In previous work...")
with a citation to the Project 1 BMC Medical Research Methodology paper, never
"Chapter 2" or "Project 1" as reader-facing references. Describe findings within
their studied conditions; verify final bibliographic details before adding the citation.

**Note (2026-10-14).** The Project 1 paper is now published (doi:10.1186/s12874-026-02968-0).

**Source.** Author-confirmed 2026-09-27 (moved from AGENTS.md 2026-10-14).

---

## 2026-09-27 — SRS as an informative benchmark for a conditional recommendation

**Decision.** SRS serves as an informative benchmark against which the
systematic designs are judged. Recommendations stay evidence-based and
conditional: state the conditions and criterion under which a design improves
on SRS, report SRS's real strengths, and do not discount the benchmark to reach
a preferred conclusion. The allocation-risk pilot scope was not approved by
this clarification.

**Why.** The author did not want readers to infer from average MSE alone that
allocation strategy is inconsequential. At queen/τ = 1 the
Incidence-Guided-to-SRS MSE ratio was about 0.63 under control-only spillover
and 1.08 under both-arms spillover, which supports a control-only benefit, not
universal superiority.

**Source.** `docs/findings/2026-09-27-ch2-check-findings.md` ("Author
clarification: SRS benchmark and recommendation (2026-09-27)").

---

## 2026-09-26 — SRS is a naive benchmark, not a proposed design

**Decision.** Add Simple Random Sampling (Design 9) as a benchmark: complete
randomization with exactly N/2 treated, ignoring incidence and geography. It is
a naive benchmark, not a proposed design. Its Z seed key includes d = 9, so the
other designs' draws are unchanged.

**Why.** There was no naive baseline: Checkerboard continues Chapter 2 but is
not a generic benchmark (2026-09-25 TODO).

**Source.** ROADMAP.md (the 2026-09-25 "TO DO before the CTJ" item);
`docs/plans/2026-09-24-simulation-revision-spec.md` §3 (design table and note);
`code/03_designs.R` (Design 9 comment).

---

## 2026-09-25 — Chapter first; Thesis_Chapter.qmd is the only edited file

**Decision.** The long-form chapter is the primary source (set in the 2026-09-18
plan); the CTJ manuscript
and SI are condensed from it by hand. The SpatialCRT chapter `.qmd` (now
`paper/thesis_chapter/Thesis_Chapter.qmd`, formerly
`Dissertation_Chapter.qmd`) is the only file anyone edits; the bios-dissertation
Chapter 3 copy is generated by the sync script.

**Why.** The chapter draft was already close to a full argument and a better
starting point than the ~8-page CTJ draft. Project 1's CTJ-first order was an
exception because its chapter came from an accepted paper.

**Replaces.** An early draft of the 2026-09-18 plan that made the CTJ the
primary source (corrected in the same session).

**Source.** `docs/plans/2026-09-18-manuscript-unification-plan.md` (header and
"Decision:", ~line 154); AGENTS.md "Dissertation Chapter 3 Sync"; ROADMAP.md
~line 733 (step 7 / sync).

---

## 2026-09-25 — Timing: chapter and CTJ as soon as possible

**Decision.** Both the chapter and the CTJ submission are wanted as soon as
possible. Step 5 (CTJ derivation) starts once steps 3, 7 and 4 are done, not
after the prelim deadline. Only the CTJ application section waits on the
real-data results. The chapter-first order stands.

**Why.** Not recorded in the source.

**Note.** On 2026-10-02 the author said old deadline dates no longer govern work.

**Source.** `docs/plans/2026-09-18-manuscript-unification-plan.md` ("Timing revised
(author, 2026-09-25)", ~line 484).

---

## 2026-09-25 — SUD numerator: Habib's corrected counts (23,523)

**Decision.** The application numerator is `final_county_sudden$num_obs`
(23,523 deaths, 2018–2021). The 21,147 counts behind `habib_temporal_2026` are
kept only as the reconciliation column `deaths_habib2026`.

**Why.** Habib's corrected case filtering no longer excludes heart-failure
deaths as presumed non-sudden, because adjudicated sudden cardiac death overlaps
non-negligibly with heart-failure patients (Habib, personal communication,
2026-09-25).

**Replaces.** The earlier choice of the 21,147 counts.

**Source.** `application/AGENTS.md` ("Real NC SUD data", ~line 46);
`application/README.md` ("Numerator decision", ~line 127).

---

## 2026-09-25 — Legal county boundaries; neighbors across water

**Decision.** Application contiguity weights use full legal county boundaries
(`cb = FALSE`), so counties that meet across water are neighbors. Compared with
the shoreline-clipped boundaries (`cb = TRUE`) of the earlier synthetic runs,
this adds 3 pairs: College of The Albemarle–Martin CC, Carteret–Pamlico and
Carteret–Beaufort County CC.

**Why.** An education intervention can plausibly benefit nearby counties
across a river or sound.

**Source.** `application/AGENTS.md` ("Contiguity weights", ~line 72);
`application/README.md` ("Boundaries", ~line 204).

---

## 2026-09-24 — Fix the simulation (M1–M8) and re-run rather than patch

**Decision.** Fix the simulation and re-run it before any committee or journal
review and before rewriting Methods/Results. Fixes M1–M8 include matched
incidence surfaces, 100,000 people per Poisson cluster, random tie-breaking,
key-based seeds, one noise column per fit, flagged aliasing, and a non-oracle
estimator as a sensitivity analysis. All April 2026 numbers are superseded.

**Why.** The step 0 chapter review found results that rest on simulation
oversights, not wording: designs and outcomes used unrelated incidence
surfaces, Poisson clusters were too small, seeding was not reproducible,
iteration counts overstated precision, aliasing was silent, and stale
checkpoints loaded silently.

**Replaces.** The April 2026 runs (archived in
`results/archive/pre_revision_20260924/`).

**Source.** `docs/reviews/2026-09-23-step0-review-findings.md` ("Decisions
(author, 2026-09-23/24)"); `docs/plans/2026-09-24-simulation-revision-plan.md` ("Context",
"Decisions recorded"); `docs/plans/2026-09-24-simulation-revision-spec.md`.

---

## 2026-09-24 — Aliased Checkerboard × rook estimates kept and flagged (M6)

**Decision.** Keep the aliased estimates (Checkerboard × rook, where WZ = 1 − Z
and the fit returns τ − γ) and flag them in exhibits (`N_Aliased`,
`Z_WZ_rank_deficient`). Explain the anomaly in Results/Discussion.

**Why.** So all six designs appear in every ranking; the anomaly connects to
Chapter 2's checkerboard paradox.

**Source.** `docs/plans/2026-09-24-simulation-revision-plan.md` (decisions table, M6;
B2 M6/M7, ~line 138); `docs/plans/2026-09-24-simulation-revision-spec.md` (Flags, ~line 150).

---

## 2026-09-24 — Balanced Quartiles and High Incidence Focus treat exactly N/2

**Decision.** High Incidence Focus treats exactly N/2 (part of fix M3). After the
B4 pilot, the author also decided that Balanced Quartiles treats exactly N/2:
floor(25/2) = 12 per quartile, plus one more in 2 randomly chosen quartiles.

**Why.** Before the revision, High Incidence Focus (`> median`) treated 14–50 on
tied Poisson incidence, and Balanced Quartiles treated 48 (`round(25/2)` = 12
per stratum).

**Source.** `docs/plans/2026-09-24-simulation-revision-spec.md` §3 (~line 74);
`docs/plans/2026-09-24-simulation-revision-plan.md` ("Context" item 2; B2 M3);
`code/03_designs.R` (Design 6 comment).

---

## 2026-09-24 — Six-design set kept for detailed grid rankings

**History (author, 2026-10-14).** 6 designs at first, expanded to 8, cut back to 6 because some were similar, then the SRS benchmark added. Now the CTJ main text has 6 + SRS and the thesis chapter 8 + SRS (6 in the detailed grid rankings).

**Decision.** Keep the six-design set for detailed grid rankings; Saturation
Quadrants and Balanced Halves stay out of them. The full benchmark and
application cover eight strategies plus SRS.

**Why.** After the re-run, Saturation Quadrants and Balanced Halves remained
statistically indistinguishable from their retained counterparts.

**Source.** ROADMAP.md "Prior State (as of 2026-09-24)" (~line 815); AGENTS.md
"Research Focus and Framing".

---

## 2026-09-24 — Research framing: compare designs

**Decision.** The research question is which designs estimate τ well under
heterogeneous incidence (low MSE, valid coverage). It is not which estimator is
better, and not whether knowing prior incidence is necessary.

**Why.** Not recorded in the source.

**Source.** AGENTS.md "Research Focus and Framing"; ROADMAP.md "Prior State
(as of 2026-09-24)".

---

## 2026-09-24 — Queen primary, rook sensitivity, never pool silently

**Decision.** Queen contiguity is primary; rook is reported separately as the
Chapter 2 continuity case and a sensitivity analysis. Incidence modes and
neighbor types are never pooled for headlines; pooled numbers appear only when
labeled, with splits wherever a conclusion changes. The CTJ main text leads with
the most application-realistic mode.

**Why.** The step 0 review found pooled averages that hid the split: a 52.7%
coverage figure averaged queen (~0.92) and rook (~0.13). Why queen rather than
rook is primary is not recorded.

**Source.** `docs/reviews/2026-09-23-step0-review-findings.md` (findings 26–27;
"Decisions"); AGENTS.md Critical Invariants #2; README.md "Key Design
Decisions" #4.

---

## 2026-09-24 — Design 1 full name: "Block Stratified Sampling (Checkerboard)"

**Decision.** `DESIGN_FULL_NAMES` labels Design 1 "Block Stratified Sampling
(Checkerboard)" in the design-sample figures. Elsewhere it is "Checkerboard";
the chapter states once that the two names are used interchangeably.

**Why.** It combines Chapter 2's name with this study's name.

**Source.** `docs/TODO.md` ("Known minor code issues").

---

## 2026-09-24 — One self-contained bibliography

**Decision.** Use one bibliography file in this repo,
`paper/SpatialCRT_IncidenceDesign.bib` (renamed from `IncidenceDesign_shared.bib`,
three diverging entries corrected from the master). Step 7 checks agreement with
the master entry by entry before porting.

**Why.** Not recorded in the source.

**Replaces.** The plan's pointer to bios-dissertation's `prelim/references.bib`.

**Source.** `docs/plans/2026-09-18-manuscript-unification-plan.md` ("Step 0 progress
(2026-09-24)", ~line 254); `docs/plans/2026-09-24-simulation-revision-plan.md` (decisions
table, Bib).

---

## 2026-07-02 — No DIM-vs-MLE comparison in manuscript prose

**Decision.** Manuscript prose (CTJ and chapter) never compares DIM with MLE and
does not mention DIM; ML estimation is presented as the estimator. Code and
reports may still mention DIM as a baseline.

**Why.** DIM was an internal validation baseline, never a scientific comparator.

**Source.** AGENTS.md "Research Focus and Framing"; date from the Claude memory
note `feedback_no_dim_comparison_no_shorthand` (manuscript-planning session,
2026-07-02).
