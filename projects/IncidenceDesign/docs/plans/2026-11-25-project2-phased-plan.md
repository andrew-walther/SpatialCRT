# Project 2 after the prelim submission: phased plan

Author: Andrew Walther. Written 2026-11-25. This is the master checklist for the Project 2 work
before the oral prelim (2026-12-07) and the CTJ submission. Work sequentially, one item at a time,
checking each off before the next (author, 2026-11-25); only long-running code or simulations
run in the background while the next item proceeds. Each phase runs in a fresh session started from this file, and check items off here as they finish. Detail for
Phases 1 to 3 is in `docs/plans/2026-11-25-morning-review-guide.md`; the overall brief is
`docs/prompts/2026-11-24-project2-post-prelim-prompt.md`.

**Notation reminder (since 2026-11-15, DECISIONS).** The documents write β for the intervention
effect on **Z** (swept over 0.8, 1, 1.5, 2, 3 on the grid) and τ for the incidence coefficient on
**X** (pre-specified: τ = 1 on the grid, τ = 0 in the NC case study with τ = 1 as a sensitivity).
ψ is the spillover magnitude. The code still uses the older names (`tau` on Z, `beta` on X,
`gamma` for spillover).

**Standing priorities for the CTJ (author, 2026-11-25).** First, Dr. Lin's feedback is fully
addressed. Second, wording and clarity. Third, length. Keep the key message in the main text so
readers do not have to look far; move to the SI only content that is not essential to it.

**Model use in every session brief (author, 2026-11-25).** Opus 5.5 where reasoning matters most:
drafting sign-off wording, structure and design decisions, statistical judgment, and the
independent reviews of main deliverables. Sonnet 5.5 subagents (`model: "sonnet"`) for lower-reasoning
work: reading and extracting text, targeted edits applied from an approved list, builds, text and
PDF comparisons, searches, table syncs and checker runs. Check a Sonnet result before relying on it.

**Codex prompts to paste (whenever convenient; the session running this plan reminds you at each
phase start until they are started).**

| Prompt | Runs | Status |
|---|---|---|
| `docs/prompts/2026-11-25-incidence-coefficient-codex-prompt.md` | Phase 4 (τ grid and estimation), in a Codex session started in this folder | [x] pasted 2026-11-25 |
| `docs/prompts/2026-11-25-igs-sensitivity-rerun-codex-prompt.md` | Phase 3a rerun and 3b figure, in a Codex session started in this folder | [x] pasted 2026-11-25 |

When Andrew says a Codex run has started or pastes its report or commit message, check the
named commit and outputs, tick the matching items, and log them here, so the record stays in one
place.

## Phase 0: setup and version record

- [x] Version record for the CTJ (tags, plus one clean PDF for the version Dr. Lin annotated):

  | Version | What it is | Body words | Where |
  |---|---|---|---|
  | V5 markup | Your markup of the 10-25 reading draft | — | `revision_markup/CTJ_Reading_Markup_V5.pdf` |
  | Pre-FCL | The draft Dr. Lin annotated (after your V5 revisions) | 3,498 | commit `8e36959` (2026-11-04), tag `ctj-pre-fcl`; his annotated copy is `revisions/2026-11-10-paper-2-revisions/2026-11-11-Project-2-Manuscript-Feedback-FCL.pdf` |
  | V6 | After Dr. Lin's feedback was applied (too long) | 3,853 | `revision_markup/CTJ_Reading_Markup_V6.pdf`, commit `4386bcb`, tag `ctj-v6` |
  | V7 | After round 6 (cuts and wording) | ≤ 3,500 | made at the end of Phase 1 |

  Done 2026-11-25. Dr. Lin's annotated PDF matches the `CTJ_Reading.pdf` committed at `8e36959`
  word for word; `a685b4a` (2026-11-10, first proposed here) came after it and only rounds MSE to
  three decimals and bolds matrices (same 3,498 words). Author decision: tag `8e36959` as
  `ctj-pre-fcl` and use it as the 1h diff base; tag `4386bcb` as `ctj-v6`; both tags pushed. The
  clean copy is that committed PDF, extracted unchanged (a rebuild would pick up today's
  figures): `revisions/2026-11-10-paper-2-revisions/CTJ_Reading_preFCL_clean.pdf`.
- [x] Push the SpatialCRT commits and the `prelim-submitted-2026-11-15` tags. Checked 2026-11-25:
  already done. SpatialCRT, bios-dissertation and multiomicsGEP were level with origin, and the
  tag was on the remote in all three.
- [x] Phase 4 runs on the side in Codex (table above); started 2026-11-25.
- Context, no action: the Chapter 3 sync hook keeps writing to `bios-dissertation/prelim/` until
  after the exam (2026-12-07); `thesis/` takes over later. Dissertation dates: BIOS abstract
  2027-03-04, thesis 2027-03-18, final exam 2027-04-09 (backup 04-07), ETD 2027-04-14.

## Phase 1: CTJ round 6 (Dr. Lin's feedback, wording, then length)

- [x] 1a. Read the audit (`docs/reviews/2026-11-24-ctj-revision-audit.md`) and the word-growth
  map. Decide the open feedback items: a reference for the Wald interval (FCL-024), and anything
  in the audit you want handled differently. Also decide what remains of
  `docs/plans/2026-10-22-ctj-structure-and-prose-review-plan.md` (not started): fold its open
  items into the 1c list, or drop it and say so in the plan's status.
- [x] 1b. Cut list (`revisions/2026-11-24-round6-word-limit/signoff.md`): accept, reject or edit
  each ID. Recommended set 3,385 words; about 115 words of slack. **Stop for decisions.**
- [x] 1c. Human-writer list (`human-writer-review.md`, HW-01 to HW-51): accept, reject or edit.
  Note HW-28 (β = 2 values are 0.1596 and 0.1604, so "nearly equal" stays). **Stop for decisions.**
- [x] 1d. Apply the accepted items (`cuts.py --apply`, then the HW edits), move cut text to the
  SI or `trimmed-material.md`, regenerate `Figure_Legends.md`, update the exhibit test if HW-28
  changes, rebuild both layouts, run `ctj_wordcount.py`, `test_manuscript_exhibits.R` and
  `verify_manuscripts.py`. Body ≤ 3,500. Also confirm the current Clinical Trials counting rule
  (what counts toward the 3,500) and the abstract limit (the checklist says 425 words) and check
  the abstract against it.
- [x] 1e. Chapter 3 carry-over as a second sign-off list. **Stop for decisions**, then apply and
  let the hook re-render the prelim copy.
- [x] 1f. Fix "application" in the cover letter and submission checklist (D-3).
- [x] 1g. You read the rebuilt `CTJ_Reading.pdf`. When you approve it, copy it as
  `CTJ_Reading_Markup_V7.pdf` and commit.
- [x] 1h. Review package for Dr. Lin, so he can review V7 without rereading the whole paper
  (author request, 2026-11-25). In `revisions/2026-11-24-round6-word-limit/for-lin/`:
  1. the updated manuscript (`CTJ_Reading.pdf`, V7);
  2. a marked-up diff from the draft he annotated (`8e36959`, tag `ctj-pre-fcl`) to V7, made with `latexdiff`
     (available in TinyTeX: `~/Library/TinyTeX/bin/universal-darwin/latexdiff`) on
     `CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`, built in the reading layout; check that
     tables, captions and math survive the diff markup (fall back to `--math-markup=0` or a
     per-section diff if a build breaks);
  3. a revision report (PDF, from a short Markdown or LaTeX source) that walks through the
     changes: for each of his comments (FCL IDs, in manuscript order) what he asked, what changed
     and where (section and page in V7), with the audit's "partial" items explained; then the
     additional revisions from further review (framed that way, not as committee questions,
     which he has no context for); then the length cuts, grouped by section, saying where each removed
     passage now lives (SI section or `trimmed-material.md`); and the word count (3,498 → 3,853
     → V7). Built only from what was approved and applied in V7 (the
     accepted cut, HW and audit decisions, checked against the diff), never from proposals that
     were rejected or edited, and never retyped. Human-writer pass. **Stop for your read before anything is sent**; you send it.

## Phase 2: SI curation (after Phase 1, because the cuts move text into the SI)

Decided 2026-11-25: one source with a switch. `Supplementary_Information.tex` builds the curated
SI that is submitted; `Supplementary_Information_Full.tex` turns everything on, so nothing built
before or later is lost and anything can move between the two. Target for the curated SI: about
30 pages, with room for new content (the DeSurv SI is 37 pages in all).

- [x] 2a. Add the switch and the full-version wrapper; both build.
- [x] 2b. Link the CTJ's 12 hard-coded SI numbers to SI labels (`xr-hyper`) before anything moves.
- [x] 2c. Correctness and rule fixes from the SI review: the "recommended design" heading and
  wording, SRS in every design test, stale statements, "identified" → "identifiable" (FCL-028),
  figure-file issues (old labels, minus signs as hyphens), ratio decimals. **Decision on ratio
  decimals** (3 per the MSE rule, or 2 like main Table 2).
- [x] 2d. Bring the SI tables up to date with Chapter 3 (four changed tables, two missing). The
  verifier must pass.
- [x] 2e. Build the curated outline (methods before results, one numbered list); everything else
  stays in the full version. **Stop to review the curated table of contents.**
- [ ] 2f. Place the text moved out of the CTJ in Phase 1.
- [ ] 2g. Prose review: human-writer pass, then a statistician and a clinician reader on the
  curated SI.
- [ ] 2h. Rebuild both versions, run the checkers, record the page counts.
- [ ] 2i. Chapter 3 alignment: resolve `paper/ctj_manuscript/notes/si-ch3-alignment-log.md` (42
  SI edits not yet in Chapter 3) plus this phase's SI changes, as one sign-off list. **Stop for
  decisions**, then apply and let the hook re-render the prelim copy.

## Phase 3: Incidence-Guided Saturation sensitivity study

- [x] 3a. Review the spec and findings (`docs/plans/2026-11-24-igs-sensitivity-spec.md`,
  `docs/findings/2026-11-24-igs-sensitivity-findings.md`). Open choices:
  - Saturation sets. Each region gets its own level, evenly spaced from 50% + s down to 50% − s,
    so four regions at s = 0.3 get 80/60/40/20 and four regions at s = 0.2 get 70/57/43/30. The
    widest tested was s = 0.4 (90 to 10). s = 0.5 would give 100/67/33/0 for quadrants: one region
    fully treated and one fully control. Recommended (author leaning, 2026-11-25): add s = 0.5
    as a boundary case, to show whether going past 90/10 helps or hurts rather than to propose
    it. It is eight more variants (four region counts, guided and unguided), about 20 minutes.
    Expectation to check: under control-only spillover the gain from wider spreads levels off or
    reverses, because a fully treated region contributes no control clusters and a fully control
    region no treated ones; under both-arms spillover it gets worse, as the wider spreads already
    did.
  - Repeated levels for many regions (for example nine regions sharing the four quadrant levels)
    instead of one level per region. Add?
  - Run the North Carolina map too?
  **Decisions; rerun if anything changes (about 80 minutes).**
  Decided 2026-11-25: add s = 0.5 (8 variants, IDs 131–138); no repeated-level variants; grid
  only (NC later only if the findings call for it). The rerun goes to Codex (prompt in the table).
- [x] 3b. (Decided 2026-11-25: x-axis labels 60/40 … 100/0, keep the circles with a legend; Codex
  does it with the rerun.) Figure: a plainer x-axis (the highest and lowest saturation, for example "60/40" to
  "90/10", instead of the spread s), and either explain the circles (they mark the original design
  and Saturation Quadrants) in a legend and caption or drop them.
- [ ] 3c. SI section and figure (curated SI; both the method and the findings, from the findings
  report), Chapter 3 appendix mirror, and the optional CTJ
  sentence if the word count allows. **Stop for wording.**
- [ ] 3d. Update the Q-16 answer in `presentation/prelim-qa-prep.md` and make the backup-slide
  figure.

## Phase 4 (side, Codex): the incidence coefficient τ

Plan: `docs/plans/2026-11-25-incidence-coefficient-plan.md`. Kickoff prompt:
`docs/prompts/2026-11-25-incidence-coefficient-codex-prompt.md`.

What τ means: it is a regression coefficient, not a correlation. Incidence X is scaled to [0, 1]
on the grid, and the residual SD is 1, so τ = 1 means the highest-incidence cluster's
outcome sits one residual SD above the lowest one's (before the spatial multiplier), about the
size of the treatment effect β = 1. τ = 0 means no association. τ is not bounded by 1 and could be
negative. So far τ was only set in the data generation (1 on the grid; 0, with 1 as a sensitivity,
in the NC case study). The grid fit always estimated it freely, because X is a covariate, but the
estimate was never reported. In the NC case study X enters the fitted model only in the τ = 1
sensitivity (`rs_model()` in `application/code/real_sud_simulation.R`).

- [x] 4a. Codex writes the spec; you approve it. Proposed τ values {0, 0.25, 0.5, 1, 2}, about
  15 to 20 minutes each.
- [x] 4b. Grid at each τ, with τ = 1 reproducing production exactly. Expected: β̂ barely changes
  with τ, because both fits include X, so τX is absorbed except through the estimate of ρ. A
  near-null result here is itself the answer (the τ = 1 choice did not drive the grid
  conclusions).
- [x] 4c. Estimating τ when both prior incidence and trial outcomes are known: bias, coverage and
  power by design. Rough expectation, to be checked: the SE of τ̂ is near 0.35 for designs whose
  treatment is unrelated to incidence and about twice that for High Incidence Focus (its Z is
  X above the median, so X has little independent variation given Z). The extra τ values give a
  power curve for detecting τ.
- [x] 4d. Findings reviewed; you decide placement. Decided 2026-11-25: curated SI, mirrored in
  the Chapter 3 appendix; drafted in the Phase 2 session after 2e, with a CTJ pointer sentence of
  about 15 words if it fits (about 65 words of slack remain). This reopens the 2026-11-14 decision that
  estimating τ is future work; record the new decision in DECISIONS.md.
- [ ] 4e. SI Methods subsection (author request, 2026-11-25), drafted by the main session after 4d,
  not by Codex: "Estimating the incidence coefficient after a trial". What τ is, that after a trial
  the investigator has both prior incidence and outcomes so τ can be estimated within the same
  model, what it tells them (the association between prior incidence and the later outcome,
  adjusted for treatment and spillover), how the design affects that precision (with the 4c
  results), and the limits. Curated SI if the findings are useful, full SI otherwise; Chapter 3
  mirror. The SI holds the full discussion of τ (author, 2026-11-25); the CTJ gets at most a
  sentence pointing to it, if any. **Stop for wording.**

## Phase 5: prelim presentation (after Phases 1 to 3)

- [ ] 5a. Write the presentation prompt (include Q&A prep on τ: what the coefficient means on
  the X ∈ [0, 1], σ = 1 scale, that it is estimated in every grid fit, and the Phase 4 results): Parts I to III, the script, backup slides, the Chapter 2
  two-reviewer question review, rehearsal with Dr. Lin (2026-12-02 to 12-04) and Part III at the
  Rashid Lab meeting (2026-12-03).
  The prompt must also carry these author requirements:
  - Exam: 2026-12-07, a 2-hour block with about 60 minutes of talk and about 90 minutes of exam
    in all; Parts I to III, about 20 minutes each.
  - The presentation script is the engine: slides carry only key details, and each backup slide
    focuses on one result. The talk may use the current project state.
  - Part I: update and refine the prelim slides first; then make a standalone copy of the
    finished Part I for Dr. Lin's conference, with a new title and agenda.
  - Part II: keep a frozen copy of the deck as presented at the SUDDEN meeting
    (`presentation/sudden_talk_2026-10-29.qmd`); the prelim Part II is a separate copy that is
    updated and refined further.
  - Part III is presented at the Rashid Lab meeting on 2026-12-03.
  - Before the slides use them, settle the case-identification counts in
    `presentation/data/case_identification_counts.csv` (author: Gan's method, with Habib's counts
    later corrected by reclassifying cardiac cases) and cite the right source (TODO item).
  The prompt must cover a Project 1 question bank as well: the prelim review produced
  anticipated questions only for Project 2 (the Chapter 3 committee-style review C-01 to C-31 and
  the Q-01 onward list in `presentation/prelim-qa-prep.md`), so Chapter 2 gets its own
  two-reviewer question review, answers and backup slides.
- [ ] 5b. Run it.

## Log

Record each phase's completion here with the date and commit.

- 2026-11-25. Phase 0 done (`b0f7891`): tags `ctj-pre-fcl` (`8e36959`) and `ctj-v6`
  (`4386bcb`) pushed; clean pre-FCL PDF saved; earlier pushes verified. Codex Phase 4 prompt
  not yet pasted.
- 2026-11-25. Phase 4a–4c done in Codex (`3f92ffe`): spec, smoke, pilot and full sweep at τ ∈ {0,
  0.25, 0.5, 1, 2} (145.8 min). τ = 1 reproduces every stored production β̂ metric exactly (28,800
  scenario and 288,000 surface rows). Coordinator spot-checked the High Incidence Focus vs SRS τ̂ SE
  and power against `performance.csv`. 4d (placement) waits on the author.
- 2026-11-25. Phase 1 done (`568a565` to `ca20672`; V7 `e9fee02`, review package `ca20672`).
  Body 3,498 → 3,853 → 3,435 words; abstract 422 of 425. V7 copy matches the current
  `CTJ_Reading.pdf`; V1–V6 markups untouched; the diff base is `ctj-pre-fcl`. Exhibit test passes.
  `verify_manuscripts.py` still fails on SI-vs-Chapter 3 table drift (`tab:designs` now differs
  after the 1e carry-over), which Phase 2d resolves. The package waits on the author's read; he
  sends it.
- 2026-11-25. 4d decided (curated SI, Chapter 3 mirror; DECISIONS updated). 4e moves into the
  Phase 2 session.
- 2026-11-25. 2a done (`2b47ab4`): curated/full switch; both versions build at 65 pages, text
  identical to the previous SI. 2b (`3affdff`) checked: printed CTJ text unchanged, but the word
  counter dropped the 12 `\ref{SI-...}` numbers (3,435 → 3,423); fix requested before ticking.
- 2026-11-25. 2b done (`3affdff`, counter fix `18b81bf`): 13 SI references live; body back to
  3,435 words.
- 2026-11-25. 2c–2e done (`e031057` to `77839ec`): SI fixes, tables match Chapter 3, curated SI
  44 pages and unabridged SI 64 pages from one source (above the ~30-page target). CTJ still
  3,435 words with no unresolved SI references; exhibit test and verifier pass. Case-study
  production manifest re-recorded for the code/16 hash only (`93ffeb3`, DECISIONS); exhibits
  reproduce byte for byte. Next: Phase 2 Part B (2f, 4e, 2g–2i).
- 2026-11-25. 3a–3b done in Codex (`198f188`): s = 0.5 variants (IDs 131–138) added; variants
  101–130 and Designs 3, 8, 9 reproduce the 11-24 run, and the 09-27 production run, exactly.
  Figure relabeled 60/40 … 100/0 with a legend for the circles. At 100/0 the two-region version
  turns worse under control-only spillover and about 8 times SRS under both-arms spillover.
  Next: 3c (SI section) after Phase 2.
