# Project 2 after the prelim submission: phased plan

Author: Andrew Walther. Written 2026-11-25. This is the master checklist for the Project 2 work
before the oral prelim (2026-12-07) and the CTJ submission. Work one phase at a time, each in a
fresh session started from this file, and check items off here as they finish. Detail for
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

**Codex prompts to paste (whenever convenient; the session running this plan reminds you at each
phase start until they are started).**

| Prompt | Runs | Status |
|---|---|---|
| `docs/prompts/2026-11-25-incidence-coefficient-codex-prompt.md` | Phase 4 (τ grid and estimation), in a Codex session started in this folder | [ ] pasted |

## Phase 0: setup and version record

- [ ] Version record for the CTJ (tags, plus one clean PDF for the version Dr. Lin annotated):

  | Version | What it is | Body words | Where |
  |---|---|---|---|
  | V5 markup | Your markup of the 10-25 reading draft | — | `revision_markup/CTJ_Reading_Markup_V5.pdf` |
  | Pre-FCL | The draft Dr. Lin annotated (after your V5 revisions) | 3,498 | commit `a685b4a` (2026-11-10); his annotated copy is `revisions/2026-11-10-paper-2-revisions/2026-11-11-Project-2-Manuscript-Feedback-FCL.pdf` |
  | V6 | After Dr. Lin's feedback was applied (too long) | 3,853 | `revision_markup/CTJ_Reading_Markup_V6.pdf`, commit `4386bcb` |
  | V7 | After round 6 (cuts and wording) | ≤ 3,500 | made at the end of Phase 1 |

  Proposed: local tags `ctj-pre-fcl` (`a685b4a`) and `ctj-v6` (`4386bcb`), and a clean rebuild of
  the pre-FCL reading PDF saved as `revisions/2026-11-10-paper-2-revisions/CTJ_Reading_preFCL_clean.pdf`
  (the markup folder stays for your annotated copies). **Decision: yes or no.**
- [ ] Push the six unpushed SpatialCRT commits (overnight reviews, cut list, sensitivity study,
  this plan) and the `prelim-submitted-2026-11-15` tags. **Decision.**
- [ ] Phase 4 runs on the side in Codex (table above); start it whenever convenient.

## Phase 1: CTJ round 6 (Dr. Lin's feedback, wording, then length)

- [ ] 1a. Read the audit (`docs/reviews/2026-11-24-ctj-revision-audit.md`) and the word-growth
  map. Decide the open feedback items: a reference for the Wald interval (FCL-024), and anything
  in the audit you want handled differently.
- [ ] 1b. Cut list (`revisions/2026-11-24-round6-word-limit/signoff.md`): accept, reject or edit
  each ID. Recommended set 3,385 words; about 115 words of slack. **Stop for decisions.**
- [ ] 1c. Human-writer list (`human-writer-review.md`, HW-01 to HW-51): accept, reject or edit.
  Note HW-28 (β = 2 values are 0.1596 and 0.1604, so "nearly equal" stays). **Stop for decisions.**
- [ ] 1d. Apply the accepted items (`cuts.py --apply`, then the HW edits), move cut text to the
  SI or `trimmed-material.md`, regenerate `Figure_Legends.md`, update the exhibit test if HW-28
  changes, rebuild both layouts, run `ctj_wordcount.py`, `test_manuscript_exhibits.R` and
  `verify_manuscripts.py`. Body ≤ 3,500.
- [ ] 1e. Chapter 3 carry-over as a second sign-off list. **Stop for decisions**, then apply and
  let the hook re-render the prelim copy.
- [ ] 1f. Fix "application" in the cover letter and submission checklist (D-3).
- [ ] 1g. You read the rebuilt `CTJ_Reading.pdf`. When you approve it, copy it as
  `CTJ_Reading_Markup_V7.pdf` and commit.
- [ ] 1h. Review package for Dr. Lin, so he can review V7 without rereading the whole paper
  (author request, 2026-11-25). In `revisions/2026-11-24-round6-word-limit/for-lin/`:
  1. the updated manuscript (`CTJ_Reading.pdf`, V7);
  2. a marked-up diff from the draft he annotated (`a685b4a`) to V7, made with `latexdiff`
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

- [ ] 2a. Add the switch and the full-version wrapper; both build.
- [ ] 2b. Link the CTJ's 12 hard-coded SI numbers to SI labels (`xr-hyper`) before anything moves.
- [ ] 2c. Correctness and rule fixes from the SI review: the "recommended design" heading and
  wording, SRS in every design test, stale statements, "identified" → "identifiable" (FCL-028),
  figure-file issues (old labels, minus signs as hyphens), ratio decimals. **Decision on ratio
  decimals** (3 per the MSE rule, or 2 like main Table 2).
- [ ] 2d. Bring the SI tables up to date with Chapter 3 (four changed tables, two missing). The
  verifier must pass.
- [ ] 2e. Build the curated outline (methods before results, one numbered list); everything else
  stays in the full version. **Stop to review the curated table of contents.**
- [ ] 2f. Place the text moved out of the CTJ in Phase 1.
- [ ] 2g. Prose review: human-writer pass, then a statistician and a clinician reader on the
  curated SI.
- [ ] 2h. Rebuild both versions, run the checkers, record the page counts.

## Phase 3: Incidence-Guided Saturation sensitivity study

- [ ] 3a. Review the spec and findings (`docs/plans/2026-11-24-igs-sensitivity-spec.md`,
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
- [ ] 3b. Figure: a plainer x-axis (the highest and lowest saturation, for example "60/40" to
  "90/10", instead of the spread s), and either explain the circles (they mark the original design
  and Saturation Quadrants) in a legend and caption or drop them.
- [ ] 3c. SI section and figure (curated SI), Chapter 3 appendix mirror, and the optional CTJ
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

- [ ] 4a. Codex writes the spec; you approve it. Proposed τ values {0, 0.25, 0.5, 1, 2}, about
  15 to 20 minutes each.
- [ ] 4b. Grid at each τ, with τ = 1 reproducing production exactly. Expected: β̂ barely changes
  with τ, because both fits include X, so τX is absorbed except through the estimate of ρ. A
  near-null result here is itself the answer (the τ = 1 choice did not drive the grid
  conclusions).
- [ ] 4c. Estimating τ when both prior incidence and trial outcomes are known: bias, coverage and
  power by design. Rough expectation, to be checked: the SE of τ̂ is near 0.35 for designs whose
  treatment is unrelated to incidence and about twice that for High Incidence Focus (its Z is
  X above the median, so X has little independent variation given Z). The extra τ values give a
  power curve for detecting τ.
- [ ] 4d. Findings reviewed; you decide placement. This reopens the 2026-11-14 decision that
  estimating τ is future work; record the new decision in DECISIONS.md.
- [ ] 4e. SI Methods subsection (author request, 2026-11-25), drafted by the main session after 4d,
  not by Codex: "Estimating the incidence coefficient after a trial". What τ is, that after a trial
  the investigator has both prior incidence and outcomes so τ can be estimated within the same
  model, what it tells them (the association between prior incidence and the later outcome,
  adjusted for treatment and spillover), how the design affects that precision (with the 4c
  results), and the limits. Curated SI if the findings are useful, full SI otherwise; Chapter 3
  mirror. **Stop for wording.**

## Phase 5: prelim presentation (after Phases 1 to 3)

- [ ] 5a. Write the presentation prompt (include Q&A prep on τ: what the coefficient means on
  the X ∈ [0, 1], σ = 1 scale, that it is estimated in every grid fit, and the Phase 4 results): Parts I to III, the script, backup slides, the Chapter 2
  two-reviewer question review, rehearsal with Dr. Lin (2026-12-02 to 12-04) and Part III at the
  Rashid Lab meeting (2026-12-03).
- [ ] 5b. Run it.

## Log

Record each phase's completion here with the date and commit.
