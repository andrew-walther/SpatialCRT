# Prompt: Project 2 work after the prelim submission (CTJ, SI, sensitivity study)

Author: Andrew Walther. Written 2026-11-24; replaces the Project 2 parts of
`2026-11-19-post-prelim-planning-prompt.md`. The bios-dissertation organization has its own
prompt (`bios-dissertation/thesis/notes/2026-11-24-thesis-organization-prompt.md`). The
slides, the script, the Chapter 2 question review and the rehearsals come after this work
and get their own prompt.

Paste the fenced block into a session started in `SpatialCRT/projects/IncidenceDesign/`.

```text
Goal: bring Paper 2 (the Clinical Trials manuscript, CTJ) and its Supplementary Information
(SI) to a submittable state, and finish the Incidence-Guided Saturation sensitivity study,
before the prelim slides are built. The oral prelim is Monday 2026-12-07, with rehearsals
with Dr. Lin on 12-02 to 12-04, so this work should be finished, or cut back to fit, by about
2026-11-29. Estimate in real hours. Use the brainstorming skill to plan, then work
through the plan, stopping where noted.

Order (author, settled): (1) CTJ revision check and further revisions, then the word cut;
(2) SI review, curation and prose review; (3) the sensitivity study's results into the SI
and a backup slide. Step 3's compute already ran (see below), so its write-up can slot in
wherever it fits.

Read first: AGENTS.md (standing author rules, notation β on Z, τ on X, ψ for spillover;
design names; MSE to 3 decimals; t + df intervals with Wald beside), DECISIONS.md,
docs/TODO.md, ROADMAP.md, memory notes on CTJ revisions (edits go to CTJ_Manuscript.tex and
CTJ_Exhibits.tex, rebuild both CTJ_Manuscript.pdf and CTJ_Reading.pdf, check Chapter 3 both
ways; never touch revision_markup/).

Inputs prepared overnight on 2026-11-24 (read these before planning):
- docs/reviews/2026-11-24-ctj-revision-audit.md: every FCL item (Dr. Lin's 60 annotations)
  and every meeting item checked against the CTJ and Chapter 3. 81 done, 10 partial,
  0 missing or deviating. Partial items include the Wald citation (FCL-024), "identified"
  left in the SI (FCL-028), "application" left in the cover letter and checklist (D-3),
  ", so" joins (W-03), the five-configuration justification (S-03) and the binary-spillover
  mention (Q-01). Both checkers fail on stale references, not wrong numbers:
  paper/ctj_manuscript/figures/Figure_Legends.md still uses the old notation, the verifier
  looks for a removed label tab:nc-summary-sens, and test_manuscript_exhibits.R searches
  for τ = 2 and τ = 3 where the text now says β.
- docs/reviews/2026-11-24-si-claim-inventory.md and 2026-11-24-si-review.md: the SI section
  by section (what points to it from the CTJ, what it supports, length), a proposed abridged
  outline, where the full-length version should live, correctness and prose findings.
- docs/plans/2026-11-24-igs-sensitivity-spec.md (provisional spec), code/29_igs_sensitivity.R,
  code/29_igs_variant_designs.R, code/30_igs_sensitivity_summary.R, and
  results/igs_sensitivity_2026-11/full/ (run log, reproduction check, summary/ tables and
  figure).
- paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/word-limit-trim-plan.md (T1 to
  T14) and docs/plans/2026-10-22-ctj-structure-and-prose-review-plan.md (not started).
- CTJ_Reading_Markup_V6.pdf is a frozen copy taken after Dr. Lin's feedback was applied. It
  is not marked up; the TODO line asking for a markup of it is wrong.

Part 1: CTJ
1a. Fix the checkers first (Figure_Legends.md from the master captions, the removed label,
    the τ→β search strings) so every later edit is checked. These are tooling and legend
    fixes, not content.
1b. From the audit's partial items plus the structure and prose review plan, build one list
    of further revisions. From the trim plan, build a per-section word budget with a single
    scripted count (look up the current Clinical Trials counting rule and the abstract limit,
    425 words in the checklist; the trim plan's 3,838 counts math as one word, the
    checklist's older 3,498 excluded math, so they are not comparable). Target 3,500 with a
    margin. T1 reverses Dr. Lin's request FCL-032; flag it.
1c. Present both as one sign-off list (ID, source, location, current → proposed text, words
    saved), in the Chapter 1 round 5 format. STOP for Andrew's decisions. Apply the accepted
    items, rebuild both layouts, run the checkers, then propose the Chapter 3 carry-over as a
    second list. Cut text goes to the SI with a pointer or to
    paper/thesis_chapter/notes/trimmed-material.md, never silently deleted.
1d. After the cut, one reading-flow pass so the paper reads as one piece.
Submission blockers that do not gate this work (IRB determination and consent wording,
permission to publish the 58-college maps, coauthor approvals, Dr. Simpson's answers on
the control colleges and the pilot): keep them on one list in docs/TODO.md.

Part 2: SI
2a. From the claim inventory, decide with Andrew where the full-length SI lives (a new file
    or Chapter 3's Appendix B) and the abridged outline and target length. STOP for that
    decision.
2b. Build the abridged SI; keep the full version complete. Every CTJ pointer must land on
    the right abridged section; sections moved in by the CTJ cut (Part 1) must be there.
2c. Correctness fixes from the review, then the prose review (human-writer pass, then a
    statistician and a clinician reviewer reading the abridged SI cold). Resolve
    paper/ctj_manuscript/notes/si-ch3-alignment-log.md along the way.
2d. Run verify_manuscripts.py (with --pending-si for candidate copies) and
    test_manuscript_exhibits.R; rebuild the SI PDF.

Part 3: sensitivity study
3a. Review the spec's provisional choices with Andrew (regions 2, 4, 10, 25; spreads 0.1 to
    0.4; linear saturation sets; grid only; pooled over effect sizes) and the results. If he
    changes the design, rerun 29 (about an hour) and 30.
3b. Write it up: one SI subsection with one figure or table (the questions are whether the
    variants keep the original design's standing against SRS, and whether guidance still
    beats its unguided twin), one CTJ pointer sentence only if the word budget allows,
    a Chapter 3 appendix mirror, and a backup-slide figure for committee question Q-16
    (update its answer in presentation/prelim-qa-prep.md).
3c. Add the scripts to the AGENTS.md file map and the README.

Throughout: numbers come from saved results or executed code, never typed by hand; no DIM;
American English; avoid colons in prose. Commit after each logical step and push once
reviewed. Update ROADMAP.md, docs/TODO.md and AGENTS.md at the end.
```
