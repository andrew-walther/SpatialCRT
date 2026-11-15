# Chapter 3 (incidence design) review and CTJ revision: session prompt

Written 2026-11-14. Author: Andrew Walther. Paste the block below into a new session started in
`SpatialCRT/projects/IncidenceDesign`, with
`bios-dissertation/prelim/project-proposals/project2-incidence/draft/` added to the session (or
the reverse). It is
modeled on the Chapter 4 plan,
`multiomicsGEP/docs/plans/2026-11-14-chapter4-review-session-plan.md`, so the two reviews run
the same way.

---

```text
I'm reviewing Chapter 3 of my prelim (Project 2, sampling design for spatial cluster randomized
trials under heterogeneous incidence) and revising the CTJ manuscript from the feedback on its
draft. I want the insights from that feedback carried into Chapter 3. Work as described below.

## Read first

1. This project's AGENTS.md (standing author rules), DECISIONS.md, ROADMAP.md, docs/TODO.md.
2. bios-dissertation/AGENTS.md and ROADMAP.md (rendering, hooks, prelim status).
3. The Chapter 4 review plan, ../../../multiomicsGEP/docs/plans/2026-11-14-chapter4-review-session-plan.md
   (paths relative to this project). Use its structure (phases, stops, review output formats)
   and its general decisions. Its PDAC-specific content does not apply.
4. The meeting decisions and revision list, docs/meetings/2026-11-14-ctj-feedback-meeting-notes.md
   (decisions D-1 to D-5; revision items M-, S-, A-, I-, W-; open questions Q-). These come from
   my review meeting with Dr. Lin. Its decisions are settled. Its original notes are at the end
   of the file, and my notes govern where the Zoom summary differs.
5. The CTJ feedback register, paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/
   (README.md, feedback.json with FCL-001 to FCL-060, feedback.html). Ask me whether I have an
   exported notes JSON from feedback.html to add.

## Where to work and how edits are synced

Chapter 3, working location:
/Users/ajwalther/GithubProjects/bios-dissertation/prelim/project-proposals/project2-incidence/draft/
is where the chapter is read, rendered and checked.
- Make your edits in the source, paper/thesis_chapter/Thesis_Chapter.qmd in this project (and
  in the figure scripts for exhibits built by code). Commit and push in SpatialCRT as you work,
  after each approved group. The post-commit hook re-renders the draft and commits it in
  bios-dissertation. Check each change in the re-rendered draft.
- Never edit or commit the draft folder by hand: the bios-dissertation pre-commit guard refuses
  it, and the next sync would overwrite it.
- I will not edit the draft .tex by hand until the revamped chapter is complete. When I do, my
  edits are carried back with the prelim-edit-sync skill (prelim/full-draft/tools/prelim_review.py:
  start, diff, prompt, verify, finish ch3).

CTJ manuscript, working location:
/Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/paper/ctj_manuscript/manuscript_ctj/
- Edit CTJ_Manuscript.tex and CTJ_Exhibits.tex there (CTJ_Reading.tex only selects the reading
  layout). Rebuild both CTJ_Manuscript.pdf and CTJ_Reading.pdf after each edit. The supplement is
  ../supplementary_info/Supplementary_Information.tex. A hook copies the committed manuscript to
  bios-dissertation/manuscripts/project2-incidence/.
- When I approve a version as final for a revision round, copy CTJ_Reading.pdf to
  revision_markup/CTJ_Reading_Markup_V<n>.pdf, with n one more than the highest existing
  number (V1 to V5 exist now, so the next is V6). Never edit or overwrite an existing markup copy.
  The markup copies are tracked in git: commit and push the new copy, and tell me its name.

Both: render with RStudio's Quarto 1.10 (bios-dissertation/AGENTS.md), and after every render
compare page counts with the committed PDFs. After any edit to one of Chapter 3 or CTJ, say
whether the other needs it too.

## Context

- Project 2 is complete: the grid simulation, the North Carolina sudden unexpected death
  application and the methods extensions. Chapter 3 is the long-form canonical version, and CTJ
  is condensed from it by hand. Chapter 3 presents a finished study with a short Future
  Considerations. That differs from Chapter 4, which is a proposal with continuing work.
- The standalone draft is 95 pages. In the full prelim (as of 2026-11-05) it is 51 body pages and
  34 appendix pages, with 27 figures and 28 tables. Results have nine subsections (about 5,000
  words) and the Application has eight (about 5,300 words). The Discussion's Future
  Considerations is about 140 words. The supplementary material is about 13,000 words, including
  about 3,900 on participation, power and sensitivity in the NC application.
- Pending and approved (DECISIONS 2026-11-11, TODO "Next up"): reparameterize the outcome model
  so that β multiplies the treatment indicator Z (the intervention effect), τ multiplies the
  incidence covariate X, and spillover is γ S(Z) with S(Z) = WZ. Apply it consistently to
  Chapter 3, CTJ, the supplement, figure and table labels and code comments. Check what this
  means for Chapter 2's notation (Project 1 is published, but the prelim Chapter 2 should match)
  and ask me before touching Chapter 2.
- Also pending: state once in CTJ (around line 145) and the supplement (around lines 158 and
  197) that Checkerboard is the earlier study's Block Stratified Sampling design. This was done
  in Chapter 3 on 2026-11-11.
- FCL-050 is partly addressed: the three-decimal rule is settled, and the pooled-result /
  Table 2 question is still open.
- Earlier Chapter 3 prose reviews are in bios-dissertation/prelim/full-draft/notes/
  human-writer-review/ (ch3/, round2/ch3/, round3/ch3/; see RESOLVED.md and round2/GO-NO-GO.md).
  Do not reopen settled items.
- The written prelim is due 2026-11-15/16. Confirm the date with me.
- Priority: Chapter 3 ready for the prelim is the deliverable. The targeted revisions are made
  first in CTJ and then carried into Chapter 3, so CTJ gets its first revision round on the way.
  Do not spend effort on the CTJ 3,500-word limit now.
- Settled in the meeting (see the notes file): the β/τ/γ reparameterization, a pre-specified
  incidence coefficient whose value is given in the simulation settings, "Application"
  renamed to a "Case Study" (a simulation based on the North Carolina data), cited t and
  normal intervals, and estimating the incidence coefficient after a trial as future work.

## Principles carried over from the Chapter 4 plan

- Curate, don't dump. The body carries the main argument. Detail a reader can consult goes to
  the appendix. Hold out anything that doesn't change what the reader understands, and record
  it in paper/thesis_chapter/notes/trimmed-material.md. Never delete it. Avoid bloat in the
  body and the appendix alike, but don't drop results that are valuable.
- Page count is a guide, not a limit. The aim is focus and less white space from floats that
  stand alone on a page.
- Headings must be clear and specific to a reader scanning the table of contents, with no
  internal labels or undefined abbreviations.
- Tense: completed work in the past tense, future work in the future tense.
- The appendix is one numbered list of sections, methods sections first and then results
  sections, with no "Part" headings. The CTJ supplement's Part I/Part II layout is a reference
  for ordering only.
- Follow the standing author rules in AGENTS.md: the one design-name set, t + df intervals with
  Wald beside them, SRS in every design test, Checkerboard in the main comparison, "Chapter 2"
  allowed in the chapter but never in CTJ, no repository links in the manuscript, incidence
  guidance framed as stated there, and the three-decimal rule.
- Writing: apply the human-writer skill. No mannered prose. Avoid semicolons and dashes as
  punctuation. Define every abbreviation at first use. Title Case headings.
- "Fitted" only as an adjective before a noun ("the fitted model", "the fitted coefficients").
  The passive verb is "fit": "the model is fit by maximizing...", never "is fitted by" or "is
  fitted to". The writing review flags every misuse.

## Math accuracy

Every equation, symbol and formula in Chapter 3, CTJ and the supplement must be correct for the
current code and results. This applies to the text as it stands and to every edit.
- Each symbol is defined at or before its first use, with the same meaning across the three
  documents (after the β/τ/γ reparameterization, if approved).
- Math taken from outside work (estimators, the t interval and its degrees of freedom, the Wald
  interval, spillover models, design definitions) is cited to its source.
- Math that describes this project's own model, data-generating process, designs or estimators
  is checked against the code that produced the results: code/00_mathematical_specification.Rmd,
  code/02_incidence_generation.R, code/03_designs.R, code/04_estimation.R, the interval scripts
  (code/26_grid_tdf_intervals.R, code/27_interval_comparison.R) and the case-study code in
  application/. Check parameter values against the settings actually run, not the comments.
- The Phase 3b subagent writes one row per equation or symbol: location, the math as written,
  its source (citation, or file and line of the code), status (verified / cited / wrong /
  undefined / inconsistent / uncited), and a proposed fix. Wrong or inconsistent math is fixed
  before Phase 5, through the Phase 4 sign-off lists. If the code and the text disagree and it is
  unclear which is right, ask me. Never change the code to match the text.
- An edit that adds or changes math names its source in the sign-off list item.

## Work plan

Phase 1. Set-up. Confirm the deadline and sync state (git status in SpatialCRT and
bios-dissertation, draft PDFs dated after the last source commits). Record the baseline: commit
hashes and page counts for the chapter, CTJ, the supplement and the full prelim.

Phase 2. Plan, content map and feedback triage (STOP for me at the end).
 a. Write docs/plans/2026-11-14-chapter3-review-session-plan.md in this project, modeled on
    the Chapter 4 plan: goals, context, decisions, content map, phases, done-when. Add a row to
    docs/plans/README.md.
 b. Content map: for every Methods, Results, Application and Discussion subsection, and every
    figure and table, propose keep in the body, shorten, merge, move to the appendix, or hold
    out. Propose clear headings and the appendix order.
 c. Feedback triage: one disposition table covering the meeting items (M-, S-, A-, I-, W-, Q-)
    and FCL-001 to FCL-060, merging duplicates (an FCL annotation and a meeting item about the
    same passage become one row listing both IDs). Columns: IDs, issue (short), where it applies
    (Chapter 3, CTJ, supplement), the proposed change, and a status (to do / already addressed /
    proposed decline with reason / after the prelim). Use the register's provisional Chapter 3
    destinations only as hints, and check each against the current chapter text. Check the code
    before writing A-03 (how the blocks were built on the North Carolina geography).
 d. Decisions for me: whether the reparameterization lands before the prelim given that it touches figures
    and tables; Chapter 2 notation; which items wait until after the prelim; Q-02 (Q-01 is settled); the
    scope of Future Considerations.

Phase 3. Consistency check (subagent, during Phase 2). Compare every number, design name,
interval convention and notation symbol across Chapter 3, CTJ and the supplement, and against
the saved results. Output one row per mismatch to
bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/consistency-check.md.

Phase 3b. Math check (subagent, during Phase 2). See "Math accuracy" above. Output
bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/math-check.md.

Phase 4. Edits, in this order (settled): first the targeted revisions and clarifications in CTJ
and the supplement (the reparameterization if approved, then the accepted feedback items); then
carry the same changes into Chapter 3, together with the chapter's content-map restructuring.
 a. CTJ sign-off list (STOP). Before editing CTJ or the supplement, write every proposed edit to
    paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/proposed-edits.md, in the format
    of bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch1/proposed-edits.md:
    a "How to respond" note, a decisions log table at the top, items grouped by section or issue,
    and for each item an ID, the source IDs it answers (FCL-xxx, M-xx, S-xx, A-xx, I-xx, W-xx,
    Q-xx, consistency-check and math-check rows), the location (file and line), the current text
    and the proposed text. Cover every "to do" row of the Phase 2 triage table, whether it came
    from an annotation on the draft or from my meeting notes. Where an item needs a new analysis,
    figure or table rather than wording, describe the change and what it would take instead of
    giving text. Notation changes from the reparameterization get one item per passage. Items that
    change a claim or a number go in a final group that needs my decision. Rows triaged as already
    addressed, declined or after the prelim are listed at the end for confirmation, without text.
    I reply by item ID ("accept all A, reject C3, edit B2"). Record my decisions in the log table.
 b. Apply only the accepted items to CTJ and the supplement, rebuild both PDFs, and show me the
    result. Commit in SpatialCRT.
 c. Chapter 3 carry-over list (STOP). For each accepted CTJ item that applies to Chapter 3, write
    the matching Chapter 3 edit (current and proposed text from Thesis_Chapter.qmd, since the
    wording differs) to round5/ch3/carry-over-edits.md in the same format, together with the
    content-map restructuring for each section. I sign off by item ID as in (a).
 d. Apply the accepted Chapter 3 items section by section. STOP after each section so I can check
    it in the re-rendered draft. Commit in SpatialCRT after each approved group so the hooks
    re-render and each step can be reverted.

Phase 5. Independent reviews (two subagents in parallel, on the revised chapter, without this
session's conclusions):
 - Writing review with the human-writer skill: tics, typos, choppy sentences, ambiguous
   referents ("it", "this", "that"), undefined or uncited terms, jargon, overstated claims.
   Output round5/ch3/proposed-edits.md in the Chapter 1 round-5 format (grouped by issue, with
   ID, location, current text and proposed text; claim-changing items in a final group).
 - Committee-member review: inconsistencies, gaps, unsupported claims, missing citations,
   notation or number errors, any equation that is undefined, uncited or disagrees with the
   code, likely oral-prelim questions, and content that could be cut.
   Output round5/ch3/committee-proposals.md in the Chapter 1 format. Stage new references in
   round5/ch3/new-citations.bib after checking them on Crossref or PubMed.
 STOP for my accept/reject decisions by item ID.

Phase 6. Apply the accepted items, re-render the chapter, CTJ, the supplement and the full
prelim, check numbering, cross-references and page counts against the baseline, place the approved CTJ reading copy in revision_markup/ as the
next CTJ_Reading_Markup_V<n>.pdf, and update the feedback disposition table, ROADMAP.md, docs/TODO.md and DECISIONS.md here and
bios-dissertation/ROADMAP.md.

Done when: the content map is confirmed and followed; every FCL and meeting item has a recorded
disposition;
the reparameterization (if approved) is consistent across all documents; the consistency check
and the math check have no open rows; the CTJ and Chapter 3 sign-off lists have a decision on
every item and accepted items are applied; both review files have decisions and accepted items are applied; everything
renders with page counts recorded.

Whenever you stop, tell me which done-when items are still open. If the session runs long, split
after Phase 4: commit, then start fresh with this prompt plus "Phases 1 to 4 are done."
```
