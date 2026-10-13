# Chapter 3 length review: plan and kickoff prompt (2026-10-12)

## Status (2026-10-13): complete and adopted; see the ROADMAP progress log.

### Earlier status (2026-10-12)

Steps 1-3 done on a separate copy at the author's request
(`paper/dissertation_chapter/Dissertation_Chapter_concise.qmd`, branch `ch3-concise`).
Body 70-115 -> 70-109 in the prelim preview. Step 4 is recorded in ROADMAP; the
bios-dissertation records wait until the author adopts the copy (adoption steps in
the ROADMAP 2026-10-12 entry). Author decisions: 36-38 pp acceptable target, Ch2-or-shorter
is the guide; keep any result with real insight; keep the error-vs-tau appendix (B.6.2, now B.7.2)
unchanged; merge the annual NC tables listing every year for every design; keep the
NC maps at full size; SI unchanged except the merged table at adoption.

## Goal

Make Chapter 3 (Project 2) shorter and easier to follow without losing results the
author needs. Each subsection and appendix section should either support the
chapter's argument or move to the appendix. Material that adds little should be
condensed or cut, and cut text is kept in a side file. The chapter is the long-form
canonical version; it has no page limit, but length should be earned.

## Current inventory (full prelim, 2026-10-12 build, 194 pages)

Body: pages 70–114 (45 pages). Appendix B: pages 137–161 (25 pages).

| Section | Pages | Approx. length |
|---|---|---|
| 3.1 Introduction (3.1.1–3.1.3) | 70–73 | 3 |
| 3.2 Methods: outcome model, incidence, designs, estimation, simulation design | 73–82 | 9 |
| 3.3.1 Design performance under queen contiguity | 82–84 | 2 |
| 3.3.2 Formal statistical comparison (Friedman/Nemenyi) | 84–86 | 2 |
| 3.3.3 Design rankings by incidence configuration | 86–87 | 1 |
| 3.3.4 Sensitivity to ρ, γ and spillover regime | 87–90 | 3 |
| 3.3.5 Rook contiguity | 90–92 | 2 |
| 3.3.6 Non-oracle sensitivity | 92–93 | 1 |
| 3.3.7 Stability across treatment effect sizes | 93–95 | 2 |
| 3.3.8 Robustness across blocks and the recommended design | 95–96 | 1 |
| 3.3.9 SRS benchmark and allocation-specific risk | 96–98 | 2 |
| 3.4 Application (inputs, adaptations, model/MC, results, downside) | 98–109 | 11 |
| 3.5 Discussion (summary, relation to Ch. 2, limitations, future work) | 109–114 | 5 |
| B.1 Eight-design comparison | 137–138 | 1 |
| B.2 Simulation mechanics | 138–141 | 3 |
| B.3 Supplementary statistical comparisons | 141–144 | 3 |
| B.4 Rook contiguity figures | 144–146 | 2 |
| B.5 Oracle and non-oracle estimates | 146 | <1 |
| B.6 Technical notes (identification; error vs τ) | 146–151 | 5 |
| B.7 Observed-incidence application detail | 151–161 | 10 |

Recompute these from a fresh build before proposing changes; the numbers above
are a starting point.

## Hypotheses to test (not decisions)

- 3.3.2, 3.3.3 and 3.3.8 may repeat one ranking result three ways. They could merge
  into one subsection, with the test statistics moved to B.3.
- 3.3.5 (rook) and 3.3.6 (non-oracle) are sensitivity analyses. Each might become
  one paragraph in the body, with tables and figures in the appendix.
- 3.3.7 (τ stability) may duplicate B.6.2.
- 3.2.5 (simulation design) may overlap B.2 (simulation mechanics).
- 3.4.2 and 3.4.3 (geographic adaptations, MC evaluation) may carry procedural
  detail that belongs in B.7. B.7 itself (10 pages) needs its own check.
- 3.5.1 may restate results at length rather than interpret them.
- B.6 (5 pages) and B.7 (10 pages) should each be checked for material that is
  never referenced from the body.

## Constraints

- **Edit only** `paper/dissertation_chapter/Dissertation_Chapter.qmd` in SpatialCRT.
  Committing it runs the post-commit hook, which syncs Chapter 3 into
  bios-dissertation. Then rebuild the full prelim there
  (`prelim/full-draft/tools/build_prelim.py`) and run its tests.
- **Never delete cut material.** Move it to
  `paper/dissertation_chapter/notes/trimmed-material.md` with the date and the
  reason.
- **Supplement mirror.** `paper/tools/verify_manuscripts.py` requires every table
  label in the chapter to appear identically in `Supplementary_Information.tex`.
  Decide with the author whether the SI follows each chapter change; if a table
  changes in the chapter, change the SI copy too.
- The CTJ manuscript is separate and stays as it is unless the author asks.
- **No invented numbers.** Every number must trace to a results file. When text is
  condensed, keep the numbers that remain unchanged.
- **No mannered prose** (repo-wide rule). Keep the author's terminology:
  "fixed" regions, "Simple Random Sampling" in tables and figures, and incidence-guided
  saturation presented as a strong design, without downplaying incidence.
- Measure length in the prelim class format (the synced draft PDF or the full
  prelim), not the standalone article render.

## Steps

1. **Inventory (read-only).** For each subsection and appendix section, record:
   its pages, what it shows, which later sections cite it, any overlap with
   other sections, and a recommendation (keep / condense / move to appendix / cut),
   with an estimated page saving.
   → verify: the page totals match the fresh build.
2. **Author review.** Present the inventory as a table with the recommendations.
   Ask for a target length, or confirm that there is none. Ask which results must stay
   in the body. **Make no edits before approval.**
3. **Edit in approved batches,** one commit per batch, with cut text in the side
   file each time.
   → verify: the chapter renders without `??`; the hook sync reports OK; the
   verifier passes; the full prelim rebuilds and its tests pass; record the new page
   counts.
4. **Record** the outcome in `ROADMAP.md` (progress log), in IncidenceDesign
   `AGENTS.md`, and in bios-dissertation `AGENTS.md`/`README.md` page counts.

## Kickoff prompt (paste into a fresh session started in `~/GithubProjects/SpatialCRT`)

> Read `projects/IncidenceDesign/AGENTS.md` and
> `projects/IncidenceDesign/docs/plans/chapter3_length_review_plan_2026-10-12.md`.
> Carry out step 1 of the plan: a read-only, section-by-section inventory of
> Chapter 3 (`paper/dissertation_chapter/Dissertation_Chapter.qmd`) and its
> appendix, with page counts from a fresh prelim-format build, citations between
> sections, overlaps, and a keep / condense / move / cut recommendation with
> estimated savings. Test the plan's hypotheses rather than assuming them. Then
> stop and present the table for my review. Do not edit the chapter yet.
