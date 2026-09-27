# Handoff (2026-09-27): SRS framing, Chapter 2 check, chapter write-up, CTJ

Context for the next agent (Codex or Claude). Read first: `AGENTS.md`, `docs/plans/TODO.md`
(item 1 status), this file, `results/srs_benchmark/srs_benchmark_summary.txt`, and
`paper/dissertation_chapter/notes/a5_technical_notes.md` Q3.

**Continuation update (2026-09-27):** Task 1's restricted code check and the subsequent
accepted-PDF review are recorded in [ch2_check_findings_2026-09-27.md](ch2_check_findings_2026-09-27.md).
Read that note before acting on the original task list below. Verdict (c) is limited
to the numerical reproduction: current code/CSV agree, but disagree with the accepted
exhibit whose historical inputs remain untraced. Project 1 stays finalized/read-only;
the author has not approved manuscript changes or settled the Task 2 decisions.

## State
- SRS (Design 9: complete randomization, 50 of 100) is in the simulation. The re-run is
  verified, and the Design 1–8 rows are identical to the old run.
- `code/15_srs_benchmark.R` produces the SRS tables. Figures 12/14 draw SRS as gray
  reference lines.
- Commits: 0052161, 197aa62, 3366af2, a0772bf (all local, unpushed).
- The chapter is not yet updated for SRS.

## Key findings
Oracle estimator, RE = MSE / MSE_SRS.

**Pooled, queen, τ = 1.**
- IGSQ RE 0.71: the only proposed design that beats SRS significantly.
- Balanced Quartiles RE 1.02 (n.s.).
- Isolation Buffer 1.25, High Incidence Focus 1.91, 2x2 Blocking 2.49, Checkerboard 8.46.

**By spillover regime.** Share of scenarios beating SRS, over all τ.

| Design | queen both | queen control-only | rook both | rook control-only |
|---|---|---|---|---|
| IGSQ | 26% | 100% | 33% | 95% |
| Isolation Buffer | 0% | 100% | 0% | 50% |
| Balanced Quartiles | about 50% | about 50% | about 50% | about 50% |
| HIF, 2x2 Blocking | about 0% | about 0% | about 0% | about 0% |
| Checkerboard | 0 of 3,200 scenarios | | | |

- ρ, γ and τ barely matter. Min Checkerboard RE is 1.10.
- Non-oracle, queen, both arms: Checkerboard wins 36% (sensitivity only).
- Scratch script for these tables:
  `/private/tmp/claude-501/-Users-ajwalther-GithubProjects-SpatialCRT-projects-IncidenceDesign/7a4523d7-0438-4686-b5bc-7e0dd793a525/scratchpad/wins.R`
  (may be gone). Rebuild it inside `code/15_srs_benchmark.R` as a regime-breakdown section.

**Conflict with Chapter 2.** Project 1 is accepted at BMC Medical Research Methodology and
locked in. It reports that block-stratified sampling (checkerboard) beats SRS. Example: on the
3×4 grid, control-only, ρ = 0, ψ = 0.5, intervention-effect MSE is 0.0004 vs SRS 0.034.
The draft's own collinearity result (z = 1 − x ⇒ β̂ ≈ β − ψ) predicts MSE ≈ ψ² ≈ 0.25 for those
checkerboards. The 2×4 table is consistent with that result (bias ≈ −ψ); the 3×4 and 3×3
tables aren't.

## Task 1 — minimal Chapter 2 code check (do first; read-only on Project 1 files)
Code: `projects/SpillSpatialDepSim/code/SpatialSim_3x4.Rmd`. Also `SpatialSim_3x3.Rmd` and
`SpatialSim_NC_DOC.Rmd` (the 2×4 grid) for comparison.
- **Author's note (2026-09-27):** Chapter 2 simulated SUBJECTS within clusters. All subjects
  in a cluster share the cluster response plus an individual random error, and the model is
  fit at the subject level.
- The author believes (2026-09-27) that W was built at the SUBJECT level. Confirm this in the code first:
  whether the spatial weights / lagsarlm fit are built at the subject level. If W links
  subjects, the cluster-level collinearity z = 1 − x may not carry over to the fit.
- Check also:
  - how x and z (spillover indicator) are built per subject;
  - whether allocations 313 and 612 really are the checkerboards;
  - how MSE is aggregated (per allocation, averaged over the 2 BSS allocations; SRS = mean
    over all 924);
  - the σ used.
- Rerun ONLY allocations 313 and 612 plus ~5 random allocations with the existing code, and
  write output to a scratch dir (never overwrite `SpillSpatialDepSim/results/`). Compare the
  checkerboard's intervention-effect (β) MSE with the published 0.0004 and with ψ² ≈ 0.25.
- Report one of three verdicts, with evidence (file:line):
  - (a) Chapter 2 numbers reproduce, and the subject-level model explains why they differ
    from Chapter 3 (name the mechanism);
  - (b) they reproduce but conflict with the draft's text;
  - (c) they don't reproduce, which suggests an error.
- Don't edit Project 1 or bios-dissertation files. If the verdict is (c), stop and tell the
  author; it's a question for the advisor (Feng-Chang Lin).

## Task 2 — author decisions to settle with the author before writing
- **SRS framing:** keep it restrained, as a benchmark only.
  - Report by spillover regime, with queen primary; pooled numbers are labeled "pooled".
  - Descriptive prose, e.g. "relative to complete randomization, IGSQ and Isolation Buffer are
    more efficient under control-only spillover".
  - Never imply that designs beat SRS where they don't. Numbers stay visible (SRS table rows,
    reference lines).
- **Checkerboard (author is considering this):** drop it from the chapter body and move it to
  an appendix note on its Z–WZ collinearity and non-identification.
  - The trade-off is weaker continuity with Chapter 2.
  - Decide after Task 1. If dropped, the body has five designs: re-run 12–14 with
    retained_ids minus Design 1 (rank statistics change), and update exhibits and prose.
- **Chapter 2 link:** one light sentence, explaining the mechanism (collinearity, and the
  subject- vs cluster-level model if Task 1 confirms it). Never critique Chapter 2 or its
  tables in Chapter 3.

## Task 3 — chapter write-up (plan steps 4–5)
Plan: `~/.claude/plans/pasted-content-id-a3af-you-are-synchronous-seal.md` (steps 4–5).
- SRS goes into Methods, Results and Discussion, with benchmark rows plus MSE/MSE_SRS in the
  tables:
  - Table 2 and Table 3;
  - the ρ/γ/regime and robustness tables;
  - A1 and A5;
  - the A2 grid, which becomes 14,400 scenarios.
- SRS is not added to Figure 1 or Table 1.
- Copy the regenerated figures from `results/six_design_manuscript/`.
- Render with RStudio Quarto: `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`,
  with TinyTeX on PATH.
- Review in fresh passes until there is no ERROR/WARN. The brief is
  `docs/plans/step0-chapter-review-prompt.md`.
- Style rules:
  - no DIM, and no numbered design shorthand;
  - "we" / "this chapter" / "Chapter 2";
  - queen primary, and rook figures in the appendix;
  - `\label`/`\ref` cross-references;
  - default `pdf()` device, plotmath Greek, and the fixed viridis palette (end = 0.85).
- Commits touching the chapter trigger a post-commit sync to bios-dissertation. Don't edit the
  bios-dissertation copy.
- Then update AGENTS.md, README and TODO.md.

## Task 4 — CTJ manuscript (TODO item 3)
- Derive it from the finished chapter: body about 3,000–3,300 words, at most 6 exhibits, an
  unstructured abstract of at most 250 words (check SAGE's live guidelines), and cut material
  moved to supplementary files.
- It must stand alone:
  - no "Chapter 2" references;
  - link SRS to Project 1 only by citing the accepted BMC paper ("as in [Paper 1]").

## Rules
- Commit per sub-step with no AI co-author lines; script headers say "Author: Andrew Walther".
- Never push without asking.
- Stop and ask at real decision points.
