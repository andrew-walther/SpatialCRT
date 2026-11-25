# Kickoff prompt: Incidence-Guided Saturation sensitivity rerun (for Codex)

Author: Andrew Walther. Written 2026-11-25. Paste the fenced block into a Codex session started in
`SpatialCRT/projects/IncidenceDesign/`. It runs Phase 3a (rerun) and 3b (figure) of
`docs/plans/2026-11-25-project2-phased-plan.md`, while the SI work (Phase 2) runs in another
session. Phase 3c (SI text) is not part of this prompt.

```text
Read AGENTS.md (standing author rules, critical invariants, file map), DECISIONS.md, the phased
plan (Phase 3), the spec docs/plans/2026-11-24-igs-sensitivity-spec.md and the findings
docs/findings/2026-11-24-igs-sensitivity-findings.md. Then read code/29_igs_sensitivity.R,
code/29_igs_variant_designs.R and code/30_igs_sensitivity_summary.R.

Notation: the documents write β for the intervention effect on Z and τ for the incidence
coefficient on X; the code uses `tau` (on Z) and `beta` (on X).

Author decisions for 3a (2026-11-25):
- Add spread s = 0.5 as a boundary case: 8 new variants (R = 2, 4, 10, 25; guided and unguided
  twin), using the existing rule sat_r = 0.5 + s(1 − 2(r − 1)/(R − 1)). For quadrants this gives
  100/67/33/0, so one region is fully treated and one fully control; for R = 2 it treats one half
  entirely. Its purpose is to show whether going past 90/10 helps or hurts, not to propose it.
  Expectation to test, not assume: under control-only spillover the gain from wider spreads levels
  off or reverses; under both-arms spillover it gets worse, as the wider spreads already did.
- No repeated-level variants (for example nine regions sharing the four quadrant levels).
- Grid only; do not run the North Carolina map.

Author decision for 3b: in the figure, replace the spread s on the x-axis with the highest and
lowest saturation as labels ("60/40", "70/30", "80/20", "90/10", "100/0"), keep the circles,
and explain them in a legend (they mark the original design and Saturation Quadrants). No
in-figure title or note (the caption carries them); true hyphens and minus signs; text legible at
manuscript size. PDF and PNG.

Scope and rules:
- Do not edit code/01 to code/05, results/sim_data/, the CTJ manuscript, the SI, Chapter 3, or
  anything under paper/ (another session is editing the SI there). Work only in code/29_*,
  code/30_*, results/igs_sensitivity_2026-11/, and docs files for this study.
- Keep the existing variant IDs 101–130 unchanged; new variants get 131–138 and their own seed
  keys. Existing variants and Designs 3, 8 and 9 must reproduce the 2026-11-24 run (and Designs
  3, 8, 9 the 2026-09-27 production run) exactly; stop if they do not. If only the new variants
  need fitting, reuse the stored results through the existing manifest checks rather than
  refitting everything, and say which path you took.
- Check identification for every new variant (rank of [1, Z, WZ, X] for each draw, aliased
  fits, failed fits) and report it. A fully treated or fully control region is allowed.
- Adding variants changes the Holm families. Report the adjusted p-values for the original 30
  variants before and after, so any change in what is detectable is visible.
- Forked mclapply with n_cores = 10, BLAS pinned to one thread. Fail loud: never drop rows or
  warnings silently.
- Keep any new output file under 50 MB; if one would be larger, gitignore it and say so.
- Another session commits in this checkout. Stage only your own paths by name (never
  `git add -A`, `git add .` or `git commit -a`), and leave others' uncommitted files alone.
- Script headers and roxygen as in the existing scripts; "Author: Andrew Walther"; never name
  an AI tool in code, documents or commit messages.

Steps, with stops:
1. Update the spec (status, the s = 0.5 variants, the new multiplicity families, runtime) and
   stop for Andrew's approval.
2. Smoke run of the new variants, then the full run. Check the reproduction.
3. Update the summaries and the figure (3b).
4. Update the findings report: what s = 0.5 shows against the expectation, any change in the
   Holm results for the original variants, and limits. Get an independent review from a fresh
   reviewer who sees the spec, the code and the outputs but not your conclusions; write it to
   docs/reviews/2026-11-25-igs-sensitivity-rerun-review.md and apply its corrections.
5. Update the docs READMEs, the AGENTS file map row, and the ROADMAP progress log. Commit with an
   imperative summary line and an itemized body. Do not push.
Report: what ran, the reproduction check, the s = 0.5 result in a few lines, any change in the
original findings, and the figure paths.
```
