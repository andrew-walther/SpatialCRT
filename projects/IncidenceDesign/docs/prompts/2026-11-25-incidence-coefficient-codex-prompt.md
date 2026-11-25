# Kickoff prompt: the incidence coefficient τ (for Codex)

Author: Andrew Walther. Written 2026-11-25. Paste the fenced block into a Codex session started in
`SpatialCRT/projects/IncidenceDesign/`. It runs Phase 4 of
`docs/plans/2026-11-25-project2-phased-plan.md`, alongside the CTJ and SI work in other sessions.

```text
Read AGENTS.md (standing author rules, critical invariants, file map), DECISIONS.md, and the plan
docs/plans/2026-11-25-incidence-coefficient-plan.md. Follow the plan.

Notation: in the documents β is the intervention effect on Z and τ is the incidence coefficient on
X. In the code they are `tau` (on Z) and `beta` (on X). This work is about the incidence
coefficient (code `beta`).
τ is a regression coefficient (X is on [0, 1], residual SD 1), not a correlation. Proposed values
{0, 0.25, 0.5, 1, 2}; the plan lists expectations to check (β̂ nearly unchanged across τ, τ̂'s SE
depending on the design), which you test rather than assume. Do not draft manuscript or SI text;
the main session writes the SI subsection after the findings.

Scope and rules:
- Do not edit code/01 to code/05, results/sim_data/, the CTJ manuscript, the SI, Chapter 3 or any
  file another session is working on (paper/ctj_manuscript/, paper/thesis_chapter/). Work only in a
  new script code/31_incidence_coefficient.R (plus a small helper file if needed, as
  code/29_igs_variant_designs.R does), results/incidence_coef_2026-11/, and new docs files.
- Reuse the patterns of code/26_grid_tdf_intervals.R and code/29_igs_sensitivity.R: source 05 with
  sim_define_only <- TRUE, manifest-checked checkpoints under a gitignored checkpoints/ folder,
  forked mclapply with n_cores = 10, BLAS pinned to one thread, t + df intervals.
- The run at τ = 1 must reproduce the production run results/sim_data/*_20260927_000502* exactly for
  the β̂ metrics; stop if it does not.
- Fail loud: report failed, aliased and warned fits; never drop rows silently.
- Script header and roxygen as in the existing scripts; "Author: Andrew Walther"; never name an AI
  tool in code or documents.

Steps, with stops:
1. Write a short spec as docs/plans/2026-11-25-incidence-coefficient-spec.md (τ values, how τ̂ and
   its SE are captured without editing 04, metrics, figures, runtime estimate). Stop and wait for
   Andrew's approval.
2. Smoke run (05's test profile), then the pilot. Check the reproduction at τ = 1 and that τ̂ looks
   sane (bias near 0 under the oracle fit). Stop and report.
3. Full run for the approved τ values, then summaries and figures in results/incidence_coef_2026-11/.
4. Findings: docs/findings/2026-11-25-incidence-coefficient-findings.md (method, results, what
   changes, limitations). Get an independent review from a fresh reviewer that sees the plan, the
   code and the outputs but not your conclusions; write it to
   docs/reviews/2026-11-25-incidence-coefficient-review.md and apply its corrections.
5. Add rows to the docs READMEs, the AGENTS file map, the ROADMAP progress log and docs/TODO.md.
   Commit with an imperative summary line and an itemized body; no AI attribution lines. Do not push.
Report: what ran, the reproduction check, the main findings in a few lines, and the decisions
Andrew needs to make about placement.
```
