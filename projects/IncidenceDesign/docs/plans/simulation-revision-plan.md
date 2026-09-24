# Plan: Revise the IncidenceDesign simulation, re-run, regenerate results, then resume manuscript unification

## Handoff (done 2026-09-24, in the planning session; kept for the record)

1. Copy this plan verbatim to
   `projects/IncidenceDesign/docs/plans/simulation-revision-plan.md`.
2. In `docs/plans/manuscript-unification-plan.md`, insert a short "Step 0.5 (added
   2026-09-24)" pointer to that file, noting that step 0 is partly done (0a snapshot,
   review brief, findings) and that the bib location changed to
   `paper/SpatialCRT_IncidenceDesign.bib`.
3. In `docs/plans/step0-review-findings-2026-09-23.md`, fill in the "Decisions needed"
   section with the decisions table below.
4. Add a 2026-09-24 status note at the top of `projects/IncidenceDesign/ROADMAP.md`
   (revision planned, why, and a pointer to the plan).
5. Commit: "Plan step 0.5: simulation revision before manuscript rewrite".
6. Give the user a resume prompt for a fresh session, starting at Phase A and then B3's
   estimator-equivalence test first within Phase B.

Everything below is for the new session.

## Context

The step 0 chapter review (`docs/plans/step0-review-findings-2026-09-23.md`) found results
that rest on implementation oversights in the simulation, not on wording:

1. **Design/outcome incidence mismatch.** Designs use `X_matrix[,1]`; outcome replicate k
   and the MLE covariate use `X_matrix[,k]`. The columns are independent (cor ≈ 0).
2. **Poisson populations.** 1,000 people per cluster gives 0.35 expected deaths and ~5
   distinct values. High Incidence Focus (`> median`) treats 14–50 instead of 50, and
   `ntile()` breaks ties by grid position.
3. **Seeding.** No reproducible seeding: X and ε are drawn before the per-scenario seed, and
   forked workers reseed from time/PID.
4. **Iteration counts.** Deterministic designs are copied 25× against the same noise, and
   only 10 noise draws exist per (nb, ρ). `N_Valid_Est = 250` and `add_mc_ses()` overstate
   precision.
5. **Silent aliasing.** Checkerboard × rook gives WZ = 1 − Z exactly. `lagsarlm` drops
   `Spill` and returns τ − γ, and `suppressWarnings` hides the warning.
6. **Stale checkpoints.** All 25 April checkpoints load silently on a re-run.

The scratch pilot (36k fits) showed which results these fixes move and which they don't:
- The top designs (Incidence-Guided Saturation Quadrants, Balanced Quartiles) and the
  bottom (Checkerboard) are stable.
- High Incidence Focus degrades, because |cor(Z, X)| = 0.88.
- The framework is unchanged: SDM outcome, oracle `lagsarlm`, six designs, parameter grid,
  metrics, Friedman/Nemenyi.

The user chose to fix everything and re-run before any committee or journal review, then
continue the chapter-first manuscript plan.

## Decisions recorded (user, 2026-09-23/24)

| Topic | Decision |
|---|---|
| Narrative | Robust vs. fragile; rankings as supporting evidence; failure diagnostics as context; **one recommended design** for investigators (chosen after the re-run) |
| M1 | **Matched surfaces**: for surface k, designs are drawn from X_k, and the outcome and analysis use X_k |
| M2 | Simulation: **fixed 100,000 per cluster** (≈35 deaths/yr, Mirzaei base rate). Application: real county/CC populations when the data arrive |
| M8 | **Non-oracle MLE included** as a sensitivity analysis (`Y ~ Z + X`, same draws); oracle stays primary |
| M6 | **Keep Checkerboard-rook estimates, flagged**, so all 6 designs appear in every ranking; explain the anomaly (Chapter 2's checkerboard paradox) in Results/Discussion |
| Neighbors | Both, never pooled for headlines. **Queen primary**; rook is the Chapter 2 continuity case and sensitivity analysis |
| Pooling | Only when labeled, with splits wherever a conclusion changes. CTJ main text leads with the most application-realistic mode; full detail in the chapter and SI |
| Chapter 2 | Brief, explicit, cited continuity; answer "why keep checkerboard in Project 1, then extend it" |
| Runtime | SEs from the fit object; lean estimator only if it matches `lagsarlm` to numerical precision (then reported as a computational note); no Longleaf |
| Bib | One bib: `paper/SpatialCRT_IncidenceDesign.bib`, self-contained in this repo (overrides the plan's `prelim/references.bib` pointer) |
| Old outputs | Archive, then regenerate, including 07, 09 and 11 |
| ROADMAP | `projects/IncidenceDesign/ROADMAP.md`, plus the Project 2 row in `bios-dissertation/ROADMAP.md` |

**Out of scope (noted, not done):**
- DIM re-run.
- Poisson ρX = 0 config (offered to the user separately).
- Lean-engine port to the application.
- Deleting dead code: `complete_after_mle.R`; `longleaf_setup/simulation.R` just gets a
  "superseded by the 2026-09 revision" note in its README.

## Phase A — Step 0 items independent of the re-run (commits "Step 0: …")

1. **Bibliography** (one commit):
   - `git mv paper/IncidenceDesign_shared.bib paper/SpatialCRT_IncidenceDesign.bib`.
   - Repoint the `shared-refs.bib` symlink in the same commit; the chapter YAML is unchanged,
     since the underscore problem is known.
   - Repoint the CTJ `\bibliography{../SpatialCRT_IncidenceDesign}`.
   - Replace `baird_optimal_2018` and both `leung_*` keys with the master
     (`bios-dissertation/prelim/references.bib`) entries verbatim.
   - Add `moulton_covariatebased_2004`, `lesage_introduction_2009` and `moran_notes_1950`
     verbatim.
   - Move `spatialCRT.bib` (an unused 60-entry subset) into
     `archive/IncidenceDesign_Manuscripts_PreUnification_2026-09-23/` with a README line.
   - Point `paper/report/IncidenceSpatialCRT_Report.qmd`'s broken `../section_drafts/references.bib`
     at the new bib.
   - → verify: all chapter keys resolve; the 3 corrected entries match the master; CTJ
     compiles; the chapter renders without citation warnings (RStudio-bundled Quarto 1.10.18).
2. **Abstract.** Move the chapter abstract verbatim into
   `paper/dissertation_chapter/manuscript-declarations.md` (step 7 adds CTJ's Declarations).
3. **Results-independent review fixes:** citations #1–6; Chapter 2 facts and naming #7, #8,
   #62; voice #63 ("we" for actions, "this chapter" for the document, as in Chapter 2);
   Intro/Discussion mannered prose #49, #50, #60; Lenoir CC (service area = Lenoir, Greene
   and Jones counties; spillover to *adjacent* service areas). Methods and Results wait for
   Phase D. → verify: the diff touches only those sentences.

## Phase B — New step 0.5: simulation revision + re-run (commits "Step 0.5a…")

**B1. Spec.** Write `docs/plans/simulation-revision-spec.md`: before/after data-generating
process equations, seed keys, result schema, the estimand (performance averaged over
incidence surfaces), and the dependence structure below. It is the source for the 00 spec
and the chapter Methods.

**B2. Code** (`code/02–06`, `10`, `12–14`):

- **M1** (`05`): loop over surfaces k = 1..10.
  - Z_{k,j} = `get_designs(d, 25, N, X[,k], nb, coords)`.
  - Y_{k,j} = (I − ρW)⁻¹(τZ_{k,j} + spill_{k,j} + βX_k + ε_{k,j}).
  - Each of the 250 fits gets its own noise column.
  - Deterministic Checkerboard reuses Z across j, with distinct noise.
- **M2** (`05`): `pop_per_cluster <- 100000`. Replace the hard-coded 1000 in `06`
  (L269/272/313/316) with the parameter.
- **M3** (`03`):
  - Ties are re-broken randomly for every design resample via
    `rank(x, ties.method = "random")`, for designs 2, 6, 7 and 8.
  - High Incidence Focus treats exactly N/2; Balanced Quartiles uses equal rank-quartiles;
    Balanced Halves uses exact halves; Incidence-Guided Saturation Quadrants ranks quadrant
    means.
  - `is_design_deterministic()` becomes Checkerboard only.
  - Apply the same M3 change to `application/code/application_designs.R`, keeping its
    `lagsarlm`. The application is **not** re-run now.
- **M4** (`05`), keyed with `sprintf("%.2f")` fields and
  `set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")`:
  - X: `("X", mode, ρX, k)`, shared by every block in a config (Invariant 1).
  - Z: `("Z", mode, ρX, nb, ρ, γ, spill, k, d)`.
  - ε: `("eps", mode, ρX, nb, ρ, γ, spill)`.
  - This gives common random numbers across designs, τ and estimator, and independent Z and
    ε across Friedman blocks. Blocks still share X surfaces; the spec states this, and B3
    adds a surface-level robustness check.
- **M6/M7** (`04`):
  - `withCallingHandlers` records warnings. `N_Aliased` counts the aliasing warning and
    `N_Warn` counts everything else, with messages logged.
  - New design-level flag `Z_WZ_rank_deficient = rank([1, Z, WZ]) < 3`, recorded in **both**
    the oracle and non-oracle files. The non-oracle Checkerboard-rook bias has the same
    cause but raises no warning.
  - Estimates are kept (user decision).
- **M8:** fit the non-oracle model on the same data. Save it as
  `sim_results_MLEnonoracle_tau_sweep_*.rds`, which `load_latest_results("MLE_tau_sweep")`
  won't match.
- **Monte Carlo SEs:**
  - 05 also saves a long per-(scenario, k) summary table (~15k × 10 rows per estimator) and
    stores `SE_Bias`, `SE_MSE`, `SE_Coverage` and `SE_Power` computed from the 10
    surface-level means, with t₉ intervals.
  - `add_mc_ses()` returns the stored columns when they're present, which protects
    `plot_mse_vs_tau`, `plot_coverage_vs_tau` and the Report.
  - Design-difference SEs are paired at the surface level.
- **Lean estimator** `fit_sar_lag()` (`04`), replicating `lagsarlm(method = "eigen")`:
  - `eig = Re(spatialreg::eigenw(listw))` computed once per nb.
  - Interval `1/range(eig) ± .Machine$double.eps`: rook ≈ (−1, 1); queen lower bound ≈ −1.97.
  - `optimize(tol = .Machine$double.eps^0.5)` on the concentrated log-likelihood, with SSE
    from the e.a/e.b/e.c decomposition and s² = SSE/n.
  - Full (σ², ρ, β) analytic information matrix, inverted with
    `solve(tol = .Machine$double.eps)`, varb = s²·inv.
  - Replicate the warnings, or return NA plus a flag, for rho-on-bound and
    inversion-failure cases.
  - LINPACK QR (tol 1e-7) drops the same aliased column; SEs are indexed by name.
  - Intercept included; `engine = "lagsarlm"` retained.
  - If validation fails, use `lagsarlm` with `fit$rest.se` on 10 cores (~7 h).
  - If it passes, report it as a computational note in the chapter Methods/appendix and the
    CTJ SI: the same estimator via concentrated likelihood, the tolerances, and the speedup.
    Validation output goes to `results/estimator_validation/`.
- **Runner** (`05`):
  - `mclapply(mc.preschedule = FALSE)` over (config × nb × ρ) units; BLAS threads set to 1
    in workers.
  - A `NULL` or `try-error` result **stops** the run.
  - Checkpoints go to `results/checkpoints/rev_2026-09/` with a manifest (parameter hash,
    code-file hashes, package versions); a mismatch refuses to load.
  - Old checkpoints move to `results/checkpoints_pre_revision_2026-09/`, added to
    `.gitignore`.
- **Downstream code:**
  - 12, 13 and 14 split or filter by `Neighbor_Type` (queen primary, rook separate) in
    their tests, figures, `best_to_worst_6` and extract sections. Pooling appears only where
    labeled.
  - 12 takes an estimator tag and writes non-oracle output to its own directory, so it
    doesn't overwrite `six_design_comparison_report.rds`.
  - NA guards in 13 (L71) and 14 (L340/448); drop 14's hard-coded Poisson 0.20 label.
  - Remove hard-coded counts: 12 L62–65, 13 L66, 06 L257–322 (5-config overview), 08
    validators (1,920 rows).
  - Rook exhibits carry a "τ not identified for Checkerboard" note, and pooled summaries
    label or exclude those rows.
  - 06 and 08 are always called with `estimation_mode = "MLE_tau_sweep"`; the 08 defaults and
    `validate_*` are updated from `"MLE_combined"`.

**B3. Tests** (`code/tests/test_simulation_revision.R`, base-R `stopifnot`, run with Rscript):
- High Incidence Focus gives exactly 50 treated on tied X.
- Quartile/half strata are balanced, and tie assignment is uncorrelated with grid index over
  1,000 draws.
- Design k is drawn from X[,k].
- The same unit gives identical results sequentially and in parallel, and in any order.
- Oracle: aliased in 100% of Checkerboard × rook fits, 0% elsewhere (γ ≥ 0.5).
  `Z_WZ_rank_deficient` is TRUE only for Checkerboard × rook, in both files.
- **Estimator equivalence:** ≥5,000 fits over all designs × nb × ρ × spill × modes ×
  oracle/non-oracle, including queen Checkerboard and aliased fits. Pass: |Δτ̂| < 1e-6,
  |ΔSE|/SE < 1e-5, |Δρ̂| < 1e-6, |ΔlogLik| < 1e-6, and identical aliasing.
- **Known answer:** Balanced Quartiles with γ = 0.5 (oracle), plus γ = 0 with the non-oracle
  model, gives bias and coverage within 3 MC SE of 0 and 0.95.
- **Surface-level robustness:** for the key pairwise claims, a paired test over 50
  independent (config × k) blocks.

**B4. Pilot** (τ = 1; all configs; ρ ∈ {0, 0.5}; γ ∈ {0.5, 0.8}), then **stop and report**:
- Directions match the scratch pilot.
- Poisson X has a median of ≥ ~55 distinct values.
- `Mean_Treated` = 50 for High Incidence Focus.
- Aliasing only for Checkerboard × rook.
- The **8-design consolidation check** (14(a)) still shows Saturation Quadrants and
  Balanced Halves redundant. If not, stop and decide the design set with the user.
- Record the projected runtime.

**B5. Archive + full run.**
- `results/archive/pre_revision_20260924/` with a README (pre_tau_sweep convention). Move
  the superseded rds files (they're untracked, so note that in the README),
  `six_design_manuscript/`, `eight_design_supplementary/`, `mle_per_config/`, the rendered
  00/07/09/11 reports, `MLE_*` PDFs, `results/figures/design_samples_*`, and the stale
  `07_results_summary.qmd` / `09_…Rmd` sources.
- Full run: 5 configs × 2 nb × 4 ρ × 4 γ × 2 spill × 8 designs × 5 τ = 12,800 scenarios ×
  250 fits × 2 estimators. Save `sessionInfo()` and the manifest.
- 1% cross-check with `lagsarlm`.
- → verify: 12,800 rows per estimator; `N_Valid_Est` = 250; aliasing and rank-deficiency
  flags only where expected; every non-aliasing warning listed in the log.

## Phase C — Regenerate downstream (commit "Step 0.5f: …")

**Run order:** 05 → 06 → 08 → 10/11 → 12 (oracle, then non-oracle) → 13 (extended) → 14,
with 14(d), the application table, labeled stale until the application re-run → regenerate
design-sample figures → copy manuscript figures into `paper/*/figures/` → render 00, 07, 09
and 11.

**13 additions:** rook/queen splits; per-mode win counts; `Mean_Treated`; flag summary;
non-oracle headline; average-rank order per τ; the surface-level robustness result.
Note that per-τ Friedman tests share draws (common random numbers) and aren't independent
confirmations.

**Text updates:**
- `00_mathematical_specification.Rmd` (full method text).
- `11_statistical_comparisons_report.qmd` (hard-coded counts, source filename line 568).
- `paper/report/IncidenceSpatialCRT_Report.qmd` and `IncidenceDesign_ProjectSummary.qmd`:
  re-render, then fix hand-typed numbers and the specific sentences about the changed
  mechanics; no broader rewrite.
- `README.md`; project `CLAUDE.md`: Critical Invariants **2** (pooling rule vs. the decision),
  **3** (deterministic designs), **4** (seeding), **5** (`X[,1]` reversed by M1), Current
  State, bug table; root `CLAUDE.md` headline numbers; `application/README.md`.
- CTJ/SI stay untouched until step 5; a "numbers superseded" note goes in the README,
  CLAUDE.md and ROADMAP.

→ verify: all regenerated outputs are newer than the run; a grep for old headline numbers
(0.079, 0.802, 800.9) across live docs returns only intended historical mentions.

## Phase D — Resume the manuscript plan (chapter-first, as originally specified)

- **Step 0, remainder:**
  - Rewrite Methods from the spec; rewrite Results as robust vs. fragile (queen primary;
    rook sensitivity including the Checkerboard-rook anomaly; per-mode results; non-oracle
    sensitivity; one recommended design).
  - Brief Chapter 2 continuity passage.
  - All remaining findings; replace the `knitr::kable` chunk with a static table.
  - Update `step0-chapter-review-prompt.md` (new ground truth; the spec as method authority)
    and loop a fresh reviewer until clean.
- **Steps 1, 2, 3, 4, 7** as written.
  - Step 2 adds the one-college-per-county rule (Northampton → Halifax CC, Bertie → Martin
    CC, Roanoke-Chowan CC → Hertford only; already in `cc_mapping_data.R`).
  - Step 3 adds the exhibit-fate table.
  - Step 7 checks that bib entries agree with the master before porting.
- **Steps 5–6** (CTJ at ~3,000–3,300 words; retire the old CTJ/SI) remain post-deadline.

## Documentation / checkpointing

- **`docs/plans/manuscript-unification-plan.md`:** insert step 0.5 and the decisions table;
  record the bib override; mark progress.
- **`projects/IncidenceDesign/ROADMAP.md`:**
  - New 2026-09-24 status block.
  - New items: the 2026-09 revision (M1–M8), downstream regeneration, and
    **future work** (latent-risk/noisy-snapshot data-generating process; heterogeneous
    populations in the simulation; application lean-engine port).
  - Mark Completed lines 148/154/157 as superseded.
- **`bios-dissertation/ROADMAP.md`:** Project 2 status row.
- **Memory:** a project memory with the revision decisions and why.
- Commit per sub-step; offer a push at the end of each phase; never push unasked.

## Timeline

About 1–2 working days in total. Most of this fixes or updates existing material, and with
the lean estimator the full run takes minutes.

- **Day 1:** A (≈1–2 h) → B1–B4 (spec, code, tests, pilot; ≈half a day) → **pilot stop /
  your review** → B5 full run → C (regenerate outputs and update text).
- **Day 2:** D, which is the chapter Methods/Results rewrite, the reviewer loop, then steps
  1–4 and 7.

The items that set the pace are the ones waiting on people or checks, not compute: your
review at the pilot stop, the fresh-reviewer loop until it's clean, and the step 4 render.
If the lean estimator fails validation, the `lagsarlm` fallback adds one overnight run (~7 h).

## Verification

1. `Rscript code/tests/test_simulation_revision.R` exits 0.
2. Pilot checks pass, including consolidation (B4).
3. Full-run integrity checks (B5).
4. 12/13/14 regenerate; the updated `validate_*` pass.
5. Chapter renders with no missing keys or `??`; CTJ compiles.
6. A fresh step-0 reviewer returns no unresolved findings.
