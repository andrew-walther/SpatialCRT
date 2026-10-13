# Step 3: exhibit fate and body vs. appendix

**DRAFT: for author review. Nothing here has been applied to the chapter.**

Written 2026-09-24 against `Dissertation_Chapter.qmd` after the step-2 edits (working tree on top
of commit 3205c6c). The render is 40 pp.: the body runs pp. 1–35 and the references pp. 36–40.
The pre-step-2 render was 39 pp., with the body on pp. 1–34 and the references on pp. 35–39.

The test comes from plan step 3, biased toward the body per the 2026-09-23 advisor feedback. An
item goes to the appendix only if it exists to spare the reader mechanical detail. If it holds a
result the reader needs to trust or understand the central claims, it stays in the body. When in
doubt, it goes to the body.

Checksums are the first 12 hex characters of SHA-256.

---

## 1. Inventory

### 1a. The chapter's current exhibits

| Exhibit | Content | Source | Checksum status |
|---|---|---|---|
| Table 1 | The six designs: rule, clusters treated, deterministic | Inline LaTeX (settled text) | n/a (static) |
| Table 2 | MSE (SE), bias, coverage, power, average rank, clusters treated, τ = 1, queen + rook | Inline LaTeX, from `phase_d_1e_exhibits.md` §1 | n/a (static; values from `six_design_summary.txt` / extract) |
| Table 3 | MSE by ρ, queen | Inline LaTeX | n/a (static) |
| Table 4 | MSE by γ, queen + rook | Inline LaTeX | n/a (static) |
| Table 5 | MSE by spillover regime, queen + rook | Inline LaTeX | n/a (static) |
| Table 6 | Checkerboard coverage by ρ, queen + rook | Inline LaTeX | Values match `dissertation_results_extract.txt` §8 (queen 0.9408/0.9430/0.9400/0.9391; rook 0.1607/0.1564/0.1482/0.1338), checked in step 1 |
| Table 7 | Block-level MSE quantiles, queen | Inline LaTeX | n/a (static) |
| Figure 1a/b | MSE by design (ranked), queen/rook | `figures/fig_mse_by_design_ranked_6design_{queen,rook}.pdf` | **Match** `results/six_design_manuscript/{queen,rook}/fig_mse_by_design_6design.pdf` (75c03c8dda10 / 7a5ef9ac8b33) |
| Figure 2a/b | Coverage by design, queen/rook | `figures/fig_coverage_by_design_6design_{queen,rook}.pdf` | **Match** `results/six_design_manuscript/fig_coverage_by_design_6design_{queen,rook}.pdf` (4f469d6c95b5 / 85c23acdf34d) |
| Figure 3a/b | Bias–variance decomposition, Poisson ρ_X = 0.20 | `figures/fig_biasvar_6design_{queen,rook}.pdf` | **Match** `results/six_design_manuscript/{queen,rook}/fig_biasvar_6design.pdf` (6874edb91c65 / b71ea72292dc) |
| Figure 4a/b | MSE across τ, queen/rook | `figures/fig_tau_sensitivity_6design_{queen,rook}.pdf` | **Match** `results/six_design_manuscript/fig_tau_sensitivity_6design_{queen,rook}.pdf` (fe73c44d97ed / 930b57525607) |

**Orphans in `paper/dissertation_chapter/figures/` that the chapter doesn't reference.** They are
listed here so nobody reuses them by mistake:
- `fig_mse_by_design_6design_{queen,rook}.pdf` (ed2763f025a1 / 1cb32b2d06dc): these are the
  alphabetical versions. They **differ** from the current canonical
  `results/six_design_manuscript/fig_mse_by_design_6design_{queen,rook}.pdf` (279647e50748 /
  16964b2d5092), so they are stale.
- `fig_mse_by_design_6design.pdf`, `fig_coverage_by_design_6design.pdf` and
  `fig_tau_sensitivity_6design.pdf`: dated 2026-07-02, pre-revision, pooled over neighbor types.
  Stale.

The plan's note that "the bias-variance figure ... is not discussed anywhere in the current chapter
text" no longer applies. It is Figure 3, and the Results text cites it.

### 1b. `results/six_design_manuscript/{queen,rook}/si_figures/` (regenerated 2026-09-24 by `code/14_manuscript_supplement_figures.R`)

| File | Content (per the 14 script and the SI captions) | Queen / rook checksum | Page size (pt) |
|---|---|---|---|
| `si_fig_cd_diagram_rank.pdf` | Nemenyi critical-difference diagram, average Friedman rank | 4ef06f8b7e53 / 541b25cc4604 | 720×374 |
| `si_fig_cd_diagram_mse.pdf` | The same comparison, positioned by MSE | bc0d13305d9e / bbb6ac77dbda | 792×468 |
| `si_fig_mse_heatmap.pdf` | Mean MSE by design × ρ, Poisson ρ_X = 0.20 only | 1433fb1ccb62 / 5c0a9bb58f9c | 576×324 |
| `si_fig_rank_heatmap.pdf` | Average rank by design × ρ, Poisson ρ_X = 0.20 only | 492569bfc867 / f4fd5b527f85 | 576×324 |
| `si_fig_performance_pvalue_twopanel.pdf` | Mean MSE ± SE (left); pairwise Nemenyi p-value heatmap (right) | 02928f6c6561 / f0ea7854a9b3 | 792×324 |
| `si_fig_tau_mse.pdf` | Mean MSE across τ | 3ef3291afb4c / b1cc5c74bc18 | 576×360 |
| `si_fig_tau_coverage.pdf` | Coverage across τ | 1aca89df711d / 5e97daf14bcf | 576×360 |

Also in `results/six_design_manuscript/`, and not yet in any manuscript:
`{queen,rook}/fig_coverage_tau_6design.pdf` (b086a3f49173 / 2ad44f5fa914; the two-panel
coverage + τ figure made for the CTJ) and `MLE_statistical_comparisons_6design_{queen,rook}.pdf`
(a 9-page internal report). Elsewhere in `results/`: `figures/design_samples_8panel.{pdf,png}`
(311cdc5eeaa6 / f27fbfc10d12; one sample assignment per design, all 8 designs, 2026-09-24) and
`six_design_manuscript_nonoracle/` (the non-oracle figure set).

### 1c. `results/eight_design_supplementary/`

| File | Checksum | Content |
|---|---|---|
| `eight_design_summary.txt` | 772470508a5c | Friedman/Nemenyi/Wilcoxon over all 8 designs at τ = 1 for queen, rook and a labeled pooled slice: mean MSE, coverage, average rank, the two consolidation pairs |
| `eight_design_comparison_report.rds` | a4df359d1ff4 | The same, as an R object |

No figure exists. Chapter Methods (Treatment assignment designs) already quotes the consolidation
numbers from this file.

### 1d. The CTJ Supplementary Information (`paper/ctj_manuscript/Supplementary_Information.tex`)

**The whole SI predates the revision. Every number and figure in it is superseded.** Its seven SI
figures (`paper/ctj_manuscript/figures/si_figures/*`, 2026-07-03) **all differ** from the
regenerated `results/six_design_manuscript/{queen,rook}/si_figures/*`, and each is a single
pre-revision version. Its prose in S1–S4 describes the old pipeline. Nothing in it can be copied
as is. Only its structure can be reused.

| SI exhibit | Content | Regenerated equivalent |
|---|---|---|
| Table S1 | Full parameter grid (8 designs) | Chapter Methods prose (Simulation design) |
| Table S2 | Resample counts per scenario | Chapter Methods prose (K = 10 × J = 25 = 250 fits) |
| Table S3 | All eight designs | Table 1 plus the two excluded designs' rules (Methods prose) |
| Table S4 | All eight designs ranked by mean MSE, τ = 1 | `eight_design_summary.txt` (queen + rook) |
| Table S5 | Mean MSE by design within each incidence configuration | `dissertation_results_extract.txt` (per-config ranges are quoted in Results prose; no table in the chapter) |
| Table S6 | MSE by ρ | Chapter Table 3 |
| Table S7 | MSE by γ | Chapter Table 4 |
| Table S8 | MSE by spillover regime | Chapter Table 5 |
| Table S9 | Friedman χ² by τ | Chapter Results prose (queen χ² at five τ) |
| Table S10 | MSE quantiles across blocks | Chapter Table 7 |
| Table S11 | Win rate | Chapter Results prose (95/58/7 queen; 97/60 rook) |
| Table S12 | Application-scale six-design table (synthetic) | **Deleted in step 2**; `application_table_6design.txt` is labeled STALE |
| Figures S1–S2 | CD diagrams | 1b `si_fig_cd_diagram_{rank,mse}` |
| Figures S3–S4 | MSE / rank heatmaps | 1b `si_fig_{mse,rank}_heatmap` |
| Figure S5 | MSE ± SE + p-value two-panel | 1b `si_fig_performance_pvalue_twopanel` |
| Figures S6–S7 | τ MSE / τ coverage | 1b `si_fig_tau_{mse,coverage}` |
| Figures S8–S10 | Application maps | 1e |

### 1e. Application figures

| File | Checksum | Size | Data status |
|---|---|---|---|
| `application/report/figures/community_college_service_area_clusters.png` | 3621cfe04052 (identical copy in `paper/ctj_manuscript/figures/application_maps/`) | 3600×2400 px | **Real geography.** It shows the county → CC service-area mapping and doesn't depend on SUD data |
| `application/results/full/figures/synthetic_incidence_map.png` | b9183783b167 (identical CTJ copy) | 2400×1650 px | **Synthetic placeholder incidence** |
| `application/results/full/figures/kmeans_regions_map.png` | 089f9a6f06b6 (identical CTJ copy) | 2400×1650 px | The partition is chosen by coordinates plus population, cluster and county balance, not incidence (`make_kmeans_regions()` in `application_designs.R`). But the "full" run used 2024 OSBM populations, and the real run will use the data file's populations, so the map will change |

### 1f. SI prose sections S1–S4

| Section | Content in the SI | Status against the revised simulation | Where the revised version already lives |
|---|---|---|---|
| S1 Reproducibility | Script table (lists DIM, `lagsarlm`); seeding by `digest2int(paste(...))` per scenario; X generated once and only column 1 given to designs | **Wrong post-revision** (key-based `set_seed_key()` for X, Z, eps; matched surfaces M1; `fit_sar_lag()` engine) | Chapter Methods, "Simulation design, resampling, and metrics" (the seeding paragraph) and "Estimation" (engine validation) |
| S2 Parameter grid | Tables S1–S2 | Grid unchanged; resample counts changed (250 = 10 × 25, one noise column per fit) | Chapter Methods prose |
| S3 Metric formulas | MSE = Bias² + SD²; naive MC SEs (SD/√N, delta-method SE_MSE); Fail_Rate | **Wrong post-revision** (MSE = mean(τ̂ − τ)²; MC SEs from 10 surface means; t₉ intervals) | Chapter Methods metrics paragraphs |
| S4 Estimation model | Oracle `lagsarlm` equation; "non-oracle variant (not run)" | **Wrong post-revision** (the non-oracle model has now been run) | Chapter Methods "Estimation" and Results "Non-oracle sensitivity analysis" |

---

## 2. Body vs. appendix: exhibit-fate table

Fates: **Body**, **Appendix**, **Omit**, **Later-with-real-data**.

| Item | Current location | Proposed fate | Reason (one line) | Format note |
|---|---|---|---|---|
| Table 1 (six designs) | Chapter body | Body | It defines the objects being compared, so the argument can't be followed without it | Single-column; too wide for the CTJ as is (four p-columns); the CTJ needs a two-column condensed version |
| Table 2 (headline performance) | Chapter body | Body | It is the central result | Seven columns, two panels; in the CTJ, fits only as a full-width table |
| Table 3 (MSE by ρ) | Chapter body | Body | Per advisor feedback the ρ breakdown is a result, not mechanics | Small; fine in both |
| Table 4 (MSE by γ) | Chapter body | Body | Same; it also carries the rook-Checkerboard γ trend the text explains | Two panels; CTJ appendix only |
| Table 5 (MSE by regime) | Chapter body | Body | The regime is the one factor that changes the ranking, which is central to the recommendation | Small; candidate for the CTJ main text |
| Table 6 (Checkerboard coverage by ρ) | Chapter body | Body | It supports the identification argument; step 1 confirmed it | Two rows; could fold into the text for the CTJ |
| Table 7 (MSE quantiles) | Chapter body | Body | It supports the "lowest worst case" claim in the recommendation | Small |
| Figure 1 (MSE by design, ranked) | Chapter body | Body | Visual form of the headline, showing the spread across incidence configurations | Two stacked full-width panels; the CTJ needs one panel (queen) at column width, with rook in the SI |
| Figure 2 (coverage by design) | Chapter body | Body | Shows coverage validity, including the rook-Checkerboard failure | The CTJ already merges coverage + τ into `fig_coverage_tau_6design`; the chapter keeps separate figures |
| Figure 3 (bias–variance) | Chapter body | Body | Explains *why* Checkerboard fails (variance under queen, bias under rook), which the argument needs | One configuration only (Poisson ρ_X = 0.20); the caption says so |
| Figure 4 (MSE across τ) | Chapter body | Body | Supports the τ-stability claim | Two stacked panels |
| **Per-configuration MSE table** (SI Table S5 analogue) | Not in chapter; ranges in Results prose | **Body** (new compact table, 6 designs × 5 configurations, queen; rook either as a second panel or in the appendix) | Per advisor feedback, per-incidence-mode breakdowns belong in the body; the prose now gives only ranges | Needs a writer pass: the values must come from the extract, and this would be Table 3, renumbering the rest |
| **8-design comparison, compact** | Chapter Methods prose (Treatment assignment designs, the consolidation paragraph) | **Body, as it stands** | It is already written as a robustness check *of* the six-design choice (the two consolidation pairs, with MSE, Nemenyi and Wilcoxon p-values and average ranks), which is what the revised test asks for | None |
| **8-design comparison, full table** (mean MSE, coverage, average rank, all 8, queen + rook) | Not in chapter; the chapter promises it "in the appendix of this chapter" (Introduction and Methods) | **Appendix** (A1) | Readers who want to check the consolidation can, but the argument doesn't need all 8 rows; the chapter already promises this appendix | Values from `eight_design_summary.txt`; two panels, 8 rows each; about half a page |
| Rules of the two excluded designs (Saturation Quadrants, Balanced Halves) | Chapter Methods prose (one sentence each) | Appendix (A1, as two extra rows in a Table-1-style table) | Mechanical definitions of designs that aren't compared in the body | Copy the wording from the Methods sentences and `03_designs.R`; no new claims |
| CD diagrams (`si_fig_cd_diagram_rank`, queen + rook) | Not in chapter | Appendix (A3) | The Nemenyi results are fully stated in the body prose (critical difference 0.596, 11 of 15 pairs); the diagram only visualizes them | 720×374 pt, legible at full width; use the regenerated 1b files, never the CTJ SI copies; the known label-collision issue concerns `plot_cd_diagram()` in `10_...R`, but these come from `14_...R`, so check legibility |
| `si_fig_cd_diagram_mse` | Not in chapter | Omit | Duplicates the rank CD diagram with a different x-axis | n/a |
| `si_fig_performance_pvalue_twopanel` | Not in chapter | Appendix (A3), p-value panel only if it can be split; else Omit | Pairwise p-values for all 15 pairs are mechanics behind the body's summary counts | 792×324 pt; wide |
| `si_fig_mse_heatmap`, `si_fig_rank_heatmap` | Not in chapter | Omit | One incidence configuration only, and Tables 3 and 5 already give the pooled ρ × design values | n/a |
| `si_fig_tau_mse` | Not in chapter | Omit | Duplicates Figure 4 | n/a |
| `si_fig_tau_coverage` (queen + rook) | Not in chapter | **Body** (as panels next to Figure 4) or Appendix | The body claims coverage ≈ 0.94 for every design at τ = 1 but shows coverage across τ nowhere; if the claim extends to all τ, the reader needs this | Needs a writer check that the body text makes (or should make) a claim about coverage across τ; 576×360 pt |
| `fig_coverage_tau_6design` (CTJ two-panel) | Not in chapter | Omit (chapter) / CTJ main text | Built for the CTJ's exhibit cap; the chapter has room for separate figures | CTJ narrow-column version already exists |
| Friedman χ² by τ (SI Table S9) | Chapter Results prose (queen) | Body, as prose (no change) | Already stated; a table would add nothing | The rook χ² values by τ are not in the chapter; add to the appendix only if the author wants them |
| Win rate (SI Table S11) | Chapter Results prose | Body, as prose (no change) | Already stated | n/a |
| `design_samples_8panel` (one sample assignment per design) | Not in any manuscript | **Body** (next to Table 1), six-design version | A picture of each design makes Table 1's rules immediately readable, which serves understanding, not mechanics | Regenerate as a 6-panel version, or use the 8-panel figure in the appendix; `DESIGN_FULL_NAMES` labels design 1 "Block Stratified Sampling (Checkerboard)", which the chapter now explains; needs a caption from the writer |
| Non-oracle figure set (`six_design_manuscript_nonoracle/`) | Not in chapter | Appendix (A4, MSE-by-design only), or Omit | The body's non-oracle section is prose-only and self-contained; a figure is verification, not argument | Optional; lowest priority |
| `community_college_service_area_clusters.png` | CTJ SI only | **Body, Application section** ("The planned pilot design"), full width | Real geography, independent of the data; shows the 58 treatment units and the one-college-per-county partition the text describes | 3600×2400 px with a dense 58-entry legend; fine at `\linewidth` in the chapter; for the CTJ it needs a simplified legend or a cross-reference |
| `synthetic_incidence_map.png` | CTJ SI only | **Later-with-real-data** (replace with the real cluster-level incidence map) | Built from placeholder data that step 2 removed | Don't import |
| `kmeans_regions_map.png` | CTJ SI only | **Later-with-real-data** (regenerate with the real populations) | Illustrates the region-based adaptation of Incidence-Guided Saturation Quadrants described in the Planned analysis, but the partition will change with the real population data | Once regenerated, it belongs in the Application body next to the Planned-analysis paragraph |
| Real county-level incidence map | Doesn't exist yet | Later-with-real-data (Body) | New exhibit (plan step 2; stub TODO item 2) | n/a |
| Real cluster-level incidence map | Doesn't exist yet | Later-with-real-data (Body) | Replaces the synthetic map (stub TODO item 3) | n/a |
| Real-data design comparison | Doesn't exist yet | Later-with-real-data (Body; run through this test again when it exists) | Stub TODO item 4 | n/a |
| SI Table S12 / `application_table_6design.txt` | Deleted from the chapter in step 2; CTJ SI still has it | Omit (chapter) | Synthetic and STALE | CTJ keeps it until step 5 |
| SI prose S1: reproducibility and seeding | CTJ SI (stale) | **Appendix** (A2), rewritten from the spec, **not copied** | Mechanics: script pipeline, exact seed-key tuples, RNG kinds, checkpoint manifest; the body's seeding paragraph already carries the part a reader needs | Needs the writer: content from `docs/plans/simulation-revision-spec.md` and `05_run_simulation.R`; the SI text is wrong post-revision |
| SI prose S2: parameter grid | CTJ SI (stale) | Appendix (A2), one small table | The body gives the grid in prose; a table is a convenience | Values unchanged; resample counts must use 10 × 25 |
| SI prose S3: metric formulas | CTJ SI (stale) | Appendix (A2), displayed formulas | Exact formulas (including surface-level MC SE and the paired-difference SE) are mechanics; the body's prose definitions suffice for the argument | Needs the writer: the SI's formulas (MSE = Bias² + SD², naive SEs) are wrong post-revision |
| SI prose S4: estimation model | CTJ SI (stale) | Body (already there); Appendix only for the log-likelihood and information-matrix details of `fit_sar_lag()`, if wanted | The body's Estimation section already states both models and the validation | The SI's "non-oracle not run" is wrong post-revision; don't copy |

**The explicit 8-design judgment.** Keep the compact version in the body and put the full table in
the appendix. The compact version is the consolidation paragraph already in Methods. It is written
as a robustness check of the six-design choice and gives the two pairs' MSEs, test p-values and
average ranks under both neighbor types. The full 8-row table adds only the six retained designs'
numbers again (already in Table 2) plus the two excluded rows. A reader verifying the
consolidation needs it, but the argument does not. The chapter already promises "the appendix of
this chapter" for it, twice (Introduction and Methods). So the appendix must exist and must hold
this table, whatever else is decided.

---

## 3. Proposed appendix outline and page targets

**Appendix to Chapter 3** (working title; its numbering depends on the dissertation template)

- **A1. The eight-design comparison.** A Table-1-style definitions table for Saturation Quadrants
  and Balanced Halves (2 rows), and the full 8-design mean MSE / coverage / average rank table,
  queen and rook panels, at τ = 1. Source: `eight_design_summary.txt`. About 1–1.5 pp.
- **A2. Simulation mechanics.** The parameter-grid table; the seed-key tuples and RNG settings;
  the checkpoint manifest; the script pipeline (01–05, 12–14); exact metric and Monte Carlo SE
  formulas. Rewritten from the spec and the code, not from the SI. About 2 pp.
- **A3. Supplementary statistical comparisons.** Rank CD diagrams (queen, rook), and, optionally,
  the pairwise Nemenyi p-value heatmap. Regenerated 1b files. About 1–2 pp.
- **A4 (optional). Non-oracle MSE by design** (queen, rook). About 1 p.
- **A5. Technical notes on anticipated questions** (added 2026-09-25 at author request: "good items
  to cover in the appendix if applicable"). Each note is included only if an investigation
  (computation on the results `.rds`, algebra, or the Chapter 2 draft) gives a supportable answer.
  Otherwise the question goes only to the Q&A companion doc (below). Candidates:
  1. Why the oracle bias of Incidence-Guided Saturation Quadrants grows with τ (queen 0.029 → 0.047
     at τ = 0.8 → 3), and why its control-only MSE rises while Isolation Buffer's is flat.
  2. Why rook Checkerboard coverage rises with τ (0.138 at τ = 0.8 → 0.442 at τ = 3), even though
     its bias is −γ at every τ (the SE presumably grows with τ; verify).
  3. Chapter 2's 3×4 control-only block-stratified MSE ≈ 0.0004 at ψ = 0.5 (draft L972–976) vs.
     its own β̂ ≈ β − ψ collinearity result (L1350–1356). Note: this is a Chapter 2 issue; if it's
     real, the fix belongs in Chapter 2, with Chapter 3 at most cross-referencing it.
  4. The rook Checkerboard regime gap (0.477 vs 0.483; paired p = 0.009 over 80 blocks), which is
     identical by construction, so Monte Carlo error only.
  About 1–2 pp. if all four hold up. Needs phase-d-writer (investigation plus new claims).

**Page estimate** (current render: double-spaced `article`, 1-inch margins; body pp. 1–35):

| Change | Body Δ (pp.) | Appendix (pp.) |
|---|---|---|
| Step 2 already applied (placeholder run deleted; planned analysis + stub added; pages shifted) | +1 (34 → 35) | n/a |
| CC service-area map in the Application | +1 | n/a |
| Per-configuration MSE table | +0.5 to +1 | n/a |
| Design-sample figure (six panels) next to Table 1 | +1 | n/a |
| τ-coverage panels, if placed in the body | +1 | (or +1 in the appendix) |
| A1 | n/a | 1–1.5 |
| A2 | n/a | about 2 |
| A3 | n/a | 1–2 |
| A4 (optional) | n/a | about 1 |

**Rough targets:** body **about 37–39 pp.** (from 35), appendix **about 4–6 pp.** (5.5 with A4).
The appendix lands near Project 1's 4-page scale, as the plan expected once most SI content stays
in the body. The later real-data backfill will add about 2–4 body pages (two maps plus the
comparison). All of this is estimated for the current `article` format. Step 4's `bios-prelim.cls`
render replaces these estimates with real counts.

---

## 4. Items that need a writer or author decision before step 3 can be applied

1. The per-configuration table, the design-sample figure's caption, and the τ-coverage placement
   each involve new numbers or claims. Route them to phase-d-writer.
2. A2 has to be written fresh from the spec and code. The SI's S1, S3 and S4 are wrong
   post-revision.
3. The author should decide whether the appendix holds a curated subset (proposed above: about 6
   exhibits) or all of the SI's 22. This draft proposes curated: most SI exhibits are already body
   tables or prose, or are single-configuration duplicates.
4. The design-sample figure uses `DESIGN_FULL_NAMES`, which labels design 1 "Block Stratified
   Sampling (Checkerboard)". It needs either a six-design regeneration or a caption note.
5. **New deliverable (author, 2026-09-25):** Q&A backup slides for the prelim presentation, plus
   a written companion document answering likely committee questions, covering the full prelim
   (literature review, Chapter 2, Chapter 3). The A5 candidates above seed its Chapter 3 part.
   Scoped as its own task after step 7; it likely lives in bios-dissertation.
