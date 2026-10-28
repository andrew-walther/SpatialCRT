# Grid allocation risk on all ten surfaces (surfaces10): findings

Date: 2026-10-27. Author decision: DECISIONS.md, 2026-10-27 (rerun on all ten surfaces).

**No conclusion changes.** Every direction the manuscripts draw from the two-surface
pilot holds on all ten surfaces. Under control-only spillover, both saturation designs
and Isolation Buffer beat SRS on mean and worst-decile MSE in all 200 fixed units. Under
both arms, the saturation designs stay slightly above SRS, Balanced Quartiles and Balanced
Halves tie it, and the SRS ordering does not flip. The largest movement is High Incidence
Focus under control-only spillover: its mean ratio rises from 1.38 to 1.51 and its
worst-decile ratio from 0.87 to 0.96 (still below 1). Its both-arms mean ratio falls
from 3.91 to 3.81.

## What ran

| Step | Command (from `code/` unless noted) | Output |
|---|---|---|
| Grid run | `Rscript 16_allocation_risk.R surfaces10 8` (single-thread BLAS) | `results/allocation_risk/surfaces10/` |
| R = 400 precision | `Rscript 18_allocation_risk_precision.R surfaces10` | `results/allocation_risk/surfaces10_precision_R400/` |
| Summaries | `Rscript code/17_allocation_risk_summary.R surfaces10` and `... surfaces10_precision_R400` (project root) | `results/allocation_risk/summary/surfaces10/`, `.../summary/surfaces10_precision_R400/` |
| Downside table and Figure 3 | `Rscript code/21_downside_by_design.R` (project root; default profile surfaces10) | `results/allocation_risk/summary/surfaces10/downside_by_design.csv`, `allocation_distribution_by_design.{csv,pdf}`; PDF copied to `paper/ctj_manuscript/figures/compact/` and `paper/thesis_chapter/figures/compact/` |

The `surfaces10` profile is the pilot with `surfaces = 1:10`: queen, τ = 1,
ρ ∈ {0, 0.5}, γ ∈ {0.5, 0.8}, both regimes, nine designs, J = 100 allocations,
R = 100 outcomes per distinct allocation, and the same seed namespace
(`allocation-risk-v1`).

**Cost.** 3,600 blocks; 280,912 distinct allocations; **28,091,200 outcome fits**, all
valid (0 aliased, 0 warnings). That took 5,123 s (85 min) on 8 workers. The pilot was
720 blocks, 5,617,600 fits and 777 s. The R = 400 step re-fit 50 blocks: 2,000,000 fits
(1,500,000 beyond R = 100) in 179 s. In all, the substantive analysis comes to
**29,591,200 distinct outcome fits**. There are 339 signed negative corrected-variance
estimates in the main run (pilot: 65) and 2 in the R = 400 step. They are retained and
counted, as before.

## Checks

- **Pilot reproduction.** All 720 pilot blocks (surfaces 1–2) match the surfaces10
  blocks with the same keys. The largest absolute difference over every numeric summary
  column is 0. All 56,176 pilot allocation rows match by assignment hash, and their MSEs
  are identical. The comparison table is
  `summary/surfaces10/pilot_unit_reproduction.csv`.
- **Precision.** For every design × regime, the MCSE of the mean MSE is at most 1.01% of
  the value (Checkerboard control-only; every randomized design is ≤ 0.21%), against the
  5% target. The coverage MCSE is at most 0.0016 against the 0.01 target. Coverage MCSE
  uses the same joint formula as `SE_Mean_MSE_Joint`, applied to per-allocation coverage
  and pooled as sqrt(Σ SE²)/B.
- **SRS self-ratio.** Exactly 1 for every metric and regime.
- **Direction vs the pilot.** Control-only Incidence-Guided Saturation/SRS stays at 0.647
  for the mean (pilot 0.649) and 0.552 for the worst decile (pilot 0.555). It is below
  SRS in 200/200 units on both measures.
- **R = 400 step.** The task kept the pilot's selection rule: surface 1, ρ = 0.5,
  γ = 0.8, designs {3, 4, 6, 8, 9}, giving 50 blocks. All of those blocks are pilot
  surfaces, so the extension reproduces `precision_R400/` exactly. It adds no new
  information beyond confirming that the R = 400 path is reproducible. Widening the
  rule to more surfaces would be a separate decision.
- **Code changes stay backward-compatible.** `21` in pilot mode, run in a scratch copy,
  rewrites the pilot `downside_by_design.csv` and `allocation_distribution_by_design.csv`
  byte-identically. `17` on `precision_R400` rewrites all its pilot summary CSVs
  byte-identically.

## Old vs new: numbers the manuscripts cite

Pooled over fixed units (pilot 40, surfaces10 200 per regime);
`summary/<profile>/pooled_fixed_units.csv` and `downside_by_design.csv`.

| Quantity | Pilot | surfaces10 |
|---|---:|---:|
| IGS control-only mean MSE vs SRS | 0.1409 vs 0.2172 | 0.1402 vs 0.2165 |
| IGS control-only q90 vs SRS | 0.1747 vs 0.2980 | 0.1733 vs 0.2962 |
| IGS control-only worst-decile mean vs SRS | 0.1936 vs 0.3486 | 0.1912 vs 0.3462 |
| IGS both-arms mean MSE vs SRS | 0.0458 vs 0.0428 | 0.0457 vs 0.0427 |
| IGS both-arms worst-decile mean vs SRS | 0.0582 vs 0.0544 | 0.0583 vs 0.0544 |
| Mean-MSE ratio IGS/SRS, control-only | 0.649 | 0.647 |
| SRS control-only worst decile vs average SRS allocation | 0.3486 vs 0.2172 | 0.3462 vs 0.2165 |
| Isolation Buffer control-only worst-decile ratio | 0.62 | 0.62 |
| Isolation Buffer both-arms mean ratio | 3.79 | 3.80 |
| High Incidence Focus both-arms mean ratio | 3.91 | 3.81 |
| Saturation designs' worst-decile ratios, control-only (SQ, IGS) | 0.53, 0.56 | 0.54, 0.55 |
| Saturation designs' worst-decile ratios, both arms (SQ, IGS) | 1.05, 1.07 | 1.06, 1.07 |
| Saturation designs' mean ratios, both arms (SQ, IGS) | 1.06, 1.07 | 1.06, 1.07 |
| Grid fixed units per regime | 40 | 200 |
| Surfaces per incidence configuration | 2 | 10 |

### Full grid columns of `tab:downside` (ratios to SRS, two decimals)

| Design | Both mean | Both worst 10% | Control mean | Control worst 10% |
|---|---|---|---|---|
| Checkerboard | 1.41 → 1.41 | 1.10 → **1.11** | 9.80 → **9.94** | 6.11 → **6.22** |
| High Incidence Focus | 3.91 → **3.81** | 3.11 → **3.05** | 1.38 → **1.51** | 0.87 → **0.96** |
| Saturation Quadrants | 1.06 → 1.06 | 1.05 → **1.06** | 0.62 → 0.62 | 0.53 → **0.54** |
| Isolation Buffer | 3.79 → **3.80** | 3.93 → **3.97** | 0.75 → 0.75 | 0.62 → 0.62 |
| 2×2 Blocking | 1.44 → 1.44 | 1.57 → 1.57 | 2.68 → **2.69** | 2.63 → **2.64** |
| Balanced Quartiles | 0.99 → 0.99 | 1.00 → **0.99** | 1.00 → 1.00 | 1.00 → **1.01** |
| Balanced Halves | 0.99 → 0.99 | 1.00 → **0.99** | 1.00 → **1.01** | 0.98 → **1.01** |
| Incidence-Guided Saturation | 1.07 → 1.07 | 1.07 → 1.07 | 0.65 → 0.65 | 0.56 → **0.55** |
| Simple Random Sampling | 1.00 | 1.00 | 1.00 | 1.00 |

The NC columns are unchanged: they come from the application, and their rows in the new
CSV are byte-identical.

### Tail-stability diagnostics (medians over fixed units)

| Slice | Regime | IGS Jaccard | SQ Jaccard | SRS Jaccard | IGS noise share | SRS noise share |
|---|---|---:|---:|---:|---:|---:|
| R100, all units (pilot → surfaces10) | control-only | 0.111 → 0.111 | 0.111 → 0.111 | 0.250 → 0.250 | 0.617 → 0.628 | 0.252 → 0.263 |
| R100, all units (pilot → surfaces10) | both | 0.053 → 0.053 | 0.053 → 0.053 | 0.053 → 0.053 | 0.938 → 0.924 | 0.960 → 0.932 |
| R400 subset (identical in both) | control-only | 0.176 | 0.176 | 0.538 | 0.382 | 0.083 |
| R400 subset (identical in both) | both | 0.053 | 0.111 | 0.111 | 0.828 | 0.753 |

"Noise share" is `Variance_Noise_Fraction`, the share of the raw across-allocation
variance of estimated risk that comes from outcome Monte Carlo noise. The pilot findings'
0.176/0.538 Jaccard and 0.382/0.083 noise shares come from the R400 subset, which is
unchanged. The manuscripts' "Jaccard 0–0.25" and "61–86%" are NC numbers, not grid
numbers.

### Figure 3 (grid panels)

The figure was regenerated in the same style, with no in-figure title. The grid panel
barely moves. For control-only IGS, the median ratio goes 0.634 → 0.635 and q95 goes
0.865 → 0.856, so most IGS allocations still sit below the average SRS allocation. For
SRS, q95 goes 1.521 → 1.520. Both-arms High Incidence Focus moves from 3.91 to 3.81, and
control-only from 1.38 to 1.51.

## Every pilot citation in the manuscripts

From `grep -n -i pilot`, run on 2026-10-27 before the other session's edits. Line numbers
may have shifted since. "education pilot" and the NC "after the pilot" / "a real pilot"
lines are excluded.

| File:line | What it cites | New value / action |
|---|---|---|
| `CTJ_Manuscript.tex:51` | Abstract: "A focused pilot simulated repeated outcomes…" | wording only |
| `CTJ_Manuscript.tex:181` | "two fixed surfaces per configuration … supporting pilot evidence" | ten surfaces; 100 allocations × 100 outcomes, 400 for selected |
| `CTJ_Manuscript.tex:220` | IGS worst decile 0.1936 vs SRS 0.3486; both arms 0.0582 vs 0.0544; "pilot values" | 0.1912 vs 0.3462; 0.0583 vs 0.0544 |
| `CTJ_Manuscript.tex:232` | IB control-only worst-decile 0.62 "in the grid pilot"; IB both-arms mean "3.7–3.8"; HIF "grid-pilot both-arms mean 3.91" | 0.62; 3.80 grid / 3.69 NC (range still 3.7–3.8); 3.81 |
| `CTJ_Exhibits.tex:86` | Fig. 3 caption: "focused grid pilot … 40 fixed surface/parameter settings per spillover regime" | 200 settings (10 surfaces × 5 configs × 4 corners) |
| `Supplementary_Information.tex:23` | "the selected-surface allocation-risk pilot" | wording only |
| `Supplementary_Information.tex:353` | section title "…allocation-risk pilot" | wording only |
| `Supplementary_Information.tex:989–992` | "focused grid pilot … two fixed surfaces … The pilot used 100 assignment draws…" | ten surfaces |
| `Supplementary_Information.tex:994–1004` | 40 units; IGS ctrl mean 0.1409 vs 0.2172, q90 0.1747 vs 0.2980, worst 0.1936 vs 0.3486; both mean 0.0458 vs 0.0428, worst 0.0582 vs 0.0544; "not the full-study grid estimates"; "from this pilot" | 200 units; 0.1402/0.2165, 0.1733/0.2962, 0.1912/0.3462; 0.0457/0.0427, 0.0583/0.0544 |
| `Supplementary_Information.tex:1247` | IB 0.62 "in the grid pilot"; IB both mean 3.79; HIF "In the grid pilot … 3.91" | 0.62; 3.80; 3.81 |
| `Supplementary_Information.tex:1255` | `tab:downside` caption: "Grid: focused pilot … averaged over 40 units" | 200 units |
| `Supplementary_Information.tex:1262` | `tab:downside` header "Grid pilot" | grid columns per the table above |
| `Supplementary_Information.tex:1288` | Fig. caption "focused grid pilot … 40 fixed surface/parameter settings" | 200 |
| `Thesis_Chapter.qmd:425` | same as SI 989–992 | ten surfaces |
| `Thesis_Chapter.qmd:427` | same numbers as SI 994–1004 | as SI |
| `Thesis_Chapter.qmd:575` | same as SI 1247 | as SI |
| `Thesis_Chapter.qmd:584` | `tab:downside` caption, 40 units | 200 |
| `Thesis_Chapter.qmd:591` | header "Grid pilot" + grid table cells (lines ~595–603) | grid columns per the table above |
| `Thesis_Chapter.qmd:619` | Fig. caption, 40 settings | 200 |
| `Thesis_Chapter.qmd:666` | "in the grid pilot, the worst 10% of SRS allocations … 0.3486, against 0.2172" | 0.3462 vs 0.2165 |
| `Thesis_Chapter.qmd:674` | "about four times worse than SRS under both arms in the grid pilot" (HIF); "saturation designs had mean MSE 4% to 8% above SRS on the grid"; "worst-decile MSE within 7% of SRS's" | 3.81 ("about four times" still fits); grid allocation-risk both-arms means are 1.06/1.07 (6–7% above). Check whether "4–8%" comes from the full grid benchmark before changing it. Worst decile 1.06/1.07 (still within 7%) |

Not pilot-sourced (left alone): `CTJ_Manuscript.tex:232` and SI/chapter "112 settings,
Jaccard 0–0.25, 61–86%" (NC tail refinement); SI 1741/1750 and chapter 1173/1182
"−0.176" (an unrelated table cell).

## Follow-ups for the main session

- `paper/tools/verify_manuscripts.py:223–227` reads
  `results/allocation_risk/summary/downside_by_design.csv` with source `grid_pilot`. Point
  it at `results/allocation_risk/summary/surfaces10/downside_by_design.csv` with source
  `grid_surfaces10` once the `tab:downside` cells are updated.
- `application/report/IncidenceDesign_Application_Study_Report.qmd:95` also reads the
  pilot `downside_by_design.csv`.
- Re-running the **pilot** profile through `16` or `18` would now stop on its manifest
  check, because adding the `surfaces10` profile changed the source hash of
  `16_allocation_risk.R`, and `18`'s own hash changed. The pilot results and summaries on
  disk are untouched and still read correctly (`17` and `21` reproduce them
  byte-identically).
