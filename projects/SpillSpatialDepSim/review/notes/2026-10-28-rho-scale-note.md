# What ρ = 0.01 means in the Project 1 simulation (note, 2026-10-28)

Documentation only. The published and corrected results stand as reported
(author decision, 2026-10-28), and the deck presents them as they are.

## What the code does

In the frozen simulation scripts (2×4: `code/SpatialSim_NC_DOC.Rmd`, lines
380–471; the 3×3 and 3×4 scripts follow the same pattern), each replication
samples 20 subjects per cluster. Every subject in a cluster shares its
cluster's intervention indicator x, spillover indicator z and neighbor set, and
subjects differ only by independent error. This is why the manuscript writes
the model for the cluster-level outcome y_i (added after Reviewer 1, comment
R1.M1, asked why subjects were sampled at all).

The weights matrix used to *generate* outcomes is subject-level and binary:
subject i has w = 1 for every subject in a rook-neighboring district, so each
subject has about 40 (corner district) to 60 (middle district) neighbors, and
the rows are not standardized. Outcomes are
y = (I − ρW)⁻¹(α + βx + ψz + ε) with that W. The model is *fitted* with
`lagsarlm(..., listw = mat2listw(W, style = "W"))`, which row-standardizes the
same matrix.

So ρ = 0.01 on the generating scale corresponds to roughly
0.01 × 40–60 ≈ 0.4–0.6 on the row-standardized scale that the fitted ρ̂ and
the manuscript's w_ij = 1/|N(i)| use. The saved estimates agree: in
`results/rho_mse/rho_mse_results_TrtNoSpill_2x4_A02_B1_P05_R001.csv` the bias
of ρ̂ is about +0.48 for most allocations, so ρ̂ ≈ 0.49. "ρ = 0.01" therefore
means moderate spatial dependence, not near-independence.

## Quick check

`review/code/rho_scale_diagnostic.R` re-runs the 2×4 process (ψ = 0.5,
Control Clusters Only, all 70 allocations, N = 20 each):

| Generating process | SRS mean MSE | BSS mean MSE | mean ρ̂ | Share of allocations BSS beats |
|---|---|---|---|---|
| As published: binary W, ρ = 0.01 | 0.302 (published 0.304) | 0.253 (published 0.165) | 0.42 | 36% |
| ρ = 0 | 0.173 (published 0.173) | 0.288 | −0.01 | 1% |
| Row-standardized W in the generator, ρ = 0.5 (same W as the fit) | 0.176 | 0.374 | 0.48 | 1% |

The SRS values reproduce the published ones. The BSS value does not
(0.253 vs 0.165), probably because the checkerboard design is rank-deficient
(z = 1 − x) and its estimates are numerically unstable; this was not pursued.

## Reading

When the generating and fitted weights are the same matrix, adding spatial
dependence leaves the SRS error close to its ρ = 0 value and BSS does not gain
an advantage. In the published setting the two differ (binary, unequal
neighbor counts vs row-standardized), and that difference is a plausible source
of the larger SRS error at ρ = 0.01 and of BSS's relative advantage there. This
is one ψ value and N = 20, so it is a lead, not a result.

## How it is used

- The deck and documents report the published/corrected results unchanged.
- The deck's backup and speaker notes state the scale honestly if asked
  ("ρ = 0.01 on the subject-level generating weights, about 0.4–0.5 on the
  row-standardized scale of ρ̂").
- A fuller check (all ψ, N = 50, both spillover types, the BSS gap) would be
  needed before any claim in either direction; discuss with Dr. Lin.
