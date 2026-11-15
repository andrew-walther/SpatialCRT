# Prelim Part II: question prep for Project 2

Author: Andrew Walther. Started 2026-11-15 from the Chapter 3 committee review (round 5), for
the oral prelim on 2026-12-07.

This file collects questions a committee member is likely to ask about Project 2 and short
answers drawn from Chapter 3, the CTJ manuscript and the code. Use it to decide which points
need a main slide, which need a backup slide, and which only need a prepared spoken answer.
Numbers come from the stored results as checked in the committee review
(`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/committee-proposals.md`,
items Q-01 to Q-15 and the C items cited below). Section names refer to Chapter 3 after the
round 5 edits. Recheck any number against the chapter before it goes on a slide.

## How to use this list

- **Main slide**: the answer is part of the core argument, so the talk should cover it
  before anyone asks.
- **Backup slide**: a figure or table helps answer it, but it is outside the main story.
- **Spoken answer**: a prepared sentence or two is enough.

| # | Question (short) | Suggested format |
|---|---|---|
| Q-01 | What is β? Is it the effect a trial would report? | Main slide (model) plus a spoken answer |
| Q-02 | Why do saturation designs beat SRS under control-only spillover? | Main slide |
| Q-03 | Why use incidence at all? | Main slide (take-home) |
| Q-04 | High Incidence Focus: four times worse than SRS on the grid, slightly better in NC? | Backup slide |
| Q-05 | What if the spillover term is wrong? | Backup slide |
| Q-06 | Why a t interval with N − p degrees of freedom? | Backup slide |
| Q-07 | The trial team doesn't know the spillover regime. What should they do? | Main slide (take-home) |
| Q-08 | Chapter 2 used rook neighbors. Why is queen primary here? | Spoken answer, backup for the rook results |
| Q-09 | ρ̂ is biased downward. Does that matter? | Backup slide |
| Q-10 | Why does power fall when more colleges join? | Backup slide |
| Q-11 | NC gives about 92 deaths per 100,000, but the grid uses 35. Why? | Spoken answer |
| Q-12 | Why not covariate-constrained randomization on incidence? | Spoken answer |
| Q-13 | Are ten incidence surfaces per configuration enough? | Spoken answer |
| Q-14 | Why is the CAR model biased when the Bayesian spatial lag model isn't? | Backup slide |
| Q-15 | Isolation Buffer treats 21 clusters and the others 50. Is that fair? | Spoken answer |
| Q-16 | Why saturations of 80/60/40/20% and four regions? | Spoken answer, possibly a future-work slide |
| Q-17 | Why is ψ tied to β in the participation analysis? | Spoken answer |
| Q-18 | Was the interval method chosen on the same fits used to report its coverage? | Spoken answer with Q-06's backup |
| Q-19 | Does the code fit the spillover term as the text describes? | Spoken answer |
| Q-20 | Does the recommendation hold if only 20 to 30 colleges join? | Main slide (take-home) |

## Questions and suggested answers

### Q-01. What exactly is β? Is it the effect a trial would report?

β is the coefficient on a cluster's own treatment before the spatial multiplier
(I − ρW)⁻¹ spreads outcomes between neighbors (Methods, Spatial Weights and Outcome
Model). It is not a marginal potential-outcome contrast or an average total effect. With
ρ > 0, a cluster's total response to its own treatment is [(I − ρW)⁻¹]ᵢᵢβ, slightly more than
β, plus spillover through ψ. The chapter compares designs on how well they estimate this one
parameter when the model is correctly specified. An average total or direct effect could be
derived from the fitted model, but the design comparison for those estimands was not studied.

*Slide note:* the model slide should label β as the direct (structural) effect and say in one
line how it differs from a total effect.

### Q-02. Why do the saturation designs beat SRS when spillover reaches only controls?

To estimate β, the model has to separate a cluster's own treatment from its spillover
exposure. Under control-only spillover the exposure is (1 − Z)∘WZ, and the information that
separates ψ from β comes from controls with different shares of treated neighbors. Saturation
creates regions with 20% and 80% treated, which gives controls a wide range of exposure.
Designs that fix the local share (Checkerboard, 2×2 Blocking) give almost every control the
same exposure, which is why they do poorly in this regime (Results, Design Performance by
Spillover Regime; Appendix A2). Saturation Quadrants does as well as Incidence-Guided
Saturation, which shows that the gain comes from the saturation structure, not from incidence.

*Slide note:* a schematic of control exposure under saturation versus Checkerboard answers
this better than a table.

### Q-03. Then why use incidence at all?

Incidence guidance costs no accuracy: Saturation Quadrants had MSE 0.090 and
Incidence-Guided Saturation 0.091 (Nemenyi p = 0.99). It also sends the highest saturations
to the highest-incidence regions. How much it redirects depends on how much regional
incidence differs: 2.1 to 3.1 percentage points of incidence share on the grid, and 0.7 points
on the NC map. Its value is in where treatment goes, not in estimation accuracy. Be ready to
concede that on the NC map the targeting gain is small.

### Q-04. High Incidence Focus is four times worse than SRS on the grid but slightly better in NC. Which do we believe?

Both, for different reasons. On the grid, τ = 1 and incidence is in the fitted model. Treating
the top half by incidence makes Z nearly collinear with X, which leaves little information
about β after adjustment. In the primary NC model τ = 0 and incidence is not fitted. When
incidence was allowed to shift the NC outcome and was adjusted for, its MSE rose 3.7 to 4.4
times under both arms. It is also a deterministic allocation, so its accuracy rests on the
model, not on randomization (committee item C-08, now explained in the Discussion).

*Backup slide:* grid ratio (3.81) next to the NC primary ratio (0.93) and the NC matched
sensitivity result.

### Q-05. The primary estimator is an oracle. What happens when the analyst gets the spillover term wrong?

Omitting the term biased every design (bias −0.121 to −0.400, coverage 0.536 to 0.907,
Appendix A10), and Incidence-Guided Saturation still had the lowest MSE. In NC, binary and
count exposure forms in both the generating and the fitted models kept the saturation
advantage under control-only spillover. A wrong form (fitting proportional exposure when the
truth is binary, or the wrong neighbor definition) was not studied, and this is the first
limitation.

### Q-06. Why a t interval with N − p degrees of freedom? Where does it come from?

With 58 or fewer clusters the Wald interval under-covered (92.6% to 93.6% across the nine
designs in NC), and the Wald test rejected a true null in 11% to 13% of trials with 15 colleges.
We tried a t quantile alone, the degrees-of-freedom factor alone, both together, a parametric
bootstrap, a bootstrap-t and profile likelihood, and adopted the method that met a 94% to 96%
coverage rule (Appendix A12). The factor √(N/(N − p)) is the divisor change used for the
unbiased residual variance in linear regression and in the HC1 estimator (MacKinnon and White
1985). In the spatial lag model it is an approximation, and its justification is the
calibration. The t quantile follows the small-sample advice of Leyrat et al. (2018).

*Backup slide:* the interval-comparison table (Wald versus t, six designs and SRS).

### Q-07. A trial team usually doesn't know the spillover regime. What should they do?

Among the six designs, Incidence-Guided Saturation was first or second under both regimes.
Choosing it costs about 8% in MSE if spillover reaches both arms (grid, 95% interval 5.6% to
11.0%) and saves about 37% if it reaches only controls (34.5% to 39.2%). In NC the both-arms
cost was within 3%. SRS is safe under both arms but loses about 37% relative to saturation
under control-only spillover. The asymmetry favors saturation unless both-arms spillover is
expected. Pair this with Q-20: the advantage shrinks when few colleges join.

### Q-08. Chapter 2 used rook neighbors. Why is queen primary here?

Under rook contiguity on the rectangular grid, Checkerboard gives WZ = 1 − Z exactly, β is not
identifiable, and the coefficient on Z estimates β − ψ. Queen contiguity keeps Checkerboard
identifiable, so all six designs can be compared on the same estimand, and corner contacts are
plausible exposure routes between service areas. Rook results are reported in full as a
sensitivity analysis, with the same two leading designs (Appendix A9).

### Q-09. ρ̂ was biased downward. Does that affect the conclusions?

In one place. The error in ρ̂ enters β̂ through the term β(ρ − ρ̂)g (Appendix A3). For designs
whose treated clusters have treated neighbors, g > 0, and the term adds variance that grows
with β. Isolation Buffer has g = 0, which is why it passes Incidence-Guided Saturation under
control-only spillover at β = 3. A finite-sample bias of this kind is known for the spatial
autoregressive estimator (Bao 2013). Squared bias stayed at or below 3% of Incidence-Guided
Saturation's MSE.

*Backup slide:* the variance of the last term against β for the two designs (0.008 to 0.061
for Incidence-Guided Saturation, below 0.0002 for Isolation Buffer).

### Q-10. Why does power fall when more colleges join, for example Checkerboard under control-only spillover?

Colleges that don't take part count as untreated in the exposure. With few participants,
controls differ a lot in their share of treated neighbors, which separates spillover from
treatment. When nearly every college takes part under an interspersed design, nearly every
control has the same exposure, and the treatment and spillover columns become nearly
collinear. For Checkerboard at β = 1.5, power is 0.72 to 0.73 with 29 colleges and 0.38 to
0.39 with 58 (`application/results/methods_ext_2026-10/summaries_v2/a3_n80_by_year.csv`). The
same mechanism explains why the saturation advantage shrinks from 22% at 58 colleges to
about 1% at 20.

### Q-11. The North Carolina data give about 92 deaths per 100,000. Why does the grid use 35?

23,523 / 25,594,321 × 100,000 = 91.9 per 100,000 person-years at ages 18 to 64, with the
corrected numerator that includes heart-failure deaths. The grid base rate of 35 comes from
published US rates (Mirzaei et al.) and was not meant to match NC. In the Poisson mechanism
the rates are rank-normalized, so the base rate only affects the number of ties. The case
study uses the observed NC rates directly. Be ready to explain why the corrected NC rate is
higher than published rates (case definition, the heart-failure deaths, the period).

### Q-12. Why not compare against covariate-constrained randomization on incidence?

Balanced Quartiles and Balanced Halves are the stratified version of that idea, and they
matched SRS on both criteria under both regimes. Constrained randomization targets balance on
incidence, which does not control the spread of spillover exposure, and that spread is what
drives accuracy here. A constrained design that also restricts exposure patterns is a
possible extension.

### Q-13. Why ten incidence surfaces per configuration? Is that enough?

Surfaces are the primary sampling unit, and all standard errors and paired tests use the 50
independent configuration-by-surface units. The pooled MSE standard errors are below 0.005 for
the leading designs, each design's MSE range across configurations lies below the next
design's minimum, and adjacent designs differed in the same direction in at least 48 of 50
units.

### Q-14. Why is the Bayesian CAR model biased when the Bayesian spatial lag model isn't?

The data come from the spatial lag model, whose mean is (I − ρW)⁻¹(α + βZ + ψS + …). A CAR
model puts the spatial structure in the errors and fits the mean α + βZ + ψS, so its
coefficient on Z is a different quantity when ρ > 0 (Wall 2004). The maximum likelihood CAR fit
to the same data shows the same bias: 0.110 under both-arms and 0.205 under control-only
spillover at ρ = 0.5, and at most 0.005 when the CAR model is fit to CAR data. The bias is a
difference in estimand, not a failure of the Bayesian fit.

### Q-15. Isolation Buffer treats 21 clusters and the others treat 50. Is the comparison fair?

It compares the designs as specified, not at equal cost (Limitations). Isolation Buffer's
treated count is set by the graph: under queen contiguity on a 10 × 10 grid, no more than 25
clusters can be treated without two treated clusters touching. Its MSE does not depend on
the regime, since WZ = 0 at every treated cluster. Its both-arms MSE was about 3.8 times that
of SRS, and treating fewer clusters is one likely reason.

### Q-16. Why saturations of 80/60/40/20% and four regions?

The designs came from design ideas developed for this trial, not from an exhaustive search.
The saturation designs used one set of saturations and four regions, and other sets and
numbers of regions (Baird et al. 2018 treat this choice as a design problem) were not studied.
Since this version already outperformed SRS under control-only spillover, varying the
saturations and the number of regions is future work. With the existing simulation code,
each alternative is one more design ID on the grid.

### Q-17. Why is ψ tied to β in the participation analysis?

The participation analysis sets ψ = κβ with κ ∈ {0.5, 0.8}, while the main study fixes ψ
separately. Tying them means an intervention with no effect also has no spillover, which the
size calculation at β = 0 needs, and it keeps the ratio of spillover to direct effect the same
along each power curve.

### Q-18. Was the interval method chosen on the same fits used to report its coverage?

Yes for the main case-study fits, and the chapter says so. The participation runs are new
simulations that played no part in the choice. In them the interval had size 4.4% to 5.6% and
coverage 94.0% to 95.6% at every number of participating colleges, which is the independent
check.

### Q-19. Does the code fit the spillover term as the text describes?

The grid code passes ψ·S(Z) as the covariate, with true coefficient 1, rather than S(Z) with
coefficient ψ (`code/05_run_simulation.R`, around lines 296 and 310). Rescaling a covariate
doesn't change β̂, its standard error or any reported metric, so every reported result is the
same either way. The participation code passes the unscaled exposure.

### Q-20. Does the recommendation hold if only 20 to 30 colleges join?

The main analysis assumes that most colleges take part. With 30 or fewer participating
colleges, the saturation designs' advantage over SRS was 5% or less, about 1% with 20 colleges,
and 2% or less when high-incidence colleges volunteer more often. If the statewide extension
expects 20 to 30 colleges, the honest message is that the design choice matters little for
accuracy at that size, and the choice can rest on targeting and logistics.

## Topics to think about beyond the written chapter

- How the estimand β relates to what the trial will report (Q-01).
- What the trial team expects about participation and spillover, which decides which result
  applies (Q-07, Q-20). Ask Dr. Simpson before the prelim if possible.
- The limits of the oracle assumption and what a misspecified exposure would do (Q-05).
- Where the project goes next: other saturation sets, other exposure forms, and estimands
  other than β (Q-01, Q-05, Q-16).
