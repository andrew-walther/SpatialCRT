# Advisor meeting with Dr. Lin — 2026-10-14

Reviewed: the CTJ manuscript draft and the NC application study report, as of
2026-10-12. Sources: Andrew's prep and meeting notes, plus the Zoom AI summary. The
action items are in [TODO.md](TODO.md), under "Advisor meeting 2026-10-14 follow-ups".
Nothing below has been acted on yet.

## What Dr. Lin said

- **Big question: should the outcome be related to incidence?** An incidence-guided
  design makes sense only if the endpoint is correlated with the incidence used to guide
  allocation. The point of targeting high-SUD areas is to reduce SUD, so education
  effectiveness should be a proxy for reducing SUD, if the intervention works.
  - Analogy: poverty can guide stratification without being measured as the outcome.
  - Quoted Goodhart's law: "when a measure becomes a target, it ceases to be a good
    measure".
  - Open choices: incidence itself, or an intermediate endpoint (education
    effectiveness), as the outcome? Target the whole state, or only high-incidence areas?
- **Not every community college will participate.** How does the design guidance change
  with half the clusters (Zoom: 58 down to about 20)?
  - Fewer clusters means a smaller sample and possibly a non-contiguous study region.
  - Does incidence-guided saturation still win under control-only spillover with half the
    clusters?
  - Dr. Simpson will ask "how many clusters do we need?" Sample size and power results
    would be a concrete, useful answer.
- **Framing.** The application is an expanded, statewide version of Dr. Simpson's
  small-scale proposal. Results for lower participation strengthen it.
- **Simplify.** Group the allocation methods in Table 1 by type (e.g. deterministic,
  purely random, incidence-driven, saturation-driven). Consider combining the design
  table and an illustrative figure into one exhibit.
- **Pooling across years.** The current population-weighted pooling (Σ deaths / Σ
  person-years) is appropriate. An equal-weight average of the four yearly rates is an
  alternative, but is unlikely to differ much, since populations were stable over
  2018–2021. Checked afterward: no meaningful difference (correction 2).
  - Andrew's notes add: once incidence is known, a design can randomize further, e.g.
    incidence-guided saturation randomizing treatment and control within each region.
    That is how the design already works (saturation set by the region's incidence rank,
    then clusters randomized within the region), so this reads as a description of the
    current design, not a new task. Confirm with Dr. Lin if he meant more.
- **Spillover modeling.** Spillover is modeled as a single effect where it exists, not
  added up per treated neighbor. See correction 1: it is proportional to the share of
  treated neighbors here, unlike Project 1's on/off indicator.
- **Estimation sensitivity (a likely committee question).** How sensitive is the design
  guidance to the estimation method: SAR versus a conditional autoregressive (CAR) model
  or a Bayesian model? This is a robustness or sensitivity analysis.
- **Statistical tests of design differences.** These may be needed. Some are already in
  the supplement (see correction 6). Andrew (2026-10-19): the tests should show rigorously
  which designs are better or worse, together with the number of clusters needed to make
  that claim.
- **Keywords.** Update them; "mean squared error" may not be needed.
- **Prelim presentation.**
  - Use the Project 1 / 2 / 3 structure. That's more common than literature review first,
    then projects.
  - Go deep on the completed projects. Keep the literature review brief: it matters for
    the written document, not the talk.
  - Expect general-science questions from the committee.
- **Final defense.**
  - Don't revisit prelim material much. Focus on what changed between prelim and defense
    (mostly Project 3).
  - Expected status at prelim: Project 1 accepted, Project 2 submitted, Project 3 under
    way.

## Prep questions with no answer recorded

These went into the meeting but neither set of notes records an answer. Carry them
forward:

1. Formal sign-off on the application design choices settled on 2026-10-02:
   - β = 0 primary, with β = 1 as a sensitivity analysis. Now tied to the outcome
     question above.
   - Yearly planning surfaces.
   - Fixed geographic adaptations.
   - Unequal treatment budgets.
2. Guidance by spillover regime (no single winner) as the paper's framing.
3. Six designs plus SRS in the main text, with the two variants in the supplement.
4. Who holds the IRB protocol (Habib or Simpson), and permission to publish 58-cluster
   aggregates and maps.
5. Citing Habib's unpublished master's paper plus "personal communication". Dr. Lin found
   the personal communication "a bit weird", so an alternative citation is needed.
6. Funding statement ("no specific grant"): anything from Dr. Simpson's pilot or Habib's
   data work?
7. Target submission date, and a backup journal.
8. Undercoverage (about 93%): keep it as a stated limitation, or correct it? The
   t-critical-value check raised coverage only to about 93.7% (TODO.md).
9. Whether the deterministic designs (Checkerboard, High Incidence Focus) stay in, framed
   as comparators.
10. "Code available on request" versus a public deposit.
11. ENAR 2027 abstract: which project?

## Corrections to the notes (checked against the repository)

1. **Spillover form: proportional here, on/off in Project 1.** The code uses
   S(Z) = γ·WZ with row-standardized W, multiplied by (1 − Z) under control-only
   spillover (`application/code/real_sud_simulation.R:158`; same form on the grid).
   - So spillover equals γ times the *share* of a cluster's neighbors that are treated.
     Full γ applies only when every neighbor is treated; a control cluster with 2 of 4
     neighbors treated gets 0.5γ.
   - Project 1 (`SpillSpatialDepSim/code/SpatialSim_NC_DOC.Rmd`, ~lines 402–416) used a
     binary indicator instead: full ψ if *any* rook neighbor was treated, else 0.
   - "Single effect, not additive per neighbor" is right; "discrete" describes Project 1,
     not this project.
   - The chapter already explains the proportional form (row-standardized WZ) but
     doesn't contrast it with Chapter 2's on/off spillover. One sentence could be worth
     adding, keeping Chapter 2 links brief.
2. **Pooling: checked, no meaningful impact.** In the application, allocation uses each
   year's rates only (2018, 2019, 2020, 2021); no pooled-period run drives allocation.
   Pooled Σd/ΣP appears in three places:
   - the CTJ map (Figure 1A);
   - the saturation regions' reference population (pooled person-years / 4);
   - the "pooled deaths/population" regional-summary sensitivity analysis.

   The equal-weight average Dr. Lin raised is already computed as `mean_of_yearly_rates`
   in `application/data/derived/`. Checked 2026-10-14 against the 58 pooled cluster
   rates:
   - The two versions differ by −0.2% to +0.5%.
   - Rank correlation is 0.9999. Only 2 clusters move, by 1 rank each.
   - No cluster changes incidence quartile, and the top-29 set is identical.

   So the choice doesn't matter here, consistent with Dr. Lin's expectation.
3. **Word counts in the prep summary are stale** (3,425 / abstract 317). Recounted
   2026-10-14 with texcount (text + headers, math excluded):
   - Body 3,362, matching the title page.
   - Declarations 266, so 3,628 with declarations. The checklist's 3,526 predated the
     later wording changes. The 3,500 limit applies to the body (author rule).
   - Abstract 393, or 394 if the inline "10×10" counts (title page says 394).

   Recount after every CTJ edit and update the title page and checklist.
4. **Ethics statement.** The CTJ draft has a highlighted IRB placeholder in Methods (added
   2026-10-14), so the prep note "no ethics statement at all" is out of date.
5. **Zoom wording.** "OSHCON" is the Zoom transcription of "Ashkan" (confirmed by
   Andrew), so the line means discrepancies in Ashkan Habib's death data. That's the
   numerator discussion: his corrected counts (23,523) versus his paper's counts
   (21,147). Zoom also spells Ashkan "Ashkahn".
6. **Existing tests and power.**
   - The supplement already reports Nemenyi and Holm-adjusted Wilcoxon tests (used to fold
     the two variants into their counterparts).
   - `code/12_six_design_statistical_comparisons.R` has the six-design
     Friedman/Nemenyi/Wilcoxon analysis.
   - The grid study records power (reported in the supplement).
   - The application doesn't record power. The stored application fits contain τ̂ and SE
     for every outcome, so power at τ = 1 can be computed without new simulations.
7. **β = 1 sensitivity already spans the full primary queen grid in the application.**
   That covers every design, year, ρ, γ and regime. The grid study itself uses β = 1. So
   an "outcome correlated with incidence" version already exists in both studies; what's
   missing is an endpoint that *is* (future) incidence.

## Follow-up reasoning (after the meeting, 2026-10-14)

Reasoning only; nothing below has been simulated.

**Would all-or-nothing (Project 1-style) spillover change the results? Probably yes.**
- With queen neighbors (about 5–8 per cluster) and about half the clusters treated,
  almost every control cluster has at least one treated neighbor.
- An on/off spillover term would then nearly equal "is a control" under control-only
  spillover, or a constant under both-arms spillover. τ would be poorly identified.
- Only control clusters with no treated neighbors separate treatment from spillover.
- Expected by design:
  - **Isolation Buffer:** not identified. The maximal independent set guarantees every
    untreated cluster touches a treated one, so spillover = γ(1 − Z) exactly.
  - **Checkerboard:** the same problem.
  - **SRS and the balanced designs:** worse, because few controls are fully isolated.
  - **Saturation designs:** likely to gain, since the 20% region leaves some isolated
    controls.
- So the ranking would probably shift further toward saturation.
- The proportional form (γ × share of treated neighbors) is identified because controls
  vary in their treated-neighbor share. It is also the more plausible exposure model.
- This is a candidate sensitivity analysis, or a prepared answer to "why not Project 1's
  binary spillover?"

**What β means, and how it depends on the outcome**
- In Y = (I − ρW)⁻¹[τZ + S(Z) + βX + ε], β sets the outcome's *baseline level* by
  incidence. It does not change the treatment effect: τ is constant across clusters for
  every β.
  - **β = 0:** a cluster's incidence has no effect on its outcome. Incidence only guides
    which clusters the designs treat.
  - **β = 1:** higher-incidence clusters have a higher outcome regardless of treatment.
- The right β depends on what Y is:
  - **Y is an education measure** (e.g. professionals' knowledge): β = 0 is plausible.
    This is the current primary framing.
  - **Y is future SUD incidence, reduced by the education:** past incidence predicts
    future incidence, so β > 0 is natural. This is Dr. Lin's "outcome correlated with
    the guiding incidence".
- Under the second framing, the existing NC β = 1 runs approximate the setting:
  rank-scaled incidence in the outcome, adjusted for in the fit, like a baseline
  adjustment. Caveats: the outcome is continuous rather than counts, and the size β = 1
  is arbitrary.
- Support for reporting both ways (option d):
  - The grid uses β = 1, and NC uses β = 0 with a β = 1 sensitivity. The regime
    conclusions hold in both.
  - In NC under β = 1, most designs change by only a few percent. High Incidence Focus
    is the exception: its MSE rises 1.4–4.4×.
- Not modeled anywhere: a treatment effect that grows with incidence (larger reductions
  in high-burden areas, e.g. a Z·X interaction). That would be option (e) if the framing
  goes that way.

## Synthesis: decisions needed before more work

- **The outcome question comes first.** It determines what the remaining analyses and
  edits look like. Three options, from least to most work:
  - (a) Keep β = 0 primary, and state explicitly that the endpoint is assumed to track
    SUD reduction. Point to the grid study (β = 1) and the application's β = 1
    sensitivity as the correlated-outcome case.
  - (b) Promote the existing β = 1 runs to co-primary or primary. This is a reporting
    change: no new simulation, but it reverses the 2026-10-02 rule (DECISIONS.md) and
    changes numbers across the chapter, CTJ and supplement.
  - (c) Model a future-incidence endpoint (e.g. a count outcome). This is a new model and
    new runs.
  - (d) Report both ways, β = 0 and β = 1 side by side (Andrew's leaning, 2026-10-14).
    This needs no new simulation, since both run sets exist. The cost is more exhibits or
    space, and a framing choice about which version leads.
  - (e) Model a treatment effect that grows with incidence (Z·X interaction). This is a
    new model and new runs; τ is currently constant by author rule.

  The choice among these hinges on what Y is: an education measure points to β = 0,
  future SUD incidence to β > 0. See "Follow-up reasoning" above.

  Not decided yet: the choice depends on how we frame the incidence–outcome
  relationship and needs careful thought. Record the outcome in DECISIONS.md.
- **Partial participation and power** are the main new analyses.
  - Partial participation means sampling subsets of colleges (about half, or about 20),
    rebuilding W on the possibly non-contiguous subset, and re-running the comparison.
  - Power is cheap for the current runs and needed again for any subset runs.
  - Andrew, 2026-10-14: sample size/power and reduced participation (fewer counties, so
    non-adjacent clusters) must be addressed in the thesis chapter, and so in the prelim
    talk (main slides or backup). Sample size/power results may also move into the main
    CTJ manuscript.
- **Manuscript edits route through the chapter.** `Thesis_Chapter.qmd` is canonical; the
  CTJ and supplement are derived from it. Run both verifiers after edits.

## Planned extensions (added 2026-10-19)

The application report (Section 8, "Planned extensions") now tells collaborators about
three planned analyses: partial participation, sample size and power, and formal design
comparisons. Points behind that plan (reasoning, nothing run yet):

- **The outcome decision comes first.** Power needs an effect size on the outcome's scale.
- **New simulations are needed.** All application runs use τ = 1 with σ = 1, a 1-SD
  effect, so power with 58 clusters is likely near 100% and uninformative. Power curves
  need smaller effects (e.g. 0.2–0.8 SD) and fewer clusters.
- **Simulation replicates are not trial sample size.** Paired tests across settings become
  significant for tiny differences if Monte Carlo replicates are added. Present the tests
  as "consistently better across settings", and use power curves (in colleges) to answer
  "how many colleges to make that claim".
- **Partial-participation gotchas:**
  - A participating college may have no participating neighbors, making it an island
    with zero spatial weight.
  - Saturation regions may have few participants.
  - Non-participating areas are never treated, so they create no spillover.
  - Still to decide: whether participation is random or informative (e.g. urban colleges
    joining more often).

