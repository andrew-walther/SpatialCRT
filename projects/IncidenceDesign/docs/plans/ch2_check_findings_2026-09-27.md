# Chapter 2 check and accepted-manuscript review — 2026-09-27

## Status and author direction

This is the durable continuation note for the minimal Chapter 2 check in
`handoff_2026-09-27_srs_ch2_audit.md`. Read it before deciding the Chapter 3 SRS
framing, Checkerboard placement, or Chapter 2 link. It supplements the earlier
draft-only discussion in `paper/dissertation_chapter/notes/a5_technical_notes.md` Q3.

- The author regards Project 1 as accepted/final and does not believe its findings
  can be updated. Keep SpillSpatialDepSim and bios-dissertation read-only; do not
  revise their results, figures, tables or prose as part of this work.
- The author asked that the accepted manuscript's allocation-level figures and
  tables be considered, rather than interpreting one reproduction check alone.
- On 2026-09-27 the author authorized saving these findings in **IncidenceDesign's
  `docs/plans/`**, so future Codex/Claude sessions can resume without relying on
  temporary files. This authorizes documentation here, not changes to Project 1.
- No Chapter 3 or CTJ manuscript changes, design exclusions, or framing decisions
  have been approved or implemented in this exchange. Further audit work is paused
  for discussion, not an automatic prerequisite that the author has agreed to pursue.

## Narrow audit verdict and its limits

**Verdict (c) for the requested numerical check:** the reported 3×4 control-only
checkerboard MSE of 0.0004 does not reproduce from the current original simulation
script. The rerun agrees with the stored CSV, at BSS MSE approximately **0.22555**.

This establishes a discrepancy between current code/stored results and the accepted
exhibit. It does **not** establish the historical source of that discrepancy or
verify that the current script/CSV was the exact input used for the final submitted
figure and table. Do not turn this limited finding into a claim that the whole
accepted study has been invalidated, or a settled explanation of the difference
between Chapters 2 and 3. The first audit answer preceded review of the accepted
PDF; the PDF review below is essential context for interpreting it.

## What the code check established

Line references below are to `../SpillSpatialDepSim/` relative to IncidenceDesign,
as inspected on SpatialCRT commit `051ceb1`.

| Question | Finding and evidence |
|---|---|
| Is W subject-level? | **Yes.** `code/SpatialSim_3x4.Rmd:481–498` builds a 240×240 binary matrix, linking each subject to all subjects in rook-neighboring districts. |
| Is the fit subject-level? | **Yes.** Lines 515–519 convert W to row-standardized weights and fit `response ~ intervention + spillover` to individual responses. |
| How are x and z assigned? | Lines 451–466 assign treatment and spillover by district. Both are constant within a district. |
| Are 313 and 612 checkerboards? | **Yes.** `combn(1:12,6)` at 367–368 gives treated districts {1,3,6,8,9,11} and {2,4,5,7,10,12}; district geometry is at 398–412. |
| Does subject replication remove collinearity? | **No in these two allocations.** Every fitted dataset satisfies z=1−x, rank([1,x,z])=2, and normalized W x=1−x to numerical precision. |
| Noise and replication? | Lines 33–45 set N=10, 240 subjects, alpha=0.2, beta=1, SD=0.1. Lines 441–448 sample locations over the whole grid: district counts vary, with 20 expected per district. Lines 504–511 generate responses with independent subject errors. At rho=0, subjects share the district fixed mean plus their errors. |
| How is MSE aggregated? | Lines 590–599 compute mean squared beta error over the 10 replicates per allocation. `code/SimEstimateAnalysisMeans.Rmd:634–638` averages all 924 allocations for SRS and rows 313/612 for BSS; SRS includes the BSS allocations. |

For these checkerboards, the subject-level mean is still

    alpha + beta*x + psi*z = (alpha + psi) + (beta - psi)*x.

Adding subjects does not restore separate identification of beta and psi. The fit
reported an aliased spillover column. MSE near psi²=0.25 is an approximation; it is
not a claim that a ten-replicate result with estimated rho must equal 0.25 exactly.

The comparison scripts also construct subject-level W and fit individual outcomes:
`SpatialSim_NC_DOC.Rmd:430–471` (160 subjects, N=50) and
`SpatialSim_3x3.Rmd:473–511` (180 subjects, N=50). They were read, not rerun.

Two additional code observations were recorded without repairs or expanded runs:

- The 3×4 control-only exposure loop stops at source district 8 (`463–467`), whereas
  its commented both-arms loop uses all 12 (`470–473`). The 3×3 script similarly
  stops at 8 (`455–459`). This omission does not break z=1−x for either 3×4
  checkerboard, because every control still has a treated neighbor among 1–8.
- The DGP uses raw binary W; fitting row-standardizes W (`3x4:498,510–519`, with
  equivalent code in the other grids). This cannot explain the checked rho=0
  discrepancy, where the DGP's spatial term is zero. Nonzero-rho effects were not audited.

## Restricted rerun results

Only allocations 313/612 and five sampled allocations were fit, at control-only
spillover, rho=0, psi=0.5, alpha=0.2, beta=1, SD=0.1, N=10: **70 fits total**.
The original fitting, exposure construction, coefficient extraction and MSE code
were retained. The saved script's rho=0.01 was overridden to match the target row.

| Allocation | Rerun beta bias | Rerun beta MSE | Stored beta MSE |
|---|---:|---:|---:|
| 47 | −0.006968 | 0.000448092 | 0.000448092 |
| 102 | −0.007362 | 0.000211397 | 0.000211397 |
| 313 | −0.479063 | 0.241627745 | 0.241627756 |
| 505 | −0.001689 | 0.000350862 | 0.000350862 |
| 612 | −0.447825 | 0.209470567 | 0.209470551 |
| 628 | −0.002896 | 0.000886927 | 0.000886927 |
| 854 | −0.502212 | 0.252347770 | 0.252347771 |

- BSS average MSE **0.225549156**, bias **−0.463443586**; stored BSS MSE **0.225549154**.
- Maximum absolute MSE difference across the seven allocations: **1.60e-8**.
- Stored SRS mean over all 924 allocations: **0.03424694**, matching the manuscript's
  0.034. This was read from existing results, not obtained by rerunning all allocations.
- All 70 beta estimates were finite. Thirty aliasing warnings were retained: 20
  checkerboard fits and 10 fits of allocation 854, which also has z=1−x. The only
  other warning concerned the original Rmd's missing final newline.
- Both identities and model-matrix rank were asserted on all 20 checkerboard datasets;
  maximum numerical error in normalized W x=1−x was 2.22e-15.
- Five random IDs were selected with seed 20260927. The original seed=2024 stream
  was preserved by advancing the original coordinate/error draws for skipped
  allocations without fitting them. R 4.5.2; spatialreg 1.4-3; spdep 1.4-2; sf 1.1-0.

Stored evidence: `results/beta_mse/beta_mse_results_TrtNoSpill_3x4_A02_B1_P05_R00.csv`,
lines 314 and 613 (allocation rows 313 and 612).

Scratch artifacts remain at `/private/tmp/spatialcrt-ch2-audit-20260927/`:
`audit.R`, `executed_loop.R`, `comparison.csv`, `beta_estimates.csv`, `diagnostics.csv`,
`audit.rds`, `run.log`, `sessionInfo.txt`, and the original audit `README.md`.
They are temporary and may disappear; the findings and numeric evidence needed for
the manuscript decision are preserved in this note. No simulation code was added
to IncidenceDesign, and no Project 1 output was overwritten.

## Accepted PDF review, performed after the initial audit

The author identified these exact documents, both reviewed read-only:

1. Accepted submission:
   `/Users/ajwalther/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Lin)/Project 1 - Spillover Effects/Manuscript/Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2_Submission/Revisions_V2c.pdf`
2. Chapter 2:
   `/Users/ajwalther/GithubProjects/bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.pdf`

Reviewed the relevant methods, results and discussion text; visually inspected the
complete accepted pages 20–22 and 24, containing Figures 5–7 and Tables 1–3, and
the corresponding Chapter 2 PDF pages 25–27 and 29 (printed pages 23–25 and 27).
The relevant figures and table values agree between the two PDFs.

| Setting | Pattern presented in the accepted exhibits |
|---|---|
| 2×4 (Figure 5, Table 1) | BSS has a narrower spread across allocations, but often higher average MSE; it has lower average MSE for control-only spillover with rho=0.01. |
| 3×3 (Figure 6, Table 2) | BSS has lower reported average MSE across the presented settings. |
| 3×4 (Figure 7, Table 3) | BSS has lower reported average MSE except for both-arms spillover with rho=0.01. |

In particular, **both accepted Figure 7 and Table 3 show low control-only BSS MSE**.
The table reports 0.004 in a column labeled MSE ×10 (i.e., MSE=0.0004), with bias
−0.006 at rho=0, psi=0.5. This is not merely an isolated number in Chapter 2 prose.
The exact historical inputs to these exhibits have not been traced.

The paper's broader argument includes consistency across allocations and avoiding
poor allocations, not only minimizing average MSE. The two BSS layouts in the 2×4
and 3×4 grids are the **entire eligible set**, not two sampled examples; 3×3 has six
BSS allocations. Their restricted, often similar layouts help explain their narrow
across-allocation spread. Distinguish that spread from repeated-outcome uncertainty
within an allocation and from average MSE. A narrow spread alone does not establish
low error or a favorable worst-case value relative to another design.

This broader reading qualifies the interpretation of the audit; it does not resolve
the numerical discrepancy. Neither the final PDF nor the current code should be
silently substituted for the other as the verified source of the same exhibit.

## Author recollection and proposed cross-study framing (follow-up, 2026-09-27)

The author recalled Project 1's practical argument as BSS being a "good, not great,
but not bad" allocation: investigators could implement it reliably rather than
risk drawing an SRS allocation with poor MSE. The author wants Chapter 3 to explain
why its Checkerboard-versus-SRS findings do not simply discount Project 1's contribution.
This is a clarification of the intended scientific connection, not approval of final prose.

The defensible distinction is **average performance versus allocation risk**:

- For fixed simulation conditions and allocation a, let m(a) be the expected squared
  treatment-effect estimation error over repeated outcomes. SRS performance can be
  summarized by its mean over a, or by the upper tail/maximum of m(a) over allocations.
  A restricted design can have a higher mean but a lower upper tail than SRS; those
  two findings are mathematically compatible. This describes a possible trade-off,
  not a newly verified claim for every setting in either project.
- Project 1 enumerated the small-grid allocation sets and presented both mean MSE
  and the spread of allocation-specific MSEs. Its accepted narrative emphasized
  avoiding poor allocations, while acknowledging settings with higher BSS mean MSE.
- Project 2's primary comparisons average MSE over its simulation draws and conditions
  (with regime/configuration splits). Its SRS summary's block-level quantiles
  (`code/15_srs_benchmark.R:195–197`) are quantiles of already-averaged scenario MSEs,
  not the upper tail of conditional MSE over allocations at fixed conditions.
  Thus these summaries do not directly repeat Project 1's allocation-risk comparison.
- Project 2 also changes grid size, incidence heterogeneity, the unit of analysis,
  spillover specification and primary adjacency. Chapter 2 uses binary exposure to
  any treated rook neighbor; Chapter 3 uses a weighted treated-neighbor proportion
  (and a control-only modifier in that regime), with queen primary. Neither a BSS
  advantage nor the same ranking is guaranteed to transfer to this setting.

Do **not** claim that Checkerboard still protects against bad SRS allocations in
Project 2 without a comparison designed to assess that claim. A fixed allocation
removes allocation randomness; it does not guarantee low MSE, identification, or
protection against all uncertainty. Nor does the distinction above resolve the
separate accepted-exhibit reproduction discrepancy recorded earlier in this note.

Proposed discussion wording, not yet approved or inserted into the chapter:

> Chapter 2 emphasized the trade-off between average estimation error and sensitivity
> to the realized treatment allocation, motivating spatial restrictions as a way to
> avoid poorly performing allocations in the settings examined. This chapter extends
> the design comparison to larger grids with heterogeneous baseline incidence and
> evaluates average MSE relative to complete randomization. Checkerboard's higher
> average MSE here shows that its performance does not generalize uniformly to this
> setting; allocation consistency alone is insufficient to ensure accurate estimation.

No additional simulations or manuscript changes were authorized by this discussion.

## Author-confirmed reference rule (2026-09-27)

The author explicitly distinguished how Project 1 is referenced in the two outputs:

- **Dissertation:** direct references such as "In Chapter 2, we found..." or
  "Chapter 2 emphasized..." are appropriate because the earlier work appears in
  the same document. Prefer a precise statement of the finding and its conditions
  over an unqualified assertion that BSS is better.
- **CTJ manuscript:** the paper must stand alone. Use wording such as "In previous
  work, we found..." with an actual citation to the Project 1 BMC Medical Research
  Methodology paper. Do not use "Chapter 2", "Project 1", or assume the reader has
  the dissertation. The citation is appropriate for specific attributed findings
  even when the prose says "previous work"; that phrase alone is not a reference.
- Verify the BMC paper's final bibliographic details when preparing the CTJ citation;
  do not invent its DOI, year, volume or article number. No bibliography was changed
  in this discussion.

This writing rule is settled. The exact scientific wording, SRS framing and
Checkerboard placement remain under discussion; the proposed paragraphs above have
not become approved manuscript text merely because the reference rule is confirmed.

## Author-requested allocation-risk analysis (2026-09-27; implementation pending)

The author requested that Project 2 consider variation in MSE across allocations,
as in Project 1, for **all proposed designs and SRS**, before settling findings,
recommendations and manuscript prose. Checkerboard's poor average performance in
the larger heterogeneous-incidence setting should be assessed alongside downside
allocation risk. This is an analytical extension to plan, not a request to justify
Checkerboard regardless of the results. No implementation scope or simulation
budget has yet been approved.

### Feasibility check

The existing runner uses one outcome-noise vector for each allocation draw
(`code/05_run_simulation.R:289–310`). It summarizes 25 such draws within each
incidence surface (`200–237`) and saves scenario and surface summaries, not
allocation-specific repeated-outcome MSEs (`326–346`). The current oracle surface
file has 144,000 rows and no allocation ID or within-allocation error replicates.
Consequently, variability of those stored MSEs cannot isolate allocation risk.
Even retaining each fit's squared error alone would mix outcome noise with allocation
variation. A new nested simulation is needed for the requested distinction.

For a fixed incidence surface X and scenario, define

    m_d(a | X) = E_epsilon[(tau_hat - tau)^2 | design d, allocation a, X].

Estimate this by repeating outcomes while holding allocation a and X fixed, then
compare the distribution of m_d across draws from design d. Evaluate each design
according to its actual allocation probabilities; do not give rare and common
allocations equal weight by deduplicating draws without retaining frequencies.
There are too many balanced allocations on the 100-cluster grid for exhaustive
enumeration. A sampled maximum must not be called the true worst case.

### Proposed implementation plan for author approval

1. Specify a focused pilot with all eight candidate designs plus the SRS benchmark,
   oracle estimation, queen primary, tau=1 and separate spillover regimes/incidence
   configurations. Agree rho/gamma coverage and numbers of allocations and repeated
   outcomes before running. Use the existing model, design rules and key-seeding
   conventions; leave validated main-run files/checkpoints unchanged.
2. Repeat outcomes for each fixed surface/allocation and retain conditional MSE,
   bias, coverage, Monte Carlo uncertainty, aliasing and failures. Verify that only
   noise varies within an allocation and that aggregation recovers mean performance
   within Monte Carlo uncertainty. Do not conceal non-identification under rook.
3. Report mean allocation MSE, between-allocation variance (and SD for interpretation),
   an upper quantile (proposed 90th percentile), and mean MSE in the worst 10% of
   sampled allocations. Separate simulation error in estimated allocation MSEs from
   real allocation variation; increase replication if needed for stable tails.
   An observed maximum can be supplementary and explicitly sample-dependent.
4. Review precision and the mean-versus-tail trade-offs with the author before
   expanding to rook/tau sensitivities, revising recommendations or writing prose.
   The deterministic Checkerboard has no allocation variation conditional on X;
   this alone is not evidence of low MSE. Across-surface variation remains a separate
   question for every design.

The NC application should eventually use the same distinction. Current documented
state (`application/README.md`): the older synthetic-incidence service-area study
ran, real SUD incidence and 58-cluster weights are prepared, but the revised design
comparison on real incidence has not run. Do not describe its revised findings as
already established. Also, Project 2 changes more than incidence heterogeneity;
without a controlled comparison, do not attribute the changed Checkerboard ranking
to heterogeneity alone.

## Author clarification: SRS benchmark and recommendation (2026-09-27)

The author recalled Project 1's sequence: hold a treatment allocation fixed, sample
subjects onto the grid, assign treatment/exposure, generate responses and fit the
spatial-lag model, repeating to estimate that allocation's MSE. Comparing these
MSEs over the enumerated allocations gave the allocation-performance distribution.
This matches the inspected code. The repeated datasets varied both subject locations
(and therefore district sample sizes and subject-level W) and outcome noise; it was
not just repeated errors on an otherwise fixed subject dataset. A single fit gives
a squared error, while the allocation's estimated MSE averages squared errors over
its repeated datasets.

The author wants SRS to serve as an informative benchmark supporting a clear,
convincing recommendation among the systematic designs, rather than having readers
infer from average MSE alone that allocation strategy is inconsequential. The
proposed allocation-risk analysis addresses this concern for all designs, but its
pilot scope has **not** been approved by this clarification.

Keep the recommendation evidence-based and conditional:

- At queen/tau=1, the current IGSQ-to-SRS MSE ratio is about 0.63 under control-only
  spillover (about 37% lower average MSE), 1.08 under both-arms spillover (about 8%
  higher), and 0.71 when pooled over the studied regimes/settings (about 29% lower).
  These support a substantive control-only benefit, not universal superiority.
- SRS's good mean performance does not establish that every allocation is good;
  equally, it does not establish that its upper-tail risk is worse than a systematic
  design's. That is the open question for the proposed extension.
- Report SRS's actual strengths even if it remains a reasonable choice in some
  settings. A strong systematic-design recommendation should specify the conditions
  and criterion under which it improves on that benchmark, not discount the benchmark
  to secure a preferred conclusion.
- The unresolved Project 1 exhibit discrepancy remains separate; use its verified
  simulation structure to motivate the allocation-risk question without claiming
  that the previous numerical discrepancy has been reconciled.

No new simulation, recommendation, or manuscript wording was approved or implemented
in this exchange.

## Decisions still open for IncidenceDesign

1. Agree the scope of the new allocation-risk analysis above before settling final
   recommendations. Then agree restrained SRS benchmark prose, reported by spillover regime with queen
   primary and pooled results explicitly labeled. IGSQ's pooled advantage is not a
   universal advantage: at tau=1 under queen, RE is about 0.63 for control-only but
   1.08 for both-arms spillover. Isolation Buffer also improves on SRS under
   control-only spillover (RE about 0.74). Source: `results/srs_benchmark/srs_benchmark_summary.txt`.
2. Decide whether Checkerboard remains in the body or moves to an appendix. The
   audit verdict does not itself authorize dropping it. A five-design body would
   require updated rank comparisons, exhibits and prose.
3. Agree the one-sentence Chapter 2 link without claiming the unresolved mechanism
   explains the cross-study difference. Suggested wording, **not yet approved**:

   > Building on Chapter 2's comparison of block-stratified and simple random
   > allocation, this chapter evaluates a broader set of treatment assignment
   > designs under heterogeneous baseline incidence.

4. After author decisions, update/render/review the chapter, then derive the CTJ
   manuscript. CTJ must stand alone and cite the accepted BMC paper rather than
   refer to Chapter 2. The real-data application study remains separate pending work.

No decision to trace accepted-exhibit provenance further or contact the advisor has
been made. Do not start either action based only on this note.

## Ready-to-paste continuation prompt

Continue IncidenceDesign. Read AGENTS.md, docs/plans/TODO.md, the SRS/Chapter 2 handoff,
and docs/plans/ch2_check_findings_2026-09-27.md. The minimal audit found that current
code and CSV agree at BSS MSE 0.22555, but the accepted 3×4 exhibit reports 0.0004;
the exact source of the accepted exhibit remains untraced. Its figures and tables
also support a broader allocation-consistency argument that must be represented
fairly. Project 1 is accepted/final and stays read-only, as does bios-dissertation.
Discuss the SRS framing, Checkerboard placement and neutral Chapter 2 link with me
before changing the chapter. Do not infer that the audit has settled those decisions.
