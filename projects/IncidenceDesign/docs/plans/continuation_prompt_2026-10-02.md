# Integrated continuation prompt — IncidenceDesign, 2026-10-02

This updates the September 27 prompt with the October 2 application-study request.
The eleven September-session commits were successfully pushed through `e4599ca`;
the earlier automatic approval rejection was resolved by explicit author approval.
That permission applied to those commits, not every future push.

Copy the following into Codex or Claude with the IncidenceDesign project open:

> Continue Project 2 in
> `/Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign`.
> Work directly in the checkout. Read the root and project AGENTS.md, ROADMAP.md,
> docs/plans/TODO.md, docs/plans/ch2_check_findings_2026-09-27.md,
> docs/plans/allocation_risk_findings_2026-09-27.md, and
> docs/plans/nc_sud_application_plan_2026-10-02.md. Read application/AGENTS.md
> and README before working there. Verify current Git and project state;
> historical progress paragraphs may be superseded by later notes.
>
> Our goal is a clear, defensible treatment-design recommendation supported by
> the revised grid study, the allocation-risk findings, and a completed application
> on the 58 NC community-college service-area clusters using 2018–2021 SUD incidence.
> Use observed NC incidence as fixed input; generate no synthetic incidence for
> the primary application. Simulate treatment allocations and outcomes to assess
> estimation error under known parameters.
> The allocation-risk pilot and R=400 refinement are complete. They support
> saturation designs under control-only spillover; SRS remains competitive/better
> under both-arms spillover. They use selected settings and do not establish exact
> true tails or a unique incremental benefit from incidence guidance. Keep their
> numbers distinct from the full-study results.
>
> First, interview me on the outstanding application decisions in its design
> document. The author chose the existing continuous SAR approach on October 2,
> using observed SUD incidence to inform allocation, with τ on the simulation
> scale. Cluster ranks for allocation can be computed directly from observed
> rates: transformation is not required for rank-based assignment. Separately
> agree the SAR baseline covariate scale and regional incidence summary (the
> current saturation code averages cluster ranks). Keep the study focused. Agree which
> year's incidence informs allocation, treatment-budget/design adaptations,
> fixed geographic partition and any historical-incidence sensitivity. Present
> the concrete implementation plan and obtain the required approval before code.
> Continue independent input checks and documentation while answers are pending.
>
> Also settle SRS framing by spillover regime, Checkerboard body/appendix placement,
> and the neutral Chapter 2 link before substantive manuscript edits.
> Keep the author's Project 1 allocation-consistency argument central: BSS offered
> reasonably good accuracy and limited poor-allocation downside in its studied
> settings. Project 2 must compare both mean accuracy and allocation risk with SRS;
> neither SRS's pooled mean nor the earlier BSS recommendation settles that question.
> The application complements this unfinished alignment work.
> Preserve Project 1 as accepted/final and read-only; do not reopen the audit or claim its
> historical numerical discrepancy has been reconciled. The shared question is
> whether allocation restrictions improve accuracy and limit poor-allocation risk;
> design rankings need not transport unchanged across the two settings.
>
> After approval, implement and verify the revised real-SUD application using
> named cached weights, current design rules, key-seeded draws, warning/alias and
> failure reporting, manifest-checked outputs, and the validated ML estimator.
> Include SRS, report all four yearly settings by spillover regime, and show
> treatment counts/population shares as well as estimation performance. Avoid
> replacing older synthetic results until the revised run passes checks. Keep
> primary results distinct from any agreed earlier-year-incidence sensitivity.
> Use targeted replication for allocation-risk precision rather than an unnecessary
> broad parameter sweep. Do not assume the preferred design must win.
>
> Use the completed methods/results to revise the application's chapter section,
> with detailed adaptations, yearly results and sensitivities in its appendix.
> Render and review in fresh passes, then derive the standalone CTJ manuscript
> and supplementary material. Dissertation prose can say Chapter 2; CTJ must
> cite the verified BMC Medical Research Methodology paper and stand alone.
> Edit Chapter 3 only in IncidenceDesign and use the approved sync hook for the
> bios-dissertation copy. Do not modify Project 1.
>
> Keep restricted county-level source and derived data ignored. Authorized
> cluster-level results and aggregate statistics may be tracked. Save key
> findings in docs/plans/, progress in ROADMAP, and update README/AGENTS/TODO.
> Commit logical steps without AI co-author lines. Never push without asking.
