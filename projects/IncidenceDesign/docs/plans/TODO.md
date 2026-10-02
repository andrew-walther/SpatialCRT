# IncidenceDesign open to-dos (ordered)

Maintained by each Claude session. Update it when an item finishes or changes, and before
the session ends print the remaining items and a ready-to-paste prompt for the next one.

**Active continuation (2026-10-02):** the integrated prompt is
[continuation_prompt_2026-10-02.md](continuation_prompt_2026-10-02.md).
The author requested a short application design document and interview before
implementation: [NC SUD application plan](nc_sud_application_plan_2026-10-02.md).
Continuous SAR is approved. In the author's follow-up the primary outcome is
education: observed SUD incidence guides allocation but does not enter the
simulated outcome baseline (β = 0); τ is constant within a scenario. An explicitly
hypothetical β = 1 baseline sensitivity uses rank-scaled X and an X-adjusted fit;
the primary fit uses treatment and spillover without X. The proposed queen grid,
focused sensitivities and precision targets are accepted. The author-requested
[independent alignment review](nc_sud_independent_review_2026-10-02.md) is complete:
no scientific scope change needed. The [implementation plan](nc_sud_implementation_plan_2026-10-02.md)
addresses singleton precision, frozen geography, correct model matrices,
identification, duplicate covariance and cross-year primary reuse; awaiting approval.
Regional mean cluster incidence ranks are approved
for primary incidence-guided saturation: per-capita rates → cluster ranks → mean
regional rank → 80/60/40/20% saturation. The author requested consideration of
mean-rate and population-weighted regional-rate sensitivity; a focused proposal
is accepted in the application scope. One fixed four-region
partition across years is approved, using the existing geographic/balance criterion
and pooled person-years divided by four as reference population. The author also
approved existing design adaptations and budget differences, year-specific spatial
blocks frozen within each year, and Balanced Halves treating 14+15 (randomly
choosing the half with 15), exactly 29 overall. Next: implementation-plan approval,
code/tests, verified production, chapter/appendix and CTJ/SI. Rank-based allocation can use ranks of observed rates
directly and does not require rate normalization.
The eleven September 27 commits were pushed; earlier "not pushed" statements in
conversation history are superseded. New commits still require push authorization.

**Completion goal reaffirmed October 2:** finish the NC application; integrate
its findings into a full thesis chapter and appendix; then trim/reorganize that
long-form source into full CTJ manuscript and supplementary-material drafts.
CTJ remains a required deliverable. The latest author instruction is to disregard
the recorded deadline dates and finish the application/chapter promptly, then
derive CTJ/SI for submission soon afterward. Planning and pilot runs are intermediate milestones.

**October 2 interview decisions:** retain the completed allocation-risk grid pilot
as supporting evidence and use the NC application as the next confirmation;
no extension to all ten grid incidence surfaces is needed first. The Chapter 2
linking sentence is approved and saved in the application plan. The author does
not require a direct Project 1/2 results comparison in the body; a brief reference
may explain the extension and inclusion of SRS/BSS alongside other designs.
Checkerboard remains in Project 2's main design comparison, with detailed
identification diagnostics in the appendix (author-approved follow-up). CTJ must
stand alone and cite the verified Project 1 publication as the work being extended.
Historical statements below that the linking sentence, Checkerboard placement or
grid-extension decision are pending are superseded by this update. Application
implementation decisions remain open.

1. **Simple Random Sampling (SRS) benchmark.** Add it to the simulation, regenerate the
   results, and update the chapter. Added by the author 2026-09-25.
   Plan: `~/.claude/plans/pasted-content-id-a3af-you-are-synchronous-seal.md`.

   **Status (2026-09-27): segment A (steps 1–3) done; segment B (steps 4–5: chapter, wrap-up) next.**
   - **Read before continuing:** [Chapter 2 check and accepted-PDF review](ch2_check_findings_2026-09-27.md).
     The limited reproduction discrepancy remains unresolved; Project 1 is finalized
     and stays read-only. SRS framing, Checkerboard placement and the Chapter 2 link
     still require author decisions. The accepted exhibits' allocation-consistency
     argument must be considered alongside average MSE.
   - **Author-requested extension (2026-09-27):** evaluate variation and upper-tail
     MSE across allocations for all candidate designs and SRS before finalizing
     recommendations. Existing surface/scenario summaries cannot isolate this:
     repeated outcomes per fixed allocation are needed. **Subsequently authorized and
     completed as a nested pilot**, with targeted R=400 precision follow-up; read
     [allocation-risk findings](allocation_risk_findings_2026-09-27.md) before writing.
     Queen/oracle/tau=1, two surfaces/configuration and parameter corners only.
     IGSQ retains a control-only mean and estimated-tail benefit; SRS remains better
     under both-arms spillover. Finite-noise tail uncertainty and the plain-saturation
     results preclude claims of universal or uniquely incidence-guided superiority.
     Further analysis priorities are saved in the findings note ("Further analysis
     recommendations"); they are proposals for discussion, not approved new runs.
   - Commits: 0052161 (code, tests), 197aa62 (full run), 3366af2 (summaries, figures).
   - SRS is Design 9: complete randomization, exactly 50 of 100 treated. It is a benchmark,
     not a proposed design: shown as gray reference lines in figures, as a set-off row
     plus MSE/MSE_SRS in tables, and kept out of the six-design rank tests.
   - Run: `results/sim_data/*_20260927_000502.rds`, 14,400 scenarios per estimator.
     `results/srs_reproduction_check.txt`: all Design 1–8 rows are identical to the
     2026-09-24 run. The 12/13/14 txt/rds outputs are byte-identical, so no existing
     chapter number changes.
   - SRS numbers: `results/srs_benchmark/srs_benchmark_summary.txt` (from
     `code/15_srs_benchmark.R`). Updated figures are in `results/six_design_manuscript/`;
     none are copied to the chapter yet.
   - Stop rule not triggered. Oracle pooled MSE of IGSQ vs SRS: queen 0.091 vs 0.128
     (RE 0.71) at τ = 1, RE 0.73 pooled over τ; rook RE 0.84 / 0.87.
   - **AUTHOR DECISION NEEDED before writing the chapter: SRS beats most proposed designs.**
     - Queen, τ = 1, RE = MSE/MSE_SRS: Balanced Quartiles 1.02 (paired test vs SRS
       p = 0.31, not different); Isolation Buffer 1.25; High Incidence Focus 1.91;
       2x2 Blocking 2.49; Checkerboard 8.46.
     - Rook: Balanced Quartiles 0.98, then 2x2 Blocking 1.51, Isolation Buffer 1.56,
       High Incidence Focus 2.43, Checkerboard 5.64.
     - Only IGSQ (and the consolidated Saturation Quadrants) beat SRS significantly
       (paired p < 1e-8, better in all 50 units).
     - Under both-arms spillover (queen), SRS beats IGSQ: RE 1.08. Under control-only,
       IGSQ RE is 0.63.
     - This changes the chapter's framing (value add only for IGSQ), so agree it before
       the writer runs.
2. **Application study on the real SUD data.** Apply the designs, including SRS, to the yearly
   cluster-level incidence surfaces, then fill the chapter's Application results stub in place.
   The data are already ingested and aggregated (`application/README.md`, `AGENTS.md`).
   Design → implement/verify → run/interpret → revise chapter/appendix → CTJ/SI.
   Current design/interview document: `nc_sud_application_plan_2026-10-02.md`.
3. **CTJ manuscript (manuscript-plan step 5).** SRS note (author, 2026-09-26): the CTJ stands
   alone, with no "Chapter 2" references. Link SRS to Project 1 only by citing that paper
   (accepted at BMC Medical Research Methodology), e.g. "we compare performance to SRS,
   as in [Paper 1]". Derive it from the finished chapter: body about
   3,000–3,300 words, at most 6 exhibits, an unstructured abstract of at most 250 words (check
   SAGE's live guidelines), and cut material saved to supplementary files. Can start before
   item 2 finishes, with the application written as a proposal.
   - **Reference rule reaffirmed by author, 2026-09-27:** dissertation references can
     say "In Chapter 2..."; CTJ references must say, for example, "In previous work..."
     and cite the BMC paper. Verify its final bibliographic details at manuscript
     preparation. See the author-confirmed rule in `ch2_check_findings_2026-09-27.md`.
4. **Chapter 2 (Project 1) audit.** Minimal authorized code check and accepted-PDF review
   completed 2026-09-27; further work paused for discussion. Chapter 2 is accepted and
   final: no reformulation or critique of it in Chapter 3. Current findings and limits:
   [ch2_check_findings_2026-09-27.md](ch2_check_findings_2026-09-27.md). Earlier draft-only
   observations in `paper/dissertation_chapter/notes/a5_technical_notes.md` (Q3) must be
   read with that follow-up. No further provenance tracing or advisor contact authorized.

Smaller open items:
- 17 inline `NEEDS-AUTHOR-CONFIRMATION` comments in the chapter's Application section.
- The Q&A backup slides plus a companion document for the full prelim (all 3 projects),
  seeded by `notes/a5_technical_notes.md`.
- Commit the literature-review PDF in bios-dissertation: the author's call.
