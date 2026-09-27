# IncidenceDesign open to-dos (ordered)

Maintained by each Claude session. Update it when an item finishes or changes, and before
the session ends print the remaining items and a ready-to-paste prompt for the next one.

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
