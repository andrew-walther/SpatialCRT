# Session 3 Handoff — C2 (Spillover Definition and Identifiability)

**Date:** 2026-06-13
**Status:** C2 approved and locked.

## What was done

Two edits to `Revisions_V2b.tex`:

**§2.3.2 SpilloverEffectMethods:**
- Replaced vague "There are many ways to model a spillover effect" paragraph with a clear two-paragraph treatment
- First paragraph: acknowledges spillover can generally be a distance-decay or count-of-neighbors function, but in this study it is binary cluster-level adjacency; explains why ψ is a scalar (z_i ∈ {0,1}); explains why spillover is homogeneous within clusters (z_i is cluster-level); forward-refs §5 for richer alternatives
- Second paragraph: introduces identifiability acknowledgment — z_i can be exactly collinear with x_i for some allocations, making β and ψ non-identifiable from the model fit alone; defers the full derivation to §3 (C4's job)

**§3.2 Simulation Procedure (SampledPoints paragraph):**
- Added clarification that individual within-cluster coordinates are used only for cluster membership assignment and visualization
- States the fitted model depends on each subject only through its cluster; within-cluster variation in y_ik contributes only to error-variance estimation, not to structural parameters β, ψ, ρ
- Notes individual-location modeling as a future extension

## Key decisions

- **"nearest intervention cluster"** used (not "nearest intervention source" — reviewer's exact complaint; "source" was undefined)
- **"from the model fit alone"** qualifier retained in identifiability sentence — precise and allows C4 to elaborate without contradiction
- **"discussed in Section 5"** without "(Future Considerations)" — author preference
- **Pre-review subagent** ran twice; final verdict CLEAN (two advisory notes both confirmed non-blocking: y_ik already defined in §2.5 before §3.2; "from the model fit alone" is intentionally precise)

## Files

- `Revision 2b/Walther_SpatialCRT_LaTeX_Revisions_V2b/Revisions_V2b.tex` — updated
- `Revision 2b/Diffs/snapshots/V2b_after_C2.tex` — snapshot of final C2 state
- `Revision 2b/Diffs/pdfs/C2_diff.pdf` — per-item diff (C1 snapshot → C2 snapshot), compiled clean

## Next

**Session 4 — C4 (Estimation model + bias derivation).** Must run the C4 verification gate first: locate the R simulation code and Table 1 source to confirm whether z was included in the fitted model or omitted/dropped, and what lagsarlm() actually returned for ψ̂ under BSS/TrtNoSpill.
