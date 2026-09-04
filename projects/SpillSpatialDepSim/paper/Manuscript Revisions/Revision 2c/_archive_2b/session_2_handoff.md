# Session 2 Handoff — C1 (Estimands and causal framing)

**Date:** 2026-06-13
**Status:** C1 approved → LOCKED. Session ends here.

## What was done

**Core C1 (reviewer comment):**

1. **§2.3.1 β** — Redefined as the *direct structural coefficient*. Removed "total treatment effect" claim. TrtSpill note rewritten: β stays direct; spillover to treated units carried by ψ separately.
2. **§2.5 ρ-vs-ψ paragraph** — ψ: "causal, design-level quantity representing the indirect effect" → "structural, design-level parameter — the coefficient on z_i". `% TODO-C1` flag removed.
3. **§2.5 estimands paragraph (new, no heading)** — β/ψ/ρ are structural coefficients, not potential-outcome contrasts; new `Equation:ReducedForm` with (I−ρW)⁻¹ multiplier; E[Y(0,1)]−E[Y(0,0)] equals no single coefficient; ρ propagates even under control-only spillover; study recovers structural β.
4. **§5 Discussion Limitations (new scope note)** — Single permitted "causal" location. Goal is design guidance for structural β, not causal ID. SAR could support causal under extra assumptions; causal-centered estimands are future work.
5. **Intro + Limitations** — Three gratuitous "causal" uses scrubbed.

**Additional author-approved changes:**

6. **`(Rx.x)` tags purged** from all non-comment body lines. Rule: these tags must never appear in manuscript body.
7. **Equation consolidation** — `Equation:DataGen` removed as standalone equation; §3.1 cross-references `Equation:ReducedForm`. Simulation Procedure reference updated. Appendix C retains full equation + cross-reference.
8. **Table spacing** — `{\renewcommand{\arraystretch}{1.15}}` applied to all three result tables.
9. **σ² = 0.01 stated** — Added to §3.1 prose and Appendix C step. Source: `sd_val <- 0.1` in `04_run_simulation.R`.

## Key decisions

- **"causal" rule confirmed in practice**: grep → exactly 2 hits, both in the §5 scope note. Zero elsewhere.
- **Equation:ReducedForm is the single canonical DGP equation** throughout the manuscript.
- **σ² = 0.01** (σ = 0.1) confirmed from R code; now explicit in manuscript.
- **Reviewer-provenance tags** (`(Rx.x)` in body) are forbidden in the final manuscript. Comment-line tags (`% R1.x:`) are fine.

## Compilation state

46 pages, 0 errors, 0 undefined refs. Snapshot: `Diffs/snapshots/V2b_after_C1.tex`.

## Next session

Session 3 — C2 (Spillover definition and identifiability). Use **Opus 4.8**.

Start prompt:
```
Continue Revision 2b of the Spatial CRT manuscript. C3 and C1 are approved and locked. Read the plan at /Users/ajwalther/.claude/plans/i-d-like-to-revise-sharded-brook.md and memory at /Users/ajwalther/.claude/projects/-Users-ajwalther-Library-CloudStorage-OneDrive-UniversityofNorthCarolinaatChapelHill-UNC-Dissertation--Lin--Project-1---Spillover-Effects-Manuscript-Manuscript-Revisions/memory/MEMORY.md. Then begin the Brief for C2 (Spillover definition and identifiability).
```
