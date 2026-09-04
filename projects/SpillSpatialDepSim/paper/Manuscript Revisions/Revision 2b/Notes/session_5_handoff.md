# Session 5 Handoff — C5 (Frequentist inference under BSS)

## Status
C5 implemented, pre-reviewed (subagent: all 6 checks PASS), and approved by author.

## What was done

### Framing decision (author-confirmed)
- "Minimax" as a formal label was considered and rejected — we have not proven a worst-case optimality bound, only observed it empirically across the enumerated allocations. "Bounds worst-case estimation error across the feasible allocations" is the correct, defensible phrasing.

### Edits made to Revisions_V2b.tex

**Edit 1 — §3 inferential-scope paragraph (previously labelled % R1.3)**
Replaced ~4-sentence "design-stage tool" paragraph with a 3-paragraph block (% R1.3 / C5: Inferential scope):
- Para 1: Repositions BSS as the deliberate allocation mechanism (investigator selects a BSS-valid allocation); concedes no usable randomization distribution (2–6 allocations); reframes near-determinism as deliberate risk-averse design choice that bounds worst-case error across feasible allocations.
- Para 2: Distinguishes design-based from model-based inference; states SAR MLE conditions on the realized fixed allocation and derives β̂ sampling distribution from the stochastic outcome process — not impaired by small BSS allocation set ("does not invalidate it"). Clarifies simulation variation is outcome-level Monte Carlo, not allocation-randomization variation.
- Para 3: States the relevant question is whether the principled fixed allocation yields a low-error estimator; exhaustive enumeration answers that directly.

**Edit 2 — §5.2 Limitations "Inference constraint" paragraph (previously % R1.3)**
Expanded the 2-sentence paragraph to 4 sentences (% R1.3 / C5: Inference constraint):
- Retains: findings characterize design-space performance, not a sampling distribution from a single randomization.
- Added: explicit concession that BSS does not generate a randomization distribution and does not support permutation/randomization-based inference.
- Added: model-based SAR MLE is valid regardless (conditions on fixed allocation, derives from random error term).
- Closing: BSS near-determinism is a deliberate risk-averse design choice, not a limitation of inference.

**Wording fix (pre-review suggestion)**
"does not impair it" → "does not invalidate it" — avoids brushing against the separately acknowledged collinearity/efficiency caveats (C4).

### Snapshot
`Revision 2b/Diffs/snapshots/V2b_after_C5.tex`

### Pre-review result
All 6 checklist items PASS. One optional softening applied ("does not invalidate it"). No critical or important issues.

## Next session start prompt
```
Continue Revision 2b of the Spatial CRT manuscript. C3, C1, C2, C4, and C5 are approved and locked. Read the plan at /Users/ajwalther/.claude/plans/i-d-like-to-revise-sharded-brook.md and memory at /Users/ajwalther/.claude/projects/-Users-ajwalther-Library-CloudStorage-OneDrive-UniversityofNorthCarolinaatChapelHill-UNC-Dissertation--Lin--Project-1---Spillover-Effects-Manuscript-Manuscript-Revisions/memory/MEMORY.md. Then begin the Brief for C6 (Bias severity and buffer zones).
```
