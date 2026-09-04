# Session 6 Handoff — C6 (Bias severity and buffer zones)

## Status
C6 implemented, pre-reviewed (subagent: all 6 checks PASS + 2 minor wording fixes applied), and approved by author.

## What was done

### Framing decisions (author-confirmed)
- Severity paragraph placed in **§5.2 Limitations**, not §5.1 Relevance — avoids competing with the positive comparative headline.
- Buffer zones framed as **infeasible on our grids** (every control borders a treated cluster → excluding removes whole control arm). Correctly positioned as future work on larger layouts, not a retrofittable fix.
- Concede 50–80% bias magnitude honestly, then reframe: it is almost entirely bias not variance, and it is a confounding/identifiability property shared by equally collinear SRS allocations — not a BSS-specific defect.
- §4.4 recommendation is now **conditional on anticipated spillover mechanism** (BSS when spillover flows into controls; SRS or buffer zones otherwise).

### Edits made to Revisions_V2b.tex

**Edit 1 — §5.2 Limitations (severity quantification)**
Appended to the existing collinearity paragraph: states bias ≈ −ψ = 50–80% of β=1.0; concedes severity; contextualizes with (a) error is almost entirely bias not variance, afflicts equally collinear SRS allocations too; (b) where identifiability holds, bias ~1–2% of β (Table 1 cross-ref). Closes by pointing to §5.3 via `\ref{section:FutureConsiderations}`.

**Edit 2 — §4.4 Recommendations (conditional framing)**
"BSS is the recommended approach" → recommended subject to explicit qualification tied to expected spillover mechanism. Mentions checkerboard paradox and buffer-zone alternative for bidirectional spillover + dependence regime.

**Edit 3 — §5.3 Future Considerations (buffer-zone paragraph)**
New paragraph: buffer-zone/"fried-egg" designs (cite jarvis_spatial_2017) restore identifiability by excluding contaminated boundary controls. Explicitly states infeasibility on our compact grids. Identifies comparative BSS/SRS/buffer-zone evaluation on larger grids as key next step. Added `\label{section:FutureConsiderations}` to subsection heading.

**Post-review wording fixes (subagent + author feedback):**
- "equals −ψ" → "approximately −ψ" (consistency with ≈ in derivation)
- "corresponding SRS allocations" → "equally collinear SRS allocations" (prevents misread)
- "which we discuss below under Future Considerations" → `which we discuss in Section~\ref{section:FutureConsiderations}`

### Snapshot
`Revision 2b/Diffs/snapshots/V2b_after_C6.tex`

### Pre-review result
All 6 checklist items PASS. Two IMPORTANT wording fixes applied before final compile. Compile: 48 pages, 0 errors, 0 undefined refs.

### Page count note
48 pages (was 46 after C4). The +2 pages are from the three new paragraphs in the double-spaced document format. The Session 7b length-reduction pass will address this; per the plan, length is not optimized until all science is locked.

## Scientific nuance clarified during session (for future reference)
- TrtNoSpill → BSS generally better, but the mechanism depends on the grid: on 3×3/3×4, the geometry breaks exact collinearity so BSS bias ≈ 0; on 2×4, collinearity is exact and BSS wins only via variance/worst-case bounding.
- TrtSpill → SRS generally better (2×4, 3×4), except 3×3 where BSS still wins. The "checkerboard paradox" is worst in 3×4 + TrtSpill + ρ.
- Buffer zones are not "optimal under TrtSpill" for our grids — they're infeasible on all three grids studied. They're the principled future direction for larger grids.

## Next session start prompt
```
Continue Revision 2b of the Spatial CRT manuscript. All six reviewer comments (C1–C6) are approved and locked. Read the plan at /Users/ajwalther/.claude/plans/i-d-like-to-revise-sharded-brook.md and memory at /Users/ajwalther/.claude/projects/-Users-ajwalther-Library-CloudStorage-OneDrive-UniversityofNorthCarolinaatChapelHill-UNC-Dissertation--Lin--Project-1---Spillover-Effects-Manuscript-Manuscript-Revisions/memory/MEMORY.md. Then: (1) implement CLEAN-1 through CLEAN-5 as a single pass; (2) run the length-reduction pass (7b) folding in the 2a short-version condensations; (3) run the final verification checklist and compile the final manuscript PDF.
```
