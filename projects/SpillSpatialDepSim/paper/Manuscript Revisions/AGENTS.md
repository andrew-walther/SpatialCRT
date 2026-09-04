# AGENTS.md — Project Context: Spatial CRT Manuscript Revisions

**Project:** Revisions to "Spatial Cluster Randomized Trials — Sampling Design with Spillover Effects & Spatial Dependence"
**Journal:** BMC Medical Research Methodology
**Author:** Andrew Walther (awalther@unc.edu), UNC Chapel Hill Biostatistics
**Status:** Revision round 1 — all changes implemented; response letter finalized at V3 as of 2026-04-14

---

## Directory Structure

```
Manuscript Revisions/
├── AGENTS.md                                  ← This file
├── Walther_SpatialCRT_LaTeX_Revisions/        ← SUBMISSION PACKET (zip this folder)
│   ├── main.tex                               ← ORIGINAL (read-only, do not modify)
│   ├── main.pdf
│   ├── Revisions_V1.tex                       ← REVISED manuscript (all 22 changes applied)
│   ├── Revisions_V1.pdf
│   ├── Revisions_Annotated_V1.tex             ← Annotated diff (red text + changebar markup)
│   ├── Revisions_Annotated_V1.pdf
│   ├── references.bib
│   └── figures/
│       ├── Applications/                      (NC DCS/MCO maps + BSS assignment maps)
│       ├── BetaMSEresults/                    (MSE result plots — 6 grid×scenario + 3 combined)
│       ├── InterventionSpilloverEffect/       (TrtSpill comparison figure)
│       ├── Neighbors/                         (Rook/Queen/2nd-order contiguity figure)
│       └── SpatialGridSetups/                 (Block/, Random/, combined grid PNGs)
├── journal comments/                          ← Raw reviewer/editor comment files (Word)
│   ├── editor_comments.docx
│   ├── reviewer_1_comments.docx
│   └── reviewer_2_comments.docx
└── revision response/
    ├── response_letter_V3.docx/.pdf           ← FINAL response letter for resubmission ✓
    ├── planning/                              ← Working documents used during revision
    │   ├── Revision_Plan_Item_by_Item.docx
    │   └── Supplement_Proposal.docx
    ├── response example/                      ← Style reference for response letter
    │   ├── response_letter_example.docx
    │   └── Ding-et-al-...-2025.pdf
    └── archive/                              ← Superseded drafts
        ├── Response_Letter.docx/.pdf          (V1 draft)
        ├── response_letter_V2.docx/.pdf       (V2 draft)
        └── Annotation_Reference_Guide.docx/.pdf
```

### Submission Packet
Two items to upload to BMC:
1. `Walther_SpatialCRT_LaTeX_Revisions/` — zip the entire folder (contains revised manuscript, annotated diff, figures)
2. `revision response/response_letter_V3.pdf` — standalone response letter to guide reviewers

---

## Key Scientific Concepts

### Study Design
- **Spatial CRT**: Cluster Randomized Trial where clusters are spatial units (e.g., farms, villages). Treatment assignment uses spatial sampling strategies.
- **SRS (Simple Random Sampling)**: All valid balanced allocations equally probable. 70 / 126 / 924 possible allocations for 2×4 / 3×3 / 3×4 grids.
- **BSS (Block Stratified Sampling)**: "Checkerboard" allocation — no two same-treatment clusters share a Rook edge. Only 2 / 6 / 2 valid allocations for 2×4 / 3×3 / 3×4 grids.

### Statistical Model
- **SAR (Spatial Autoregressive) model**: y = ρWy + Xβ + ε
  - ρ = spatial autocorrelation parameter (nuisance, background)
  - ψ = spillover parameter (causal, design-level)
  - These are **conceptually distinct**: ρ is a background property of the outcome; ψ is caused by the treatment spilling across cluster boundaries
  - Fitted via `lagsarlm()` in R's `spdep` package
- **Spillover mechanisms**:
  - **TrtNoSpill**: Spillover goes only from Intervention → Control clusters
  - **TrtSpill**: Spillover goes from any treated cluster boundary to neighbors (both directions)

### Key Results Summary

| Grid | Scenario | ρ | BSS MSE | SRS MSE | Winner |
|------|----------|---|---------|---------|--------|
| 2×4 | TrtNoSpill | 0.01 | 0.165 | 0.304 | BSS |
| 2×4 | TrtNoSpill | 0 | 0.262 | 0.173 | SRS |
| 2×4 | TrtSpill | 0 | 0.262 | 0.008 | SRS |
| 3×3 | TrtNoSpill | 0 | 0.001 | 0.092 | BSS |
| 3×3 | TrtSpill | 0.01 | 0.039 | 0.103 | BSS |
| 3×4 | TrtNoSpill | 0 | ~0.0004 | 0.034 | BSS |
| 3×4 | TrtSpill | 0.01 | ~0.223 | 0.112 | SRS ← checkerboard paradox |

**Checkerboard paradox**: In 3×4 TrtSpill with ρ=0.01, BSS's rigid structure amplifies spillover-induced bias (BSS MSE ~0.223 >> SRS MSE 0.112) because every intervention cluster is surrounded entirely by control neighbors, creating systematic bias unique to TrtSpill + background autocorrelation.

**Risk-management framing (key R1.5 revision)**: BSS provides near-zero variance across allocations (only 2–6 valid allocations) → bounds worst-case allocation risk. SRS sometimes achieves lower *average* MSE but with high variance across the 70–924 possible allocations → chance of good *or* bad draws.

---

## Revision Inventory (22 items implemented in Revisions_V1.tex)

### Supplement Moves (Appendix additions)
| ID | Content Moved | Appendix |
|----|--------------|----------|
| SUPP-1 | Spatial weights matrix equations (W, W_ij) | Appendix A |
| SUPP-2 | SAR log-likelihood formula | Appendix B |
| SUPP-3 | 8-step simulation algorithm enumeration | Appendix C |
| SUPP-4 | Moran's I and Geary's C formulas | Appendix D |
| SUPP-5 | Section 2.2.2 Distance-Based Neighbors | Appendix E |

### Editor Comment
| ID | Change |
|----|--------|
| E1 | Restructured Declarations under `\section*{Declarations}` with 7 BMC-standard subsections |

### Reviewer 1 (Major)
| ID | Change |
|----|--------|
| R1.1 | Added BSS novelty statement + McCann et al. citation paragraph to Introduction |
| R1.2 | Added ρ vs ψ distinction paragraph before Section 2.6 |
| R1.3 | Added inferential scope paragraph before Section 3.2 |
| R1.4 | Strengthened Introduction framing (same paragraph as R1.1/R2.2/R2.3) |
| R1.5 | Rewrote Sections 3.4, 3.5, 3.6 with risk-management framing (BSS bounds worst-case, SRS lottery) |

### Reviewer 1 (Minor)
| ID | Change |
|----|--------|
| R1.M1 | (Addressed via R1.2) |
| R1.M2 | Added boundary-level spillover limitation paragraph in Discussion |
| R1.M3 | Added inference constraint paragraph in Limitations |
| R1.M4 | (Addressed via R1.5 risk-management rewrite) |

### Reviewer 2 (Major)
| ID | Change |
|----|--------|
| R2.1 | (Addressed via SUPP moves) |
| R2.2 | Added McCann et al. paragraph to Introduction (same as R1.1/R1.4) |
| R2.3 | Same as above |
| R2.4 | Added best-effort BSS note to Section 4.4 |
| R2.5 | Added MSE = Bias² + Variance sentence to Section 3.3 |
| R2.6 | (Addressed via SUPP-3 algorithm move to Appendix C) |

### Reviewer 2 (Minor)
| ID | Change |
|----|--------|
| R2.M1 | "prevention" → "limitation" in Section 4 |
| R2.M2 | Added total effect note (ψ interpretation) |
| R2.M3 | Updated Figure 2 caption |
| R2.M4 | Added checkerboard paradox note to Section 5.1 Discussion |
| R2.M5 | Fixed Y_ik notation |
| R2.M6 | Fixed Z_i description |

---

## Code & Scripts

All scripts were session-local and are no longer needed — outputs are finalized and saved in the project directory. For reference, the approaches used were:

- **Revision script** (Python): Loaded `main.tex`, applied 22 string replacements via exact boundary extraction (`str.find()` + slicing) to handle whitespace sensitivity, wrote `Revisions_V1.tex`
- **Response letter script** (Node.js / `docx` package): Generated `response_letter_V3.docx` — B&W formatting, minor comments numbered R1.6–R1.9 / R2.7–R2.12
- **Annotated diff script** (Python / `difflib.SequenceMatcher`): Compared `main.tex` to `Revisions_V1.tex`; text-mode changes → `{\color{red}...}%`; math/table/algorithm blocks → `\cbstart...\cbend` (`changebar` package)

---

## Known Issues / Gotchas

1. **Whitespace in LaTeX replacements**: The source file has trailing spaces on some equation lines. When matching multi-line blocks, extract the exact substring from the loaded source rather than hard-coding expected whitespace.
2. **Moran's apostrophe**: In Python raw strings, `r"Moran\'s"` contains a literal backslash. Use a non-raw string or extract from source file directly.
3. **BSS allocation counts**: 2×4 → 2 allocations, 3×3 → 6 allocations, 3×4 → 2 allocations. SRS: 2×4 → 70, 3×3 → 126, 3×4 → 924.
4. **lagsarlm() in spdep**: Fits the SAR model; `listw` argument takes a spatial weights list object from `nb2listw()`.
5. **`\float@end` redefinition errors**: Two benign `LaTeX Error: Command \float@end already defined.` warnings appear in all compiles due to a conflict between the `float` and `algorithm` packages. The PDF compiles correctly despite these warnings.
6. **Annotated diff — math/table environments**: Lines inside equation/align/tabular/algorithm environments cannot be color-wrapped with `{\color{red}...}` because LaTeX alignment tokens (`&`, `\\`) must appear at the correct group depth. Use `changebar` (`\cbstart`/`\cbend`) for those blocks instead.
7. **`latexdiff` not available**: The standard `latexdiff` Perl tool cannot be installed in user mode in this environment. The custom `create_annotated.py` script replicates its functionality using `difflib.SequenceMatcher` + environment-tracking.

---

## R Packages Used in Paper
- `spdep`: spatial weights, `lagsarlm()`, Moran's I, Geary's C
- `spatialreg`: alternative SAR fitting
- Standard simulation utilities

---

## Next Steps (if additional revisions needed)
- Increment version to `Revisions_V2.tex` for subsequent revision rounds
- Compile revised `.tex` with pdflatex and verify no LaTeX errors before submission
- Archive current V3 response letter and increment to `response_letter_V4.docx`
