# Phase 2e: proposed curated SI outline (for review)

Author: Andrew Walther. Drafted 2026-11-25. Nothing in the SI source has been marked yet; this
is the outline to approve before any `\iffullsi` block goes in. Page counts are measured from
scratch builds of the current SI (after 2c and 2d) with the listed blocks removed. They are not
estimates, except where marked "est.".

## Page budget

| Version | Pages |
|---|---|
| Today's SI (curated = full until blocks are marked) | 65 |
| A. Curated, hiding only the blocks listed below | 43 (measured) |
| B. A plus four more exhibits to the full version (E1) | 40 (measured) |
| B plus two table merges (E8) | about 38.5 (est.) |
| B plus merges plus the 2g prose pass (repeated caveats stated once, 400-word paragraphs split and cut, figure descriptions moved to captions; the review's estimates) | about 30–32 (est.) |
| New content to add: τ subsection (4e) about 1.5, IGS sensitivity about 1.5 | +3 (est.) |

The ~30-page target is reachable only with the 2g prose pass. Hiding blocks alone leaves the
methods at about 15 pages and the case-study results at about 12, because most of the length is
prose. If the final count lands above 32, the next candidates are condensing the rook section to
one paragraph (−1), condensing Section 13 (−1.5) and moving the MCMC and prior settings of 5.6
to the full version (−0.5).

## Proposed curated table of contents

One numbered list, methods before results, no Parts. Exhibit names are today's; the curated
build renumbers them. "Keeps" lists the exhibits in the curated version; everything else in that
section moves to the full version. ★ marks a target the CTJ cites (it must stay curated).

Front matter: title, overview paragraph, Contents (about 1 page; List of Figures and List of
Tables in the full version only, E2).

**Methods**

1. Simulation setting (1.1 Spatial structure and outcome model; 1.2 Heterogeneous baseline
   incidence). No exhibits. About 1.5 pages.
2. Treatment assignment designs. Keeps Tables S1 (six designs) and S2 (variants). Full only:
   Figure S1 (example grid assignments; main Table 1 shows the same maps) and the
   Checkerboard-history paragraph.
   - 2.1 The two related variants.
   - **2.2 [slot] Variants of incidence guidance** (regions and saturation spread; Phase 3 IGS
     sensitivity methods, about 0.5 page).
3. Estimation and identification (receives CTJ cut T12, the `lagsarlm` agreement, already there).
   - 3.1 Identification diagnostics.
   - **3.2 [slot] Estimating the incidence coefficient after a trial** (Phase 4e, curated, Chapter 3
     mirror; methods and its short results together, about 1.5 pages).
4. ★ Simulation design, resampling and metrics (4.1 Statistical comparisons and reporting; 4.2
   Allocation-risk study). Keeps Table S3 (parameter grid), which moves up from 4.3. Full only:
   the rest of 4.3 (seed keys Table S4, checkpoints, metric and MC SE formulas).
5. Case-study methods: 5.1 Purpose and observed inputs (T5, heart-failure correction, already
   there); ★5.2 Geographic adaptations and budgets (T6b rounding, applied); 5.3 Outcome model and
   Monte Carlo evaluation; 5.4 Reproducibility and Monte Carlo precision (T2 split-half check,
   already there); ★5.5 Partial participation and sample size; 5.6 Sensitivity analyses.
   Keeps Figure S2 (service-area map), ★Table S5 (budgets) and Table S6 (regional orderings).
   Full only: Figures S3 and S4 (main Figure 1B and 1C show them) and Figure S5 (example NC
   assignments).

**Results**

6. Grid benchmark under queen contiguity: 6.1 Design performance (Table S7); 6.2 Rankings by
   incidence configuration (text only); 6.3 Spatial dependence, spillover magnitude and regime
   (★Figure S8); 6.4 The eight-design comparison (Table S12). Full only: Figures S6, S7, Tables
   S8, S9, S10, and (variant B) S11.
7. Formal statistical comparisons. Full only in variant B: Figure S9 (critical-difference
   diagrams; the text gives every separation).
8. Stability across treatment effect sizes. Keeps Figure S10. Full only: Figure S11 and its
   paragraph.
9. ★ Estimation error and its dependence on β (today's Section 6, moved here from the methods
   because it reports results, E3).
10. Rook contiguity (text). Full only: Table S13, Figures S12–S14 and their paragraph.
11. Non-oracle sensitivity analysis. Keeps Table S14.
12. Robustness, the SRS benchmark and allocation-specific risk: 12.1 Robustness across blocks
    (text); 12.2 Comparison with SRS (Figure S15); ★12.3 Allocation-specific risk. Full only:
    Table S15 and Table S16 (the SI copy of main Table 2).
13. **[slot] Sensitivity of incidence guidance to regions and saturation spread** (Phase 3
    results, one figure, about 1 page; placed before the case study, next to the spread-versus-gain
    material).
14. Case study: accuracy, inference and reach. Keeps ★Table S18 and Figure S20 (spread versus
    gain). ★Table S17 stays unless E8 folds it into S18. Full only: Figures S16 (eight-design copy
    of main Figure 4), S17, S18, S19.
15. ★ Case study: allocation downside (receives T7, already there). Keeps Tables S19 (downside
    ratios) and S20 (incidence-summary sensitivity). Full only: Figure S21 (= main Figure 3) and,
    in variant B, Figure S22.
16. Case study: participation, power and sensitivity: 16.1 Colleges needed and power (Table S21,
    Figure S23); ★16.2 Interval calibration (Tables S23, S24); 16.3 Guidance under partial
    participation (Table S25); ★16.4 Sensitivity analyses (Tables S26, S28). Full only: Table S22,
    Figure S24, Table S27 and 16.5 Implications for trial planning (main Discussion covers it).

References.

Curated exhibits (variant B, before merges): 6 figures (S2, S8, S10, S15, S20, S23) and 19
tables, plus the new figures from slots 2.2/13 and 3.2.

## Decisions (reply by ID)

| ID | Decision | Recommendation and reason |
|---|---|---|
| E1 | Variant A or B | **B**: moves Figure S9, Figure S22, Table S11 and Table S17 to the full version (−3 pages). None is cited by the CTJ except S17 (see E8); the text states their results. |
| E2 | List of Figures and List of Tables in the curated SI | **Full version only**: they take 6 pages because they print whole captions; the curated Contents is enough at about 25 exhibits. |
| E3 | Section 6 (error decomposition) moves into the results, after stability across β | **Move**: it reports results, and methods-before-results is the plan. The CTJ's pointer follows the label, so it prints the new number automatically. |
| E4 | How the curated text points to full-only material | **Wrap each pointer sentence or clause in `\iffullsi`**, so the curated build prints no "??", and add one overview sentence that a full-length version with additional tables and figures is available. Every wrapped sentence comes to you as a sign-off list before it is applied (about 25 pointers). |
| E5 | Label the full version | **Title "SI Appendix (full-length version)"** with a first-page note that it is not the submission SI and that its exhibit numbers differ from the submitted SI. |
| E6 | τ subsection slot (Phase 4e) | **3.2, inside Estimation and identification**, with its short results in the same subsection, as Phase 4e describes ("SI Methods subsection"). Alternative: a results section after 9. |
| E7 | IGS sensitivity slots (Phase 3, still running) | **Methods 2.2 and results Section 13**, as the SI review proposed, next to the spread-versus-gain material. |
| E8 | Two table merges (implemented later with a sign-off list) | **Fold Table S17 into S18** and repoint the CTJ's "SI Table S17" to S18 (a CTJ edit, R10 from 2b), and **merge Tables S7 and S12** into one nine-strategy performance table. About −1.5 pages. |
| E9 | Moved CTJ text (2f) | **No new slots needed**: every Phase 1 receiver (3, 5.1, 5.2, 5.4, 15) stays curated, and T6b is already applied. 2f checks each against `trimmed-material.md`. |

## Checks for the implementation step

- Every `\ref{SI-…}` in the CTJ must point to a curated label. The verifier currently checks
  that the SI source defines each one; it should also fail if a cited label sits inside an
  `\iffullsi` block (added with the markings).
- The curated build must log no undefined references; the full build must stay identical to
  today's text apart from E3's move and E5's title.
- Before a block goes full-only, confirm it exists in Chapter 3 Appendix B or in the full SI
  (the full SI keeps everything, so nothing is lost).
