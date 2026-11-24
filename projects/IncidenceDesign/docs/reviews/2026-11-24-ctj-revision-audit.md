# Audit of the first CTJ revision round (applied 2026-11-14/15)

Written 2026-11-24. Author: Andrew Walther. Read-only audit; no manuscript, code or build was changed.

Sources checked against: the FCL register (`paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/feedback.json`),
`proposed-edits.md` (the approved text and Andrew's edits, which govern), the 2026-11-14 meeting notes,
and the triage in `docs/plans/2026-11-14-chapter3-review-session-plan.md`. Applied text read in
`CTJ_Manuscript.tex` (C) and `CTJ_Exhibits.tex` (E); carry-over read in `Thesis_Chapter.qmd` (Ch3) and
the SI source. Lines are as of HEAD `4386bcb`.

## Summary

**Freshness.** `CTJ_Manuscript.tex` was last edited 2026-11-15 06:39; `CTJ_Manuscript.pdf` and `CTJ_Reading.pdf`
were built 06:40, and `CTJ_Exhibits.tex` (03:26) is older than both PDFs. The PDFs are current. `CTJ_Reading_Markup_V6.pdf`
was added in the closing commit. The reading layout shows the equation on one line.

**Notation.** No γ and no τ on **Z** anywhere in the CTJ source, exhibits, SI source or Chapter 3 prose. The only leftovers
are in `paper/ctj_manuscript/figures/Figure_Legends.md` (see item N-1).

**Counts by status** (60 FCL rows and 37 meeting rows, 97 in all: D-1 to D-5, M-01 to M-13, S-01 to S-04,
A-01 to A-03, I-01 to I-06, W-01 to W-04, Q-01 to Q-02). Row-level statuses in the table are authoritative.

| Status | FCL | Meeting | Total |
|---|---:|---:|---:|
| done (including author edits, and FCL-043 with the move declined) | 55 | 26 | 81 |
| partial | 2 | 8 | 10 |
| missing | 0 | 0 | 0 |
| deviates from disposition | 0 | 0 | 0 |
| rejected by author (fine) | 1 | 0 | 1 |
| no change by plan (reply only) | 2 | 0 | 2 |
| not auditable (style sweep: W-01, W-02) | 0 | 2 | 2 |
| done in Chapter 3 only, CTJ optional (D-5) | 0 | 1 | 1 |
| total | 60 | 37 | 97 |

Partial rows: FCL-024, FCL-028, D-3, D-4, M-08, M-10, S-03, A-01, W-03, Q-01.

No item deviates from an approved disposition in a way that changes meaning. Wording drift is small and listed under
"Wording drift".

### Items needing attention (partial, or outside the original lists)

1. **N-1 Stale figure legends.** `paper/ctj_manuscript/figures/Figure_Legends.md` still has the old notation
   (`$\tau=1$` for the treatment effect, `$\gamma$`, "outcome-dependence", "Application mean MSE"), so Figures 2 to 4
   disagree with the exhibits. This is why the verifier fails (below). Not in any sign-off list.
2. **N-2 Both checkers fail on stale references, not on wrong numbers.** `verify_manuscripts.py` stops at the legend check, and with
   that bypassed it stops again at a label (`tab:nc-summary-sens`) that the 2026-11-15 appendix reorganization removed.
   `test_manuscript_exhibits.R` still searches the CTJ text for `$\tau=2$` and `$\tau=3$`, which are now `$\beta=2$`
   and `$\beta=3$` (the values it checks are in C 206, 0.160 and 0.160, then 0.160 and 0.182, and look right). The reconciliation has not run on the
   current text since the notation swap.
3. **M-08 / FCL-028 SI.** The SI still says "identified" for the Checkerboard result (SI 315, 317, 334, 346, 2168: "not identified",
   "treatment is identified", "β is not identified"). The plan lists S as done for T-22.
4. **D-4 / M-10 / FCL-024.** The t quantile (Leyrat) and the SE factor (MacKinnon and White) are cited. The normal (Wald)
   interval is described but not cited. The approved text did the same, so this is a gap in the approved text.
5. **D-3 / A-01.** `paper/ctj_manuscript/cover_letter/Cover_Letter_Draft.md` lines 14 and 16 still say "application". The SI
   section label and `submission_checklist.md` line 37 also still say "application" (labels and notes only).
6. **W-03.** C still has 17 ", so" or "; so" joins and E has 5 (for example C 49, 141, 147, 158, 162, 165, 185, 187). Ch3 has about one.
   Only the new text and M-16 were swept in CTJ.
7. **S-03 (CTJ).** Plan marks it done for Ch3 only. C 197 and the Figure 2 caption state that configuration changes MSE
   by under 10%, but nothing in C justifies five configurations or says none changed the ranking.
8. **Q-01 (CTJ).** Plan lists C but marks it done for Ch3 and SI only. C does not mention the earlier study's binary spillover.
9. **D-5 (CTJ).** Ch3 Future Work has it (Ch3 162, 636); C has no sentence (plan said optional).
10. **Minor consistency.** C 243 heading says "SUD education pilot", while C 49 and 147 say "pilot trial". The plan notes this as open.

## Wording drift from the approved text

- **C 158 (M-06).** Adds "as the HC1 covariance estimator does" and "although in the spatial lag model this rescaling is approximate
  and was chosen by calibration (SI Section 15.2)". It adds a caveat and does not contradict the decision. It matches Ch3 148.
- **C 151 (M-01).** Adds "an exposure mapping" with two citations, from round 5. Harmless.
- **C 147 (I-06).** "pilot study" became "pilot trial" in commit `324cd80`, which resolves the open item in the sign-off list.
- **C 212 (R-07, FCL-060).** Heading follows Andrew's edit and uses neither "best" nor "preferred" (his instruction).
- **C 215 (R-08, I-10).** Follows Andrew's edit; the abstract (C 53) says "Across the six designs and SRS". Numbers unchanged by decision.
- **C 234 (R-09).** Matches the approved paragraph.

Nothing found that reverses a decision.

## Item table

Status key: done, partial, rejected (author), no change (reply only). "Ch3" column says whether the item reached Chapter 3 when it should have.
Quotes are from the CTJ source.

### FCL annotations

| ID | Asked | Disposition | Landed (CTJ) | Ch3 | Status |
|---|---|---|---|---|---|
| FCL-001 | Citation at end of sentence | Accept (I-03) | C 143 "serve different objectives.\cite{...}" | n/a (Ch3 49-53 own wording) | done |
| FCL-002 | Who is "its"? | Accept (I-04) | C 143 rewrite, no "its" | n/a, sentence not in Ch3 | done |
| FCL-003 | "checkerboard design" | Accept (I-04) | C 143 "in a checkerboard design with neighbors defined by shared edges" | n/a | done |
| FCL-004 | Cite first paper | Accept (I-04) | C 143 ends `\cite{walther_spatialcrt_2026}` | Ch3 55 names Chapter 2 | done |
| FCL-005 | What does the sentence mean | Accept (I-05) | C 145 "In that study, we judged each design by its average estimation accuracy" | Ch3 55 | done |
| FCL-006 | Delete sentence, revise first sentence | Accept (I-05) | C 145 first sentence now covers SRS and Checkerboard | Ch3 55 | done |
| FCL-007 | Sentence belongs last | Accept (I-05) | C 145 "In this paper, we compare six designs and SRS" is last but one | Ch3 59-61 | done |
| FCL-008 | Cite for "Two allocations..." | Accept (I-02) | C 141 `\cite{baird_optimal_2018,leung_rateoptimal_2022,walther_spatialcrt_2026}` | Ch3 49 | done |
| FCL-009 | Lin's design-stage sentence | Accept (I-03) | C 143 Lin's wording, "independent-set and spatial cluster designs" | n/a | done |
| FCL-010 | Unclear checkerboard sentence | Accept (I-04) | C 143 "every neighbor of a treated cluster is a control" | Ch3 150 has the content | done |
| FCL-011 | Is this checkerboard? | Accept (I-05) | C 145 "block stratified sampling... we call this design Checkerboard here" | Ch3 61 | done |
| FCL-012 | "In this article, we aim to plan" | Edited by author (I-06) | C 147 "we also use a case study, a simulation based on observed incidence" | Ch3 59, 374 | done (author edit) |
| FCL-013 | "further challenges" | Accept (I-01) as "violates" | C 141 "Spatial interference violates the assumption" | n/a | done |
| FCL-014 | Understand allocation risk | Accept (I-07, M-09) | C 49 "how large the error can be for the particular assignment"; C 162 | Ch3 152 | done |
| FCL-015 | "pilot trial" | Accept (I-08) | C 49 "education pilot trial"; C 147 "pilot trial of education" | Ch3 59, 374 | done (heading C 243 still "pilot") |
| FCL-016 | One-line formula | Accept (M-02) | C 153 one-line equation, also in reading PDF | Ch3 75 | done |
| FCL-017 | Queen/rook sentence came from nowhere | Accept (M-01, M-10) | C 173 "Neighbors were defined by queen contiguity" in Simulation study | Ch3 158 | done |
| FCL-018 | "noise;" | Accept (M-02) | C 156 "The noise ε is independent standard normal, and" | Ch3 78 | done |
| FCL-019 | strikeout of "(so" | Accept (M-02) | C 156 parenthesis gone | Ch3 78 (one parenthetical remains: sigma = 1) | done |
| FCL-020 | "deviation," | Accept (M-02) | C 156 "in units of its standard deviation." | Ch3 78 | done |
| FCL-021 | "baseline outcome" confusing | Accept (M-04) | C 156 "sets how strongly incidence shifts the outcome" | Ch3 78 | done |
| FCL-022 | "choose" but "not estimated" | Accept (M-11) | C 173 "set by the investigator... not estimated at the design stage" | Ch3 162 | done |
| FCL-023 | Sentence too broken | Accept (M-05) | C 158 "We estimated β by maximum likelihood (ML)." then fit sentence | Ch3 146 | done |
| FCL-024 | Reference for interval claim | Accept (M-06) | C 158 `\cite{mackinnon...}`, `\cite{leyrat_cluster_2018}`, "SI Section 15.2" | Ch3 148 | partial (Wald not cited) |
| FCL-025 | Add "mean squared error" | Accept (M-07) | C 158 "mean squared error (MSE)" (C 49 also spells it out) | Ch3 59 | done |
| FCL-026 | "confidence intervals" | Accept (M-07) | C 158 "95% confidence intervals" | Ch3 148 | done |
| FCL-027 | "coverage was measured" | Accept (M-07) | C 158 "coverage was measured as the proportion" | Ch3 168 | done |
| FCL-028 | identifiable? | Accept (M-08, R-03) | C 160, 197, 220 "identifiable" | Ch3 150 done; SI still "identified" | partial |
| FCL-029 | Delete one-assignment sentence | Accept (M-09) | C 162 sentence absent | Ch3 done (T-23) | done |
| FCL-030 | Why not γ S(Z)? | Accept (P-01) | C 153 "\psi\,\mathbf{S}(\mathbf{Z})" | Ch3 75 | done |
| FCL-031 | Comma after "regimes" | Accept (M-12) | C 173 "both spillover regimes, and" | n/a | done |
| FCL-032 | Write spatial model for X | Accept (M-13) | C 173 "X=\Phi\{(I-\rho_X W)^{-1}u\}" | Ch3 91 | done |
| FCL-033 | One normal, one Poisson? | Accept (M-13) | C 173 spatial incidence normal; Poisson counts with log-rates | Ch3 91 | done |
| FCL-034 | How? | Accept (M-13) | C 173 formula and rank scaling "(r_i-0.5)/100" | Ch3 91 | done |
| FCL-035 | "surface k assignment" | Accept (M-15) | C 173 "Within surface $k$, the assignments, the simulated outcomes" | Ch3 91-97 | done |
| FCL-036 | ρ_X? | Accept (M-13) | C 173 "$\rho_X=0.2$ or $0.5$" defined in the model | Ch3 91 | done |
| FCL-037 | Incidence row is not a design | Accept (E-01) | E 21 "Example incidence surface (not a design)" | n/a | done |
| FCL-038 | "region" instead of "cluster" | Declined | E unchanged, "cluster" | n/a | rejected (author) |
| FCL-039 | Understand NC Checkerboard | Accept (M-22) | E 22, C 187 "visits clusters from most to fewest neighbors" | Ch3 795 | done |
| FCL-040 | k-means on location | Accept (M-22) | C 187 and E 25 "k-means clustering on four standardized variables" | Ch3 416, 799 | done |
| FCL-041 | "so" | Accept (M-16) | C 175 "and therefore carries simulation (Monte Carlo) error" | Ch3 done | done (other "so" remain, see W-03) |
| FCL-042 | "configuration" confusing | Accept (M-13) | C 173 "An incidence configuration is the model that generates incidence" | Ch3 91 | done |
| FCL-043 | Is this a result? | Move declined, reword (M-17) | C 177 "We also ran a separate allocation-risk study on the grid" | Ch3 172 | done (move declined) |
| FCL-044 | Isolation Buffer fair? | Reply only | already C 232 | Ch3 698 | no change (plan) |
| FCL-045 | What is "simulated"? | Accept (M-19) | C 180 "Observed incidence and geography were fixed inputs, and outcomes were simulated" | Ch3 452 | done |
| FCL-046 | "(1)" | Accept (M-20) | C 185 "Equation~(\ref{eq:grid})" | Ch3 uses "Equation \ref" (CTJ only item) | done |
| FCL-047 | Unclear τ and β values | Accept (M-21) | C 185 "The education outcome is not expected to depend on SUD incidence" | Ch3 162, 448 | done |
| FCL-048 | Section order, name | Rename, keep placement (M-18) | C 179 "Case study: sudden unexpected death in North Carolina" | Ch3 372 | done |
| FCL-049 | Move configuration sentence | Accept (R-01) | C 197 last sentence "Incidence configuration changed each design's MSE by less than 10%" | Ch3 n/a | done |
| FCL-050 | Pooled 0.091 not in Table 2 | Accept (R-02) | C 197 pooled sentence gone; "The two spillover regimes gave opposite results." | Ch3 keeps pooled in `tab:performance` | done |
| FCL-051 | What does it mean | Accept (R-02) | C 197 deleted | n/a | done |
| FCL-052 | "(SI)" unclear | Accept (R-04) | C 204, 206, 210, 215, 222 cite SI Figure S8, Sections 6, 12.3, 14, Table S17 | n/a | done (all targets exist in the SI PDF) |
| FCL-053 | Why β=1 (old notation) | Accept (M-11, M-21) | C 173 "same size as the primary intervention effect"; C 185 | Ch3 162 | done |
| FCL-054 | Understand regions | Accept (M-22, E-03) | C 187 "We ran k-means clustering on the cluster coordinates 250 times"; E 27 | Ch3 414 | done |
| FCL-055 | "spatial dependence coefficient" | Accept (E-04) | E 38 and E 73 "spatial dependence coefficient $\rho$" | Ch3 captions done | done |
| FCL-056 | "true γ?" | Accept (P-04) | E 38 "includes the true spillover covariate S(Z), with ψ estimated" | Ch3 146 | done |
| FCL-057 | N/A to "Reference" | Accept (E-05) | E 55 "Reference" (both cells) | Ch3 tab:grid-srs, no "N/A" left | done |
| FCL-058 | Both arms: no better design? | Reply only | already C 197, 238 | n/a | no change (plan) |
| FCL-059 | How sensitive is the winner | Accept (R-09) | C 234 "Incidence-Guided Saturation was not best in every setting." | Ch3 243, 599 | done |
| FCL-060 | "the best" or "preferred" | Edited by author (R-07) | C 212 heading "had lower error than SRS under control-only spillover" | n/a | done (author edit) |

### Meeting decisions and revision list

| ID | Asked | Disposition | Landed | Ch3 | Status |
|---|---|---|---|---|---|
| D-1 | Reparameterize: β on **Z**, τ on **X**, ψ S(**Z**) | Applied (P-01; symbol is ψ) | C 153 equation; no γ or τ-on-**Z** in C, E, SI, Ch3 | Ch3 75 | done (but see N-1 for Figure_Legends.md) |
| D-2 | Coefficient on incidence is pre-specified | Accept (M-11) | C 173 "set by the investigator... not estimated at the design stage" | Ch3 162 | done |
| D-3 | Rename Application to Case Study | Accept (M-18) | C 179, 212, 49 to 55; no "application" left in C or E | Ch3 372 heading | partial (cover letter lines 14, 16, SI label, checklist) |
| D-4 | Cite t and Wald intervals | Accept (M-06) | C 158 cites t quantile and SE factor; Wald in SI 15.2 | Ch3 148 | partial (Wald uncited) |
| D-5 | Estimating τ after a trial is future work | Ch3 Future Work, CTJ optional | none in C | Ch3 162, 636 | done (Ch3); CTJ absent by plan |
| M-01 | Apply D-1 incl. Checkerboard estimand | Done | C 160 "the treatment coefficient estimates $\beta-\psi$" | Ch3 150 | done |
| M-02 | One-line equation, ε in prose | Done | C 153, 156 | Ch3 75, 78 | done |
| M-03 | Define W, avoid "ith element" | Done | C 151 "with elements $w_{ij}=1/d_i$" | Ch3 71 | done |
| M-04 | "cluster-level incidence" | Done | C 156 | Ch3 78 | done |
| M-05 | Replace "baseline outcome" | Done | C 156, 51, 191, 228 "shift the outcome" | Ch3 done (only "baseline" left is Chapter 2 risk and imbalance, Ch3 53, 55) | done (SI still says "baseline incidence", SI 66, 92, 2138; SI not in the plan's scope for this item) |
| M-06 | Move value of τ to simulation settings | Done | C 173 "On the grid we set $\tau=1$" | Ch3 162 | done |
| M-07 | Move queen/rook to simulation design | Done | C 173 | Ch3 158 | done |
| M-08 | "identifiable", drop "(cannot be estimated)" | Done in C | C 160, 197, 220 | Ch3 150 | partial (SI "identified" remains) |
| M-09 | Delete one-assignment sentence | Done | C 162 | Ch3 done | done |
| M-10 | Cite t interval, normal interval, SE factor | Done except Wald | C 158 | Ch3 148 | partial (Wald uncited) |
| M-11 | Define m(a) in words | Done | C 162 "we also asked how large the MSE can be for a particular assignment" | Ch3 152 | done |
| M-12 | Define MSE at first use | Done | C 158 (and C 49) | Ch3 59 | done |
| M-13 | Split "We estimated τ..." | Done | C 158 | Ch3 146 | done |
| S-01 | SUDDEN data only in case study | Done (author wording) | C 173 "the North Carolina data only enter the case study" | Ch3 91 | done |
| S-02 | Full generating model for X | Done | C 173 | Ch3 91 | done |
| S-03 | Justify or reduce five configurations | Keep all five, state result | C 197 only the "<10%" sentence; no justification | Ch3 done | partial (CTJ) |
| S-04 | Why incidence on (0, 1) | Done | C 173 "All configurations put incidence on the $(0,1)$ scale" | Ch3 95 | done |
| A-01 | Rename and update every "application" | Done | see D-3 | Ch3 done | partial (cover letter, labels) |
| A-02 | Observed vs simulated | Done | C 180 | Ch3 452 | done |
| A-03 | Explain block construction (k-means) | Done; wording follows the code, not the meeting note | C 187; the blocks use four features, not pairwise distances | Ch3 414, 416, 799; SI 701, 731 | done (deliberate departure from the meeting notes, per plan T-35) |
| I-01 | Cite "Two allocations..." | Done | C 141 | Ch3 49 | done |
| I-02 | Define Checkerboard at first prior-work mention | Done | C 145 | Ch3 61 | done |
| I-03 | Rewrite "strict checkerboard" sentence | Done | C 143 | Ch3 150 | done |
| I-04 | "that work" to prior project | Done, self-contained | C 145 "In previous work on small grids" | Ch3 55 | done |
| I-05 | Revise "The comparison includes..." | Done | C 145 | Ch3 59-61 | done |
| I-06 | "Our application supports planning" | Done (author edit) | C 147 | Ch3 59, 374 | done |
| W-01 | Avoid parentheses | Sweep of edited items | many parentheses remain in untouched text | not audited | not auditable (style) |
| W-02 | Break long sentences | Sweep of edited items | not audited | not audited | not auditable (style) |
| W-03 | Avoid "so" as connective | Plan says done (C, 3) | C 17 and E 5 remain | Ch3 about 1 | partial |
| W-04 | Define terms, no "tractable" | Done | no "tractable" in C | Ch3 none | done |
| Q-01 | Binary spillover only where it exists | Settled; Methods sentence plus appendix | none in C | Ch3 87, 462, 637, 1519; SI | partial (CTJ) |
| Q-02 | What the incidence coefficient means | Done | C 173 "raises the outcome... by one noise standard deviation" | Ch3 162 | done |

### Other items from the sign-off list spot-checked

All of P-01 to P-04, I-01 to I-11, M-01 to M-24, R-01 to R-11, E-01 to E-06 and S-01 to S-09 that touch C or E were
traced above under their FCL or meeting IDs and found applied as approved or as edited by Andrew (I-06, I-10, M-14,
R-07, R-08, R-11). M-21 needed no change and C 185 matches. R-05: C 208 "bias $-0.121$ to $-0.400$" matches the recorded departure.
R-06: the 63% appears in C 210 and in SI (line 2364 of the extracted text). R-10: C 271 matches. Figure labels (P-03):
no Greek letters remain in the compact figure PDFs' text layers, and the code label strings no longer say "True τ" or γ
(comments in codes 14, 19, 20 still say gamma; harmless).

## Verifier output (read-only runs, project root)

### `python3 paper/tools/verify_manuscripts.py`  FAIL

```
Chapter checked: paper/thesis_chapter/Thesis_Chapter.qmd
Traceback (most recent call last):
  File ".../paper/tools/verify_manuscripts.py", line 84, in verify
    raise AssertionError(f'Figure {index}: separate legend differs from master')
AssertionError: Figure 2: separate legend differs from master
```

Cause: `figures/Figure_Legends.md` (last committed 2026-11-10, before the notation swap) has the old captions. Bypassing that one check in
a scratch copy (not saved in the repo) shows legend mismatches for Figures 2, 3 and 4, then a second stop at line 236:
`ValueError: substring not found` for the label `tab:nc-summary-sens`, which is no longer in the chapter (removed by commit `6b8da5a`,
the 2026-11-15 appendix reorganization). Everything after that check, including the table and citation reconciliation, has not run on the current text.

### `Rscript code/tests/test_manuscript_exhibits.R`  FAIL

```
Error: length(reported) == 3L is not TRUE
Execution halted
```

Cause: the regular expression at lines 57-59 searches `CTJ_Manuscript.tex` for `$\tau=2$` and `$\tau=3$`; the text now writes
`$\beta=2$` and `$\beta=3$`, so nothing matches. The test needs the symbol updated before it can check the 0.160/0.160 and
0.160/0.182 values. `git status` after both runs shows no new changes beyond the files that were already untracked.

## Not checked

Figure image content (only text layers), SI prose beyond the items named, `CTJ_Reading.pdf` pagination, and any edit made after
commit `4386bcb`. Chapter 3 was checked by source text, not by rendering.
