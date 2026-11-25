# CTJ round 6: author decisions (Phase 1a to 1d)

Author: Andrew Walther. Decisions log for round 6, written as they are made. The 1h review
package for Dr. Lin is built only from items accepted here.

## 1a. Open feedback items (decided 2026-11-25)

| ID | Item | Decision |
|---|---|---|
| A1 | Reference for the Wald interval (FCL-024, D-4, M-10) | Cite Lee (2004, *Econometrica* 72:1899–1925, DOI 10.1111/j.1468-0262.2004.00558.x), which proves the ML estimator of the spatial lag model asymptotically normal. One `\cite` after "the Wald interval" in Methods (C 158), no added words. Applied in 1d with a new bib entry in the shared bib and in `bios-dissertation/prelim/references.bib`. Chapter 3 gets the same cite (Ch3 148) in 1e. Ord (1975) was considered first, but its information-matrix result could not be confirmed from the sources read at the time. |
| A2 | `docs/plans/2026-10-22-ctj-structure-and-prose-review-plan.md` | Dropped; status set to superseded by the phased plan. Its repeated-openers check is run once on the post-HW text in 1c. |
| A3 | Audit partial items S-03 (five configurations) and Q-01 (binary spillover) | Decided in 1b as R2 and R3. |
| A4 | D-5 (CTJ sentence on estimating τ) | Left to Phase 4e. |
| A4 | "pilot" in the Continued-work heading | Kept. |
| A4 | M-08 / FCL-028 ("identified" in the SI) | Phase 2c. |
| A4 | D-3 ("application" in cover letter and checklist) | Phase 1f. |
| A4 | W-03 (", so" joins) | Handled by the HW list. |
| A4 | C-16 (HC1 clause in C 158) | Decided with T12 and C10 in 1b. |
| A5 | `verify_manuscripts.py` stops at `Supplement derivation: tab:designs` (Chapter 3 vs SI table drift, Phase 2d) | Report the failure exactly, and also run the remaining CTJ checks from a scratch copy that skips only the chapter-vs-SI table loop, until the SI and Chapter 3 are aligned. |

Side finding: the Zotero attachment for `ord_estimation_1975` is Ord, Koehler and Snyder (1997),
from a mislabeled ResearchGate download. The author has the correct 1975 copy from JSTOR.

## 1b. Cut list (`signoff.md`)

Batch 1 (decided 2026-11-25)

| ID | Decision |
|---|---|
| T14, T12, C10, C3, N2, C4, T4 | Accept as proposed |
| N8 | Reject (keep the two checkerboard sentences, FCL-010) |
| N9 | Reject (keep "rather than treating settings that share a surface as independent replications") |
| C-16 (new) | Edit C 158: "as the HC1 covariance estimator does," → "a standard small-sample correction," (MacKinnon cite kept) |
| T11 | Edit: "…to put them on the same $(0,1)$ scale as the other configurations, so that one value of $\tau$ shifts the outcome by the same amount from the lowest to the highest incidence in every configuration." (the following "All configurations put incidence…by mechanism." sentence is deleted) |
| N3 | Edit: "Each performance metric is an average over a finite number of simulated trials, and its Monte Carlo standard error (MCSE) measures the simulation error in that average.\cite{morris_using_2019}" |

Batch 2 (decided 2026-11-25)

| ID | Decision |
|---|---|
| T5, C1, T6, N7, N6, N12b, T2, T2b, C5 | Accept as proposed |
| T6b | Accept, plus SI line 733: "Half of each group, rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), is treated at random, so a one-cluster group is never treated." |
| T10 | Edit: "For Isolation Buffer, the two regimes define the same model because no treated cluster has a treated neighbor, and its small differences between regimes are Monte Carlo error." |
| T7 | Edit (absorbs HW-32): "However, each worst decile held only ten assignments. Under both-arms spillover, most of the variance in Incidence-Guided Saturation's estimated risk was Monte Carlo error (SI Section 14). Therefore, small both-arms differences do not show better downside protection." |

Author style note (2026-11-25): "However" and "Therefore" open the sentence rather than sitting mid-sentence, where the prose allows; use articles carefully and never a bare "NC" as a descriptor ("the NC case study").

Batch 3 (decided 2026-11-25)

| ID | Decision |
|---|---|
| C8, N12, T13, C6, N13, C7, R2, R3 | Accept as proposed (R2 re-checked against the results file before applying; R3's SI Section 15.4 is "Sensitivity analyses", which compares proportional, count and binary spillover in the case study) |
| T3 | Edit: "The low both-arms error of High Incidence Focus in the case study (Figure~\ref{fig:nc-mse}) came from a fixed assignment and was sensitive to incidence adjustment. Therefore, it does not support a general recommendation. Which regime applies depends on how an intervention spreads between neighboring areas, and investigators need to judge this when planning. Saturation keeps both arms in every region while varying the local share treated, which helps separate the direct effect from neighboring exposure." |

Batch 4 (decided 2026-11-25)

| ID | Decision |
|---|---|
| TAU-1 (new) | Accept. End of C 156: "The incidence coefficient $\tau$ sets how strongly incidence shifts the outcome, holding treatment and spillover fixed. With $\tau=0$, incidence guides only the allocation, as in the case study. With $\tau>0$, clusters with higher incidence have higher outcomes, as on the grid." |
| N5 | Accept with the third sentence rewritten (no opening "Because", fewer commas); wording in the 1c batch for sign-off |
| C2, C9, T1 | Reject |
| N4 | Reject (meeting item S-01: the North Carolina data enter only the case study) |
| TAU-SI (pending) | One CTJ sentence pointing to the SI sensitivity analysis over a range of $\tau$ values, written in Phase 4e once that SI section and its values exist; keep about 15 words of slack for it |
