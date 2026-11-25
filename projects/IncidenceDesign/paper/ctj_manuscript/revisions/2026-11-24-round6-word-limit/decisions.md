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

N5 final wording (author, 2026-11-25): "The investigator specifies $\tau$ according to whether incidence is expected to shift the outcome. It is not estimated at the design stage and we set $\tau=1$ on the grid simulation. Incidence ranges from 0 to 1, and this value raises the outcome by one noise standard deviation from the lowest to the highest incidence before spatial spread. That shift equals the primary intervention effect $\beta=1$. Treatment that follows incidence is then partly confounded with it, the harder case for designs that use incidence."

## 1c. Human-writer list (`human-writer-review.md`; IDs as in that file)

Settled by cut decisions (not applied): HW-04 (N8 rejected), HW-06 (= C9, rejected), HW-14 and HW-15 (N5), HW-18 (T11), HW-32 (T7), HW-35 (T3). HW-20 rebased on the kept N9 sentence (batch 2).

Batch 1 (decided 2026-11-25)

| ID | Decision |
|---|---|
| HW-01, HW-02, HW-07, HW-09, HW-10, HW-12, HW-13 | Accept as proposed |
| HW-03 | Edit: "When spillover is expected mainly in control clusters, Incidence-Guided Saturation had the lowest or near-lowest estimation error and lower allocation risk than SRS. It also gave the highest-burden region the largest treated share." |
| HW-05 | Edit: "…in the same arm.\cite{walther_spatialcrt_2026} We refer to this design here as Checkerboard." |
| HW-08 | Edit: "A trial draws only one assignment. Therefore, we evaluated both…" |
| HW-11 | Edit: "…shows an example assignment of each on the grid, and summarizes its adaptation to the case study. The supplementary information (SI) gives full specifications." |

N5 comma (author, 2026-11-25): "It is not estimated at the design stage, and we set $\tau=1$ in the grid simulation."

Batch 2 (decided 2026-11-25)

| ID | Decision |
|---|---|
| HW-16, HW-17, HW-21, HW-22, HW-23, HW-25, HW-26, HW-27 | Accept as proposed |
| HW-03b (new) | Accept: the Incidence-Guided Saturation description moves to its first abstract mention, in Results: "On the grid with treatment effect $\beta=1$, Incidence-Guided Saturation, which treats larger shares of clusters in higher-incidence regions, had a control-only mean MSE of 0.135…" (HW-03 keeps it out of Conclusions) |
| HW-19 | Edit: "These are general simulation settings. The North Carolina data enter only the case study." |
| HW-20 | Reject (N9 kept) |
| HW-24 | Edit: "Therefore, treatment budgets differed between designs (Table~\ref{tab:allocation-rules}). Rather than equalizing budgets after the fact, we report each design's number of treated clusters…" |

Batch 3 (decided 2026-11-25)

| ID | Decision |
|---|---|
| HW-30, HW-31, HW-34, HW-37, HW-38, HW-39 | Accept as proposed |
| HW-28 | Edit (keeps "nearly equal"; unrounded 0.1596 and 0.1604): "At $\beta=2$, Isolation Buffer and Incidence-Guided Saturation had nearly equal MSE (0.160 and 0.160). At $\beta=3$, Isolation Buffer was lower (0.160 versus 0.182), because it largely avoids the error in estimating $\rho$ that for most designs passes into $\hat\beta$ in proportion to $\beta$ (SI Section 6)." |
| HW-29 | Edit: as proposed, with the last sentence "Therefore, a low MSE alone does not show that a design supports valid inference." |
| HW-33 | Edit: "The two changes were made together. Therefore, their separate contributions cannot be distinguished." |
| HW-36 | Edit: "…and investigators need to judge which is more plausible when planning." |
| TH-1 (new) | Accept: "The noise $\boldsymbol{\varepsilon}$ is independent standard normal. Therefore, the intervention effect $\beta$ is in units of its standard deviation." |
| TH-2 (new) | Accept: "Therefore, a design's ranking on the grid need not carry over unchanged to an irregular geography." |

Batch 4, captions and Table 1 cells (decided 2026-11-25)

| ID | Decision |
|---|---|
| HW-41 to HW-47, HW-49 to HW-51 | Accept as proposed |
| HW-40 | Edit: "(blue = treated, white = control)" |
| HW-48 | Edit: "Designs share surfaces and outcome noise. Therefore, Table~\ref{tab:grid-srs} gives…" |
| TC-1 (new) | Accept: Checkerboard grid cell "…in a checkerboard pattern, in which no two clusters that share an edge are in the same arm (50 treated)." |
| TC-2 (new) | Accept: Isolation Buffer grid cell "…until none remain. No two treated clusters are then neighbors (about 21 treated with queen neighbors)." |

Checks on the fully edited scratch text (2026-11-25): body 3,419 words; `repeated_openers.py` (bios-dissertation human-writer review) finds no runs of three or more sentences, or two or more paragraphs, opening the same way; eight sentences open with "Therefore," and one with "However,", none adjacent to another.

All decisions are applied by `round6_apply.py --in-place` in 1d.

## 1d. Applied (2026-11-25)

`round6_apply.py --in-place` applied every accepted item; the Lee (2004) entry was added to
`paper/SpatialCRT_IncidenceDesign.bib`; SI Section 5.2 got the T6b rounding sentence;
`figures/Figure_Legends.md` was regenerated from the master captions; R2 was re-checked against
`results/grid_tdf_2026-10/grid_tdf_scenarios.rds` (oracle, queen, β = 1): Incidence-Guided
Saturation had the lowest control-only MSE in all five configurations (0.129 to 0.139). Removed
text and its new home: `paper/ctj_manuscript/notes/trimmed-material.md`. Chapter 3 candidates:
`ch3-carryover-candidates.md`.

Counts after applying: body 3,419 (texcount text + headers, math excluded) or 3,536 with each
math group counted as a word; abstract 433 (limit 425), pending an abstract decision.

Abstract and counting (decided 2026-11-25)

| ID | Decision |
|---|---|
| AB-1 | Accept: "We then applied the designs in a case study for a statewide extension of a sudden unexpected death (SUD) prevention education pilot trial in North Carolina." |
| AB-2 | Accept: "Sensitivity analyses used rook instead of queen contiguity, let incidence shift the outcome, and ordered regions by alternative incidence summaries." (the queen sentence is folded in) |
| Counting | Math and equations are excluded from the 3,500 count (the journal rule excludes abstract, references, tables and figures and does not mention equations); `ctj_wordcount.py` text + headers is the reported number |
| Title page | "Main-text word count: 3,419. Abstract word count: 422." |

Final V7-draft counts: body 3,419, abstract 422 (limit 425). Deferred to the post-submission
thesis: carry the prose revisions (not the cuts) into `bios-dissertation/thesis/` Chapter 3 and
its appendix (`docs/TODO.md`).

## 1e. Chapter 3 carry-over (`ch3-signoff.md`; decided 2026-11-25, done in this session at the author's request)

| ID | Decision |
|---|---|
| CH-1 to CH-4 | Accept as proposed (Lee cite, TAU-1, N5, T6b rounding) |
| TS-01 to TS-35 | Accept as proposed ("Therefore," at sentence start; TS-07 and TS-15 are covered by CH-3 and CH-4) |
| TS-22 | Edit (softer wording): "(Table \ref{tab:downside}). Therefore, they are the designs to consider first in that setting. This guidance assumes that most colleges take part, …" |

Applied by `ch3_apply.py` to `paper/thesis_chapter/Thesis_Chapter.qmd` (the hook regenerates the
`bios-dissertation/prelim/` copy) and to `bios-dissertation/thesis/chapter3-project2-incidence/project2-incidence.qmd`.
Lee (2004) added to `bios-dissertation/prelim/references.bib` and `thesis/references.bib`. The frozen
`prelim-submitted-2026-11-15/` copy is unchanged.

## 1f. Cover letter and checklist (decided 2026-11-25)

CL-1, CL-2 (cover letter: "case study", 3,419 and 422 words, "draws", Incidence-Guided Saturation named
for the 21–30%, Wald undercoverage) and CK-1 (checklist IRB note) accepted as proposed. The SI label
`observed-incidence-application-full-methods-and-results` is internal and left unchanged.
