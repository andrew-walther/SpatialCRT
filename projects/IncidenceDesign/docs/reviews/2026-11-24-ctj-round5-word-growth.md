# Where the CTJ body grew past 3,500 words, by revision item

Written 2026-11-24. Author: Andrew Walther. Read-only analysis; no manuscript was edited and nothing was committed.

## Lead

The CTJ body (Introduction to Conclusions, texcount text plus headers, math excluded) grew from 3,498 words at V5 (`00270ae`, 10-27) to 3,853 now (`d8639d5`), which is +355 words. V5 sat 2 words under the limit, so the body is now 353 over. Three commits caused all of it.

| Commit | What | Measured change |
|---|---|---|
| `7e51487` | Notation swap only | 0 |
| `43b6f99` | First-round revisions from Dr. Lin's feedback | +296 |
| `09803e1` | Wording fixes from the independent check | +11 |
| `324cd80` | "pilot study" to "pilot trial" | 0 |
| `d8639d5` | Chapter 3 committee review (C-01 to C-31) | +48 |

The five items that added the most words (net, all commits) are the following.

1. **FCL-022 and FCL-053 (meeting D-2, Q-02)**, +61. The incidence coefficient got its own explanation in the Simulation study ("set by the investigator", "raises the outcome by one noise standard deviation").
2. **FCL-059 (R-09)**, +55. A new Discussion paragraph, "Incidence-Guided Saturation was not best in every setting."
3. **S-04 with committee item C-31**, +35 (+18 and +17). Why incidence is on the (0, 1) scale, then a caveat added by the committee round.
4. **FCL-003 and FCL-010 (I-03, I-04)**, +31. The rewritten checkerboard sentences in the Introduction.
5. **FCL-012 and FCL-015 (I-06)**, +30. The Introduction's case-study sentence (+23) and the Continued work heading and sentence that follow from it (+7).

Several large additions are explanatory text that goes beyond the reviewer's literal request (the FCL-059 paragraph, the S-04 scale sentence, the CC-09 sidedness note and the committee caveats C-16, C-28 and C-31, about 150 words in all). Offsetting cuts (queen/rook move, writing the incidence models as formulas, deleting the one-assignment sentence, deleting the pooled-MSE sentence, one equation copy) removed about 67 words, so gross additions were near 420.

## How the counts were made

For each commit pair I took `git diff --word-diff` on `CTJ_Manuscript.tex`, rebuilt the old and new text of every hunk, and ran texcount on each side (text plus headers, so `\cite`, `\ref` and math do not count; "21.1--29.9" counts as two words, as in `paper/tools/ctj_wordcount.py`). Within rewritten paragraphs I split old and new text into sentences and counted each one, then gave each sentence to the item that caused it. Where one sentence serves two items I split it by judgment, so item totals carry about +/- 3 words of allocation error. Moves between sections count as a removal in one place and an addition in the other, and the item gets the net.

Reconciliation against the measured counts.

| Step | Measured (ctj_wordcount.py) | Hunk-by-hunk texcount | Item totals below |
|---|---:|---:|---:|
| `7e51487` to `43b6f99` | +296 | +288 (2.7% low; one display-equation hunk is mis-counted by texcount, about 13 words) | +296 |
| `43b6f99` to `09803e1` | +11 | +11 | +11 |
| `324cd80` to `d8639d5` | +48 | +48 | +48 |
| Total | +355 | +347 | +355 |

Section totals from the full-body texcount agree with the hunk sums to within 5 words in every section (Introduction +68, Outcome model and estimation -53, Simulation study +150, Case study methods +49, Results case study +17, Implications +57 for the first round). The item totals for the first round were assigned so that the section totals match; the independent check on them is the hunk sum, which is within 3%. `00270ae` to `7e51487` is notation only. I saw only symbol swaps, rounding of MSE values to three decimals (no word change) and the declarations, which lie outside the body.

## Table by item

"Essential" means the reviewer's comment cannot be answered without it. "Shorten" means a shorter version still answers it. Judgments are mine and brief. The trim plan (`word-limit-trim-plan.md`) T-numbers are given where they exist.

| ID | What it asked | Words | Where | Essential or shorten |
|---|---|---:|---|---|
| FCL-022, FCL-053 (D-2, Q-02, PE M-11) | "choose" the incidence coefficient but "not estimated"; why 1; what it means | +61 | Simulation study (moved from Outcome model, 97 new against 36 old) | Shorten. Dr. Lin needs "set by the investigator" and the one-SD meaning (about 45 words). The confounding and "may need adjusting" sentences repeat the old text. |
| FCL-059 (PE R-09) | How sensitive is the winning design | +55 | Implications | Shorten. Two of its four sentences restate Results (the 8% and the 22% to 5%, trim plan T13). Keep one sentence. |
| S-04 and C-31 | Why incidence is on (0, 1); then the committee's endpoint caveat | +18, +17 (35) | Simulation study | Shorten. One sentence can carry both ("sets the same shift between the lowest and highest incidence"). Trim plan T11 already proposes deleting it. |
| FCL-003, FCL-010, FCL-001, FCL-009 (I-03, I-04) | Unclear checkerboard sentence; citation placement | +31 | Introduction | Essential in substance. The two new sentences explain the checkerboard in plain words, but the second could be merged into the first (about 12 words). |
| FCL-012, FCL-015 (I-06, I-08) | Rephrase "we aim to plan"; "pilot trial" | +23 Introduction, +7 Continued work (30) | Introduction; Continued work heading and sentence | Shorten the Introduction part. "statewide extension" is repeated in three places (abstract, Introduction, Continued work). |
| FCL-039, FCL-040, FCL-054 (M-22, A-03) | Explain NC Checkerboard, k-means blocks, regions | +23 | Case study methods | Essential. The reviewer could not follow the three rules. Trim plan T6 would move the rounding rule to Table 1 instead. |
| C-16 | Committee: factor attributed to HC1 | +23 | Outcome model | Shorten. Keep the caveat "approximate and chosen by calibration (SI Section 15.2)" and drop "as the HC1 covariance estimator does". |
| FCL-024 (D-4, PE M-06) | Reference for the interval claim | +20 | Outcome model | Essential. Two citations and a pointer to the SI answer a direct reviewer request; the Wald definition "q = z" could be dropped. |
| FCL-042 | "configuration" is confusing | +20 | Simulation study | Essential. Defines configuration and surface, which the rest of the section uses. Could be 12 words. |
| CC-09 (PE M-24, not Dr. Lin) | One-sided grid power against two-sided case-study power | +20 | Case study methods | Shorten. Trim plan T2 already moves detail to the SI. One clause ("two-sided, unlike the one-sided grid power in SI Section 4") does it. |
| FCL-005, FCL-006, FCL-007, FCL-011 (I-02, I-05) | What does this sentence mean; define Checkerboard at first use | +16 | Introduction | Essential. Trim plan T14 would delete "We do not compare numerical results" (10 words). |
| FCL-052 (R-04) | Replace each bare "(SI)" with a section or table | +11 | Results (6 places) | Essential and already minimal (+1 to +2 per place). |
| M-03 (meeting) | Define **W** without "the ith element" | +8 | Outcome model | Essential; the definition is needed. |
| D-3, A-01, FCL-048 | Rename "application" to "case study" | +8 | 8 places in the body | Unavoidable (two words replace one). Using "case study" as an adjective in headings costs the same. |
| W-04 (no item ID found) | Plain-language gloss of the elementwise product | +7 | Outcome model | Shorten. "which multiplies the two vectors cluster by cluster" can go if the symbol is explained once. |
| FCL-025, FCL-026, FCL-027 (M-07, PE M-07) | Spell out "mean squared error", "confidence", "was measured" | +7 | Outcome model | Essential (definitions at first use). |
| M-05, FCL-021 (PE M-23) | Replace "baseline outcome" in the sensitivity analysis | +7 | Case study methods | Essential; the new wording is longer but clearer. |
| R-08, I-10 | State which designs the coverage range covers | +5 | Results case study | Essential for accuracy. |
| C-28 | Committee: saturation recommendation needs the participation condition | +5 | Implications | Essential for the claim. |
| S-01, FCL-035 | SUDDEN data only in the case study; "Within surface k" | +4 | Simulation study | Essential. |
| C-13 | Committee: exposure mapping citation | +3 | Outcome model | Essential (adds two citations). |
| FCL-060 | Heading "best" or "preferred" | +3 | Results case study heading | Essential (heading). |
| FCL-018, FCL-019, FCL-020 (M-02) | Noise and "units of its SD" in prose | +3 | Outcome model | Essential. |
| FCL-043 | Is the allocation-risk paragraph a result | +2 | Simulation study | Essential. |
| FCL-028 (M-08, R-03) | "identifiable" | net +2 (-7 in Methods, +4, +3 and +2 in Results) | Outcome model; both Results sections | Essential. |
| W-03 and independent check (09803e1) | "so" to "which means" or "and therefore"; "before spatial spread"; "a direct effect:" | +11 | Outcome model +4, Simulation study +5, Case study methods +2 | W-03 part (+4) is a style rule, the rest answers the independent check; all could be undone for words. |
| FCL-023 | Split the estimation sentence | +1 | Outcome model | Essential. |
| R-05 | "t-based" coverage label | +1 | Results | Essential. |
| FCL-016 | One-line equation | about -5 | Outcome model | Cut. |
| FCL-017 (M-01, M-10, meeting M-07) | Move queen/rook to Simulation study | -11 (+28 new, -39 old) | Outcome model to Simulation study | Cut; the move shortened the text. |
| FCL-032 to FCL-036, S-02 (M-13) | Write the models for the spatial and Poisson incidence | -19 (formulas replace prose) | Simulation study | Cut; formulas are shorter than the prose they replaced. |
| FCL-014, FCL-029 (M-09, M-11) | Define allocation risk; delete the one-assignment sentence | -17 | Outcome model | Cut. |
| FCL-050, FCL-051 (R-02) | Pooled 0.091 is not in Table 2 | -15 | Results simulation study | Cut. |
| Sum of items | | +355 | | |

The abstract and declarations changed too (the abstract by about +27 words) but are outside the body count and not listed.

## Passage detail

Quotes are 12 words or fewer. Words are texcount words, plus for added, minus for removed.

### FCL-022, FCL-053 (D-2, Q-02): the incidence coefficient, +61

- Removed from Outcome model, -36: "We set tau=1 in the grid study so that incidence is prognostic" and "Investigators can choose tau to suit their own outcome."
- Added to Simulation study, +27: "The incidence coefficient is set by the investigator according to whether incidence is expected to shift the outcome" (answers FCL-022, D-2).
- Added, +5: "On the grid we set tau=1."
- Added, +32: "Incidence lies between 0 and 1, so this value raises the outcome by one noise standard deviation" (answers Q-02, FCL-053).
- Added, +33: "With incidence shifting the outcome, treatment that follows incidence is partly confounded with it" and "The value may need adjusting for a particular outcome." This restates the removed text, so the net of the pair is about +2.
- Later, `09803e1`, +3: "before spatial spread".

### FCL-059 (R-09): new Discussion paragraph, +55

- Implications, +55: "Incidence-Guided Saturation was not best in every setting." (8), "Under both-arms spillover its MSE was 8% above that of SRS" (14), "its advantage narrowed as beta grew, and at beta=3 Isolation Buffer had lower MSE" (15), "its advantage over SRS fell from 22% with 58 colleges to 5% with 29" (18). The last two repeat Results (the 22% to 5% figure appears in Results and again here; trim plan T13).

### S-04 and C-31: scale of incidence, +35

- Simulation study, `43b6f99`, +18: "All configurations put incidence on the (0,1) scale, so one value of tau means the same shift in every configuration."
- Simulation study, `d8639d5`, +17: "although the spread between those limits differs somewhat by mechanism" and "between the lowest and highest possible incidence".

### FCL-003, FCL-010, FCL-001, FCL-009, FCL-004 (I-03, I-04): Introduction checkerboard and design-stage sentences, +31

- Removed, -21: "In a strict checkerboard with edge neighbors, for example, a cluster's exposure to treated neighbors is fully determined by its own treatment."
- Added, +29: "For example, in a checkerboard design with neighbors defined by shared edges, every neighbor of a treated cluster is a control and every neighbor of a control is treated."
- Added, +24: "A cluster's share of treated neighbors is then fixed by its own arm, and the treatment effect cannot be separated from the spillover effect." The first sentence ("Design-stage methods, including covariate-constrained randomization ...") changed by about -1 and moved its citations to the end (citations do not count).

### FCL-005, FCL-006, FCL-007, FCL-011 (I-02, I-05): Introduction prior-work paragraph, +16

- Added, +7: "arm; we call this design Checkerboard here" (FCL-011, I-02).
- Added, about +9 net: "In that study, we judged each design by its average estimation accuracy and by the risk" replaces "That work weighed ... against the risk of a poor realized allocation" (+14 new, -14 old but split into longer sentences); "In this paper, we compare six designs and SRS on the mean squared error" replaces "The comparison includes ... (SRS) ..." (the definition of SRS moved earlier, net 0); "We do not compare numerical results between the two studies." became its own sentence (+10, the earlier clause was about -3). Allocation among FCL-005, -006 and -007 is not separable.

### FCL-012, FCL-015 (I-06, I-08): case-study sentence and its echoes, +30

- Introduction, +23: "In this paper, we also use a case study, a simulation based on observed incidence, to recommend a design for a statewide extension of a pilot study of" replaces "Our application supports planning an education intervention for" (new 39, old 16). `324cd80` changed "study" to "trial" at no word cost.
- Continued work, +7: heading "a statewide extension of" (+4) and "an allocation for the statewide extension" replacing "a pilot allocation" (+3).

### FCL-039, FCL-040, FCL-054 (M-22, A-03): NC design adaptations, +23

- Checkerboard, +8: "visits clusters from most to fewest neighbors and gives each cluster the arm held by fewer" replaces "alternates arms between neighboring clusters as far as the map allows" (26 against 18).
- Incidence-Guided Saturation regions, +11: "We ran k-means clustering on the cluster coordinates 250 times and kept the partition" replaces "chosen from 250 k-means candidate partitions" (33 against 22).
- 2x2 Blocking, +4: "k-means clustering on four standardized variables (location, that year's incidence rank" (59 against 55).

### C-16, FCL-024 (D-4, M-06), C-13: interval and exposure-mapping sentences

- FCL-024, Outcome model, +20: "The factor replaces the ML variance divisor N with N-p" with the MacKinnon citation, "the t quantile allows for the small number of clusters" with the Leyrat citation, "In SI Section 15.2, this interval covered beta more reliably than the Wald interval", and "which uses q = z = 1.96 without the factor" (58 new against 38 old for the interval sentences).
- C-16, Outcome model, +23: "as the HC1 covariance estimator does, although in the spatial lag model this rescaling is approximate and was chosen by calibration (SI Section 15.2)."
- C-13, Outcome model, +3: "an exposure mapping" (two citations add no words).

### FCL-042, S-01, FCL-035: Simulation study incidence definitions, +24 net of the incidence models

- FCL-042, +20: "An incidence configuration is the model that generates incidence over the grid; an incidence surface is one draw from it."
- FCL-032 to FCL-036, S-02, -19: the three prose sentences on spatial and Poisson incidence (85 words) became formulas and one sentence (66 words), for example "spatial incidence X = Phi{(I - rho_X W)^-1 u}, with u independent standard normal".
- S-01, FCL-035, +4: "the North Carolina data only enter the case study" replaces "rather than a calibration to the application", and "Within surface k," starts a new sentence.

### CC-09 (PE M-24): sidedness of power, +20

- Case study methods, `43b6f99`, +17 and +1 ("t-based"): "a two-sided test at level 0.05; the grid power in SI Section 4 is one-sided at level 0.025".
- `09803e1`, +2: the parenthesis was recast as "This is a two-sided test ..." as a separate sentence.

### FCL-052 (R-04): SI pointers, +11

Each "(SI)" became a named target, adding 1 to 2 words in six places: "(SI Figure S8)", "(SI Section 6)", "(SI Section 12.3)", "SI Table S17", "(SI Section 14)" and "SI Table S5".

### M-03, FCL-018 to FCL-020, W-04, FCL-016, FCL-023, FCL-025 to FCL-027: Outcome model paragraph

- M-03, +8: "with elements w_ij = 1/d_i when clusters i and j are neighbors" and "where d_i is the number of neighbors of cluster i" (42 words against 34).
- FCL-018 to FCL-020, +3: "The noise is independent standard normal, and the intervention effect is therefore in units of its standard deviation."
- W-04, +7: "which multiplies the two vectors cluster by cluster."
- FCL-016, about -5: the second (reading-layout) copy of the equation was removed; texcount counted words inside that branch.
- FCL-023, +1: "We estimated beta by maximum likelihood (ML)." and "The fitted SAR lag model included an intercept ..." (28 against 27).
- FCL-025 to FCL-027, +7: "the mean squared error (MSE)", "coverage was measured as the proportion of 95% confidence intervals" (27 against 20).

### Cuts

- FCL-017, -11 net: "Neighbors were defined by queen contiguity (clusters sharing an edge or a corner)" (+28 in Simulation study) replaces the two sentences in the Outcome model (-39).
- FCL-014, FCL-029, -17: "A design with only one possible assignment has no spread, although its error can still be high." (-18) deleted; "A trial realizes one assignment, so we also measured" became "For each design, we also asked how large the MSE can be" (about +1).
- FCL-050, FCL-051, -15: "Pooled over both regimes, Incidence-Guided Saturation had a mean MSE of 0.091" deleted; "The two spillover regimes gave opposite results." added (+7).
- FCL-028, net +2: Methods -7 ("treatment and spillover cannot be separated" and "(cannot be estimated)" deleted), Results +4 ("under rook contiguity Checkerboard's was not identifiable"), +3 ("every design whose effect was identifiable"), +2 ("Checkerboard's effect was identifiable").

### Smaller additions

- D-3, A-01, FCL-048 rename, +8, about +1 per use: "case study" for "application" in the Designs table sentence, the Case study methods, three Results and Discussion sentences, Implications (two places) and Conclusions.
- M-05 and FCL-021 (PE M-23), +7: "let incidence shift the simulated outcome (tau=1) and added incidence as a covariate in the fitted model."
- R-08 and I-10, +5: "Across the six designs and SRS and both regimes."
- FCL-060, +3: heading "had lower error than SRS under control-only spillover" replaces "is strong".
- C-28, +5: "controls and most colleges take part" in the Implications recommendation.
- FCL-043, +2: "We also ran a separate allocation-risk study on the grid, using ..."
- W-03 and the independent check (`09803e1`, +11): "so" became "which means" and "and therefore" (+4); "a direct effect:" (+2, Outcome model); "before spatial spread" (+3, Simulation study); "This is a" (+2, Case study methods).

## Notes for the trim

- The three largest first-round additions (FCL-022/-053, FCL-059, S-04 with C-31) total 151 words, and each has a shorter version that still answers the reviewer. Cutting them to about 85 words would save about 65.
- Trim plan items T11 (S-04, 25 words), T13 (FCL-059 repeat, 20 words), T14 (10 words), T2 (CC-09 detail, 90 words) and T8 (Introduction) line up with the items above and are the cheapest cuts because they restate text already elsewhere in the body.
- Nothing here suggests an addition was wrong. The words come from the reviewer's requests to define terms, give sources and state the model, so the body was at the limit before any definition could be added.
