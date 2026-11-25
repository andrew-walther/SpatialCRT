# Phase 2e, E4: pointers to unabridged-only material (sign-off)

Author: Andrew Walther. Drafted 2026-11-25. Variant B (E1) moves 29 exhibits and four text
blocks to the unabridged version. The curated text refers to them in 40 places (found by
building the curated selection and listing every undefined reference). Each item below gives
the current text and the curated text. The unabridged build keeps the current text unchanged.
In the source each item becomes `\iffullsi <current>\else <curated>\fi`. Reply by ID.

Five of these exhibits are copies of the article's own exhibits (Figure 1B and 1C, Table 2,
Figures 3 and 4). There the curated text points to the article's exhibit instead of dropping the
pointer (X1).

## Batch A: mechanism, front matter and methods

| ID | Where | Current | Curated |
|---|---|---|---|
| X1 | Preamble (Decision) | — | **Recommend:** the SI imports the article's labels with `xr-hyper` (prefix `M-`, `\externaldocument[M-][nocite]{../manuscript_ctj/.build/CTJ_Manuscript}`), so "Figure 4 of the article" prints the article's live number, like the CTJ's SI references. Build order on a clean checkout becomes SI, CTJ, SI. The alternative is hard-coded article numbers, which you have ruled out. |
| X2 | Overview paragraph (curated only) | "This supplement gives the full methods, results and supporting diagnostics behind the article. Part I …" | "This supplement gives the methods, results and supporting diagnostics behind the article. Part I …" and, after the Part II sentence: "An unabridged version with additional tables and figures, such as per-parameter tables, diagnostic plots and the details of the random-number streams, is available from the authors." |
| X3 | Title and first-page note (unabridged only, E5) | "SI Appendix" | Title "SI Appendix – Unabridged". New first sentence of the overview: "This unabridged version contains every table and figure prepared for the study. It is not the supplement submitted with the article, and its table and figure numbers differ from those the article cites." |
| P1 | 2, designs (SI 163) | "Table S1 describes the six retained designs, and Figure S1 shows one assignment from each, and from SRS, on a single illustrative incidence surface." | "Table S1 describes the six retained designs, and Table~1 of the article shows one assignment from each, and from SRS, on a single illustrative incidence surface." (The article's thumbnails come from the same surface, `code/25`.) |
| P2 | 4, metrics (SI 429) | "Section 4.3 gives the parameter grid, the seed keys, and the exact formulas for these metrics and standard errors." | "Table S3 gives the parameter grid." |
| P3 | 5.1 (SI 656) | "…and report years separately (Figure S3)." | "…and report years separately (Figure~1C of the article)." |
| P4 | 5.2 (SI 707) | "…one fixed four-region partition in all years (Figure S4)." | "…one fixed four-region partition in all years (Figure~1B of the article)." |
| P5 | 5.2 (SI 745) | "Figure S5 shows one example assignment from each strategy." | (sentence omitted) |
| P6 | Table S6 caption (SI 779) | "…using the fixed region IDs in Figure S4." | "…using the fixed region IDs in Figure~1B of the article." (The article's panel B labels Regions 1–4, the same IDs.) |
| P7 | Error decomposition (SI 1003) | "…averaged over scenarios (Section 4.3), which rose…" | "…averaged over scenarios, which rose…" |

## Batch B: grid results

| ID | Where | Current | Curated |
|---|---|---|---|
| P8 | 7.1 (SI 1023) | "…and Checkerboard (1.087) (Table S7, Figure S6)." | "…and Checkerboard (1.087) (Table S7)." |
| P9 | 7.1 (SI 1025) | "Coverage of the 95% intervals was between 0.946 and 0.950 for all six designs (Figure S7)." | "Coverage of the 95% intervals was between 0.946 and 0.950 for all six designs." |
| P10 | 7.1 (SI 1027) | "…computed from the 50 configuration-by-surface means (Section 4.3), and those of bias…" | "…computed from the 50 configuration-by-surface means, and those of bias…" |
| P11 | 7.2 (SI 1096) | "Table S8 gives each design's MSE within each of the five incidence configurations. Under queen contiguity, the mean-MSE order within every configuration…" | "Figure~2 of the article shows each design's MSE within each of the five incidence configurations. Under queen contiguity, the mean-MSE order within every configuration…" |
| P12 | 7.3 (SI 1156) | "…no design's MSE changed by more than 0.02 (Table S9)." | "…no design's MSE changed by more than 0.02." |
| P13 | 7.3 (SI 1160) | "…with no trend in ψ (Table S10)." | "…with no trend in ψ." |
| P14 | 7.3 (SI 1228) | "The spillover regime was the one factor that changed the queen ranking (Table S11, Figure S6)." | "The spillover regime was the one factor that changed the queen ranking (Table~2 of the article)." |
| P15 | 7.3 (SI 1245) | "In the full benchmark that includes SRS (Table S16), however, its both-arms MSE…" | "In the full benchmark that includes SRS (Table~2 of the article), however, its both-arms MSE…" |
| P16 | 8 (SI 1371) | "…separated 17 of the 21 pairs (Figure S9 shows the critical-difference diagrams)." | "…separated 17 of the 21 pairs." |
| P17 | 8 (SI 1416–1427), whole paragraph | "Figure S9 gives the Nemenyi critical-difference diagrams … Under rook contiguity, β is not identifiable for Checkerboard, and its estimates are kept and flagged." | (paragraph omitted; the non-separated pairs are stated in the text above it and in Section 10) |
| P18 | 9 (SI 1532) | "Figure S11 shows coverage across β for each design." | (sentence omitted; the paragraph after it, about Figure S11's bands, is already unabridged-only) |
| P19 | 10 (SI 1580–1583) | "…Checkerboard (0.480) (Table S7; Figure S12), and the same order held within each of the five incidence configurations (Table S8)." | "…Checkerboard (0.480) (Table S7), and the same order held within each of the five incidence configurations." |
| P20 | 10 (SI 1585) | "…whose coverage was 0.162 (Figure S13)." | "…whose coverage was 0.162." |
| P21 | 10 (SI 1612) | "…accounted for 0.438 of its MSE of 0.480 (Figure S14 shows the decomposition for the Poisson configuration with ρ_X = 0.20)." | "…accounted for 0.438 of its MSE of 0.480." |
| P22 | 10 (SI 1623) | "…under rook contiguity it was between 0.146 and 0.171 (Table S13)." | "…under rook contiguity it was between 0.146 and 0.171." |

## Batch C: robustness, SRS benchmark and case study

| ID | Where | Current | Curated |
|---|---|---|---|
| P23 | 12.1 (SI 1779) | "Table S15 gives the distribution of MSE across the 160 queen blocks at β = 1." | "We also examined the distribution of MSE across the 160 queen blocks at β = 1." |
| P24 | 12.1 (SI 1824) | "…and it had the lowest worst-case block MSE (Table S15)." | "…and it had the lowest worst-case block MSE." |
| P25 | 12.1 (SI 1828) | "The scenario-level worst cases in Table S15 summarize variation across parameter settings…" | "These scenario-level worst cases summarize variation across parameter settings…" |
| P26 | 12.2 (SI 1839) | "Table S16 separates the regimes and includes all eight strategies." | "Table~2 of the article separates the regimes for the six designs and SRS, and Figure S15 adds the two variants." |
| P27 | 13 (SI 1911) | "Figure S16 compares the eight designs with SRS by year and regime. Setting-specific comparisons are in Figure S18." | "Figure~4 of the article compares the six designs with SRS by year and regime." |
| P28 | 13 (SI 1915) | "…by 21.1–29.9% across the four observed years (Table S17)." | "…by 21.1–29.9% across the four observed years (Table S18)." (E8 folds S17 into S18, in both versions; the CTJ's "SI Table S17" pointer moves to S18 then.) |
| P29 | 13 (SI 1933) | "Table S5 and Figure S17 give the complete budgets and shares." | "Table S5 gives the complete budgets and shares." |
| P30 | 13 (SI 1984) | "Table S18 and Figure S19 give MSE, bias and coverage for every design in every year, and Table S5 gives the budgets." | "Table S18 gives MSE, bias and coverage for every design in every year, and Table S5 gives the budgets." |
| P31 | 14 (SI 2073) | "Figure S22 shows the ρ=0.5 setting." | (sentence omitted) |
| P32 | 14 (SI 2113) | "Figure S21 shows the full allocation-to-allocation distribution for the six designs and SRS. For each allocation in Figure S21, MSE is estimated…" | "Figure~3 of the article shows the full allocation-to-allocation distribution for the six designs and SRS. For each allocation in that figure, MSE is estimated…" |
| P33 | 15.1 (SI 2189) | "Power under control-only spillover rose slowly beyond about 40 colleges (Figure S24). In Figure S24 the two saturation designs nearly coincide, and the year-to-year range is narrower than the lines except near 58 colleges." | "Power under control-only spillover rose slowly beyond about 40 colleges." |
| P34 | 15.1 (SI 2205) | "…and 0.51 and 0.64 with Incidence-Guided Saturation (Table S22)." | "…and 0.51 and 0.64 with Incidence-Guided Saturation." |
| P35 | 15.4 (SI 2356) | "Tables S26 and S27 give each strategy's MSE relative to SRS under each model." | "Table S26 gives each strategy's MSE relative to SRS under each combination of data-generating and fitted model." |

Not listed: the four text blocks that are unabridged-only in their entirety (Section 4.3 apart
from Table S3, the paragraph on Figure S11's bands, the rook-figures paragraph, and 15.5) and the
29 floats themselves. Those carry no curated text, so they need no wording. The E3 move of the
error-decomposition section and the E8 merges come with the implementation, the merges as their
own sign-off list.

## Decisions (author, 2026-11-25)

- The curated SI and the CTJ are companions; the curated SI is self-contained and never refers to
  the unabridged version, which is treated as if it did not exist (material is cut instead).
- X1 accepted (the curated SI cites the article's exhibits with live numbers). X2 rejected (no
  availability sentence). X3 accepted, worded as a collection of everything prepared for the study.
- P1–P35 accepted. Two more pure pointer sentences found at build time are omitted from the curated
  text under the same rule: P36 "Figure S13 corresponds to Figure S7." and P37 "Figure S14
  corresponds to Figure S8." (rook section).
- E8: Table S17 is a subset of S18, so it becomes unabridged-only and the CTJ's pointer moves to S18.
