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
