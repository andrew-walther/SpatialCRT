# Phase 2d sign-off: SI tables brought up to date with Chapter 3

Author: Andrew Walther. Drafted 2026-11-25. `verify_manuscripts.py` checks every Chapter 3
appendix table body (`Thesis_Chapter.qmd`) against the SI table with the same label. A
word-level comparison of all chapter tables against the SI (after 2c) finds only the items
below. Every chapter cell here was approved in an earlier round (the 2026-11-15 review and the
round 6 carry-over), so 2d copies the chapter's text into the SI. Captions are not compared and
are not changed here. Reply by ID.

## Tables whose cells differ (punctuation only; no number changes)

### TD1. `tab:designs` (Table S1): four rows take the chapter's sentences in place of semicolons

| Row | Current (SI) | Proposed (= Chapter 3) |
|---|---|---|
| Incidence-Guided Saturation | "…from highest to lowest mean; within each quadrant, that share of its 25 clusters is treated at random" | "…from highest to lowest mean. Within each quadrant, that share of its 25 clusters is treated at random" |
| Isolation Buffer | "…until none remain (a greedy random maximal independent set); no two treated clusters are neighbors" | "…until none remain. The result is a greedy random maximal independent set, in which no two treated clusters are neighbors" |
| 2×2 Blocking | "…25 non-overlapping 2×2 blocks; two of the four clusters in each block are treated at random" | "…25 non-overlapping 2×2 blocks, and two of the four clusters in each block are treated at random" |
| Checkerboard | "…are treated, so treatment alternates from cluster to cluster; under rook contiguity every neighbor…" | "…are treated, which makes treatment alternate from cluster to cluster. Under rook contiguity every neighbor…" |

### TD2. `tab:parameter-grid` (Table S3): two rows

| Row | Current (SI) | Proposed (= Chapter 3) |
|---|---|---|
| Incidence configuration | "iid; spatially autocorrelated, ρ_X ∈ {0.20, 0.50}; Poisson, ρ_X ∈ {0.20, 0.50}" | "iid, spatially autocorrelated with ρ_X ∈ {0.20, 0.50}, and Poisson with ρ_X ∈ {0.20, 0.50}" |
| Poisson mechanism | "Base rate 35 per 100,000; P = 100,000 per cluster" | "Base rate 35 per 100,000, and P = 100,000 per cluster" |

### TD3. `tab:excluded-designs` (Table S2): two rows

| Row | Current (SI) | Proposed (= Chapter 3) |
|---|---|---|
| Saturation Quadrants | "…without reference to incidence; within each quadrant, that share…" | "…without reference to incidence. Within each quadrant, that share…" |
| Balanced Halves | "…into two strata of 50; 25 clusters are treated at random in each stratum" | "…into two strata of 50, and 25 clusters are treated at random in each stratum" |

### TD4. `tab:nc-annual` (Table S18): minus sign

The bias −0.007 is typed `-0.007` (prints a hyphen) in four cells; the chapter has `$-0.007$`
(true minus). Proposed: `$-0.007$` in all four cells, as in the chapter.

## Tables the SI lacks

### TD5. `tab:psi`: rename the SI label

The SI already has this table (Table S10, MSE by spillover magnitude ψ) under the old label
`tab:gamma`, with a body identical to the chapter's. Proposed: rename `tab:gamma` → `tab:psi`
in its `\label` and its one `\ref` (SI 1160). Nothing printed changes.

### D2. `tab:adjacent-pairs` (Decision)

The chapter reports the queen 50-unit adjacent-pair tests as a table of the six designs; the SI
reports them in prose, which after 2c T2 also places SRS in the chain (as the "SRS in every
design test" rule requires). Copying the chapter table now would bring a six-design table
without SRS into the SI.

Recommendation: **keep the SI prose and mark `tab:adjacent-pairs` as pending in the verifier**
(a default entry in `PENDING_SI`, printed on every run, with a comment pointing to 2i). In 2i,
the chapter table gains the two SRS rows (Incidence-Guided Saturation against SRS; SRS against
Balanced Quartiles), and you decide then whether the SI also carries the table. Reason: no SI
content goes backward on the SRS rule, and the verifier still reports the gap loudly rather than
passing silently.

Alternative: add the table to the SI now with the two SRS rows. The verifier would then fail on
this table until 2i updates the chapter, unless the verifier also learns an "SI ahead of
Chapter 3" exception.

## After applying

Rebuild both SI versions; `python3 paper/tools/verify_manuscripts.py` must pass (with the D2
pending notice if D2 is accepted); run `test_manuscript_exhibits.R` and the word count. The CTJ
does not use these tables, so it needs no change. Chapter 3 needs nothing from TD1–TD5; D2 adds
the SRS rows to the chapter table in 2i.

## Decisions

- Author, 2026-11-25: TD1–TD5 accepted; D2 accepted as recommended. Applied with `apply_2d.py`
  (rows taken from `Thesis_Chapter.qmd`) and a default `PENDING_SI` entry for
  `tab:adjacent-pairs` in `verify_manuscripts.py`. The verifier now also resolves the CTJ's
  `\ref{SI-...}` against the SI's labels (a later check that the 2d drift had hidden). Result:
  `verify_manuscripts.py` passes, with the pending notice.
