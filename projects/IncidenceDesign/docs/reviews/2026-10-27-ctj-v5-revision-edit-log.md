# CTJ V5 revision: edit log

Created 2026-10-27. Source review: Andrew's notes (pasted in session) and
`paper/ctj_manuscript/manuscript_ctj/revision_markup/CTJ_Reading_Markup_V5.pdf`
(100 annotations). Decisions: `DECISIONS.md` (2026-10-27). Files changed:
`CTJ_Manuscript.tex`, `CTJ_Exhibits.tex`, `figures/Figure_Legends.md`,
`code/19_manuscript_exhibits.R`, `code/20_manuscript_figure_revision.R`,
`paper/tools/verify_manuscripts.py`, `paper/SpatialCRT_IncidenceDesign.bib`,
`Supplementary_Information.tex`, `Thesis_Chapter.qmd`.

## Word counts

| | Before | After |
|---|---:|---:|
| Body (Introduction to Discussion, texcount text + headers) | 3,489 | 3,499 |
| Abstract | 411 (410 by the same texcount method) | 404 |

The full IRB statement and its placeholder moved to an "Ethical considerations"
declaration; a short IRB sentence stays in the application Methods, as Clinical Trials
requires.

## V5 notes and markup, one row each

Section labels refer to the revised manuscript.

### Front matter, abstract, keywords

| Note | Action | Where |
|---|---|---|
| "treatment effects can contaminate clusters that neighbor treated clusters" | Used: "can spill over and contaminate clusters that neighbor treated clusters" | Abstract, Background |
| Criteria (MSE; error under unfavorable allocation) | Kept, with the colon replaced by a comma ("on two criteria, the average…") | Abstract, Background |
| "simulated continuous education outcomes" unclear | "simulated a continuous education outcome that does not depend on incidence" | Abstract, Methods |
| Last Methods sentence unclear | Two sentences: queen contiguity for primary; sensitivity analyses named as actions (rook, added baseline, alternative region ordering) | Abstract, Methods |
| "treatment effect of one" → τ = 1 | "with treatment effect τ = 1" | Abstract, Results |
| Precision-target sentence unclear | Cut from the abstract (Results text now names the targets) | Abstract, Results |
| Conclusions: "designs minimize estimation error", budget/proportion treated, "interval calibration" | Rewritten: saturation "had the lowest estimation error"; budget = number treated and population share; "whether confidence intervals reach their nominal coverage" | Abstract, Conclusions |
| Keywords | Unchanged (resolved in V5) | — |

### Introduction

| Note | Action | Where |
|---|---|---|
| Subsections? | No (decided); no outline paragraph in CTJ; outline added to Chapter 3 instead | — |
| Colon after "own treatment" → "rather, adjacent populations…" | Two sentences, second opens "Rather," | P1 |
| "maximal interspersion…" = identifiability issue? | Yes; added a plain example (strict checkerboard with edge neighbors: exposure fully determined by own treatment) | P2 |
| "larger grid" → "larger study region" | Done | P3 |
| Strict alternating allocation → BSS / checkerboard | "block stratified sampling, an alternating allocation that keeps clusters of the same arm apart"; "the alternating allocation, here called Checkerboard" (not claimed identical) | P3 |
| Add another comma | Serial comma "additional strategies, and an irregular…" | P3 |
| Strike "consequently", "therefore" | Both removed | P4 |
| Application intro length | Two paragraphs merged and trimmed (application motivation + incidence strategies) | P4 |

### Methods

| Note | Action | Where |
|---|---|---|
| Model before simulation; outline (model → designs → simulation → application) | New order: Outcome model, estimation and allocation risk; Treatment assignment designs; Simulation study; Application | Methods |
| "on a 100-cluster grid" not needed | Heading removed with the restructure | — |
| Write model name; "Outcomes followed the SAR model" | "Outcomes followed the spatial autoregressive (SAR) lag model" with every symbol glossed after Eq. 1 | Model |
| Cite queen/rook; "queen was primary" wording | Defined with citation (`bivand_applied_2008`); "Queen contiguity was used for the primary analysis and rook contiguity for a sensitivity analysis" | Model |
| Why β = 1? | Prognostic-incidence rationale (harder case; investigators choose β); application β = 0 stated in the application section | Model; Application |
| "what is this???" (SAR filter + normal CDF) | Plain description: smooth normal values across queen neighbors, map to (0, 1), so nearby clusters have similar incidence | Simulation study |
| Cite base rate 35/100,000 | Cited Mirzaei et al. 2019 ("between published US SUD rates for working-age men and women") | Simulation study |
| "ten independent surfaces"? | "ten independent incidence surfaces" per configuration, shared by every setting and design | Simulation study |
| Mention Table 1 earlier | Table 1 is the first sentence of "Treatment assignment designs" | Designs |
| "…;" → "and"; "proportions of treated clusters"; "ranked incidence groups" | "Balanced Quartiles spreads treatment evenly across ranked incidence groups, whereas…"; "the same local proportion of treated clusters" | Designs P2 |
| Don't start a sentence with "Because" | Removed; no sentence in the manuscript opens with "Because" | Simulation study |
| Explain/cite Monte Carlo error | Defined; cited Morris, White and Crowther 2019 (verified in PubMed, DOI 10.1002/sim.8086; added to bib) | Simulation study |
| "concentrated likelihood" mentioned anywhere? | Replaced by plain wording: "Our own ML routine reproduced the estimates of `spatialreg::lagsarlm` to within 10⁻⁶ in validation" | Model |
| Too many short sentences ("This assumes… A grid sensitivity omitted spillover.") | Joined: "This fit assumes that the analyst knows the form of spillover, so a grid sensitivity analysis refit the model without the spillover term." | Model |
| Interval was written as Wald (τ̂ ± 1.96 SE) | Corrected to the reported t + df interval, with its formula and symbols defined | Model |
| Define the allocation-risk equation | Every term defined (assignment a, fixed incidence X, expectation over outcome noise ε); link to design MSE | Model P3 |
| "pilot" wording; "Its findings support pilot evidence…" | Removed; grid allocation risk rerun on all 10 surfaces (see DECISIONS 2026-10-27) | Simulation study P3 |
| Repeated-outcomes paragraph rewrite (use "because"/"so", not a colon) | Rewritten into the allocation-risk paragraph | Model P3 |
| Heading "Application – Sudden unexpected death in North Carolina" | "Application: Sudden unexpected death in North Carolina" (subsection-title colon kept, as in the Results headings you specified) | Application |
| Drop colons ("presumed non-sudden:") | Rewritten as "…excluded as presumed non-sudden, giving 23,523 deaths…" | Application P1 |
| SEER source | From Habib's data file; cited SEER population data (new bib entry, URL verified) | Application P1 |
| 2×2 Blocking sentence clean-up ("where" instead of colon) | "…became spatial blocking, where 14 groups of nearby clusters were formed…, and a random half of each group was treated." | Application P3 |
| "0.01; these" | "0.01. These MCSEs…" | Application P6 |
| (reviewer) power test | Settings stated (ρ = 0.5, γ = 0.8τ); two-sided rejection ("interval excluded zero") | Application P6 |
| What is the "extension"? | Removed as a concept: "For allocation risk… we simulated 400 outcomes for each of the same assignments, so that each assignment's risk was estimated precisely enough to rank." | Application P6 |
| Rules/policies should be clarified (rounding) | "Halving an odd-sized group requires rounding, so treatment budgets differed slightly…" | Application P3 |
| Choppy short sentences / colons and semicolons (general) | Methods rewritten throughout; no prose semicolons remain; colons only in subsection titles and the Discussion sentence you kept | Methods |

### Results

| Note | Action | Where |
|---|---|---|
| Headings | Your two headings, verbatim | Results |
| "by under 10%" | "by less than 10%" | Simulation P1 |
| "against 0.1285 for SRS" | Done | Simulation P1 |
| Drop "separately" | Sentence rewritten ("Rook results and further breakdowns are in the SI.") | Simulation P1 |
| "poorer" → "worse" | Done | Simulation P2 |
| "lowest; SRS" | Two sentences | Simulation P3 |
| Meaning of the omitted-spillover sentence | Rewritten to state the result and its point (lower MSE yet invalid inference) | Simulation P4 |
| Figure 3 not on full data? / "grid pilot" | Grid allocation risk rerun on all ten surfaces; "pilot" removed | Simulation P5; Figure 3 |
| Precision targets? | Named where they are used: MCSE of mean MSE ≤ 5% of its value; MCSE of coverage ≤ 0.01 | Application P1 |
| Run totals (7,598,000 …) | Cut from the main text (in SI S3) | — |
| "2,000 outcomes" sentence | Cut from the main text (in SI) | — |
| Drop "descriptive" | Done | Application P2 |
| "preserves the broad regime distinction" | "The NC map kept the regime distinction but changed several grid results." | Application P3 |
| "Spatial blocking is close…" | Done | Application P3 |
| Transfer sentence rewrite | "A design's ranking on the grid therefore need not carry over unchanged to an irregular geography." | Application P3 |
| Trim "(detailed displays in the supplementary material)" and abbreviate SI | "supplementary information (SI)" defined once, then "SI" | Throughout |
| 61–86% range: what is it? | Range over the eight year-by-ρ settings; moved to SI under the word limit, summarized as "most of the variance in estimated risk came from Monte Carlo error" | Application P4 |
| "saturation designs' advantage" | Incidence-guided saturation's own numbers (22% at 58 colleges → 5% at 29) | Application P5 |
| Budgets non-sentence | "Treatment budgets differed by design and year (SI)." Budget detail kept | Application P6 |
| MSE ranges = intervals? | No; each range now says what it spans (years, designs and regimes, year-by-ρ settings) | Throughout |
| "but its grid-pilot both-arms mean was 3.91 times…" | Moved to SI with the other designs' tail ratios (word limit) | — |

### Discussion

| Note | Action | Where |
|---|---|---|
| Keep the implications colon | Kept | Implications P1 |
| Plainer limitations | Rewritten in plain terms (wrong neighbor network, varying effects, samples of individuals, equal cost, fixed incidence) | Limitations |
| Heading "Continued work for the SUD education pilot" | Done | — |
| Define education endpoint and staffing resources | Listed among what investigators must define | Continued work |
| Conclusions without jargon | Rewritten | Conclusions |

### Exhibits

| Note | Action |
|---|---|
| Table 1 caption "representative incidence data" | "a representative spatial surface (ρ_X = 0.5) drawn for illustration" |
| Table 1 incidence row; Balanced Quartiles row | Your wording |
| Figure 1 "planning" → "incidence"; panel A title | Caption "Observed North Carolina sudden unexpected death (SUD) incidence"; panel A "Observed 2018-2021 SUD rates (pooled)" |
| Figure 1 panel C (yearly maps) | Added on a shared color scale; legible in both layouts |
| Table 2 caption wording; SRS "50 clusters" | "Grid benchmark at τ = 1 with queen contiguity"; "spatial lag model that includes the true spillover term"; "randomly selects 50 of the 100 clusters for treatment" |
| Table 2 paired comparison (decided) | New "Ratio to SRS (95% CI)" column per regime (t₄₉ over 50 surface units; `results/srs_benchmark/manuscript_paired_ratios.csv`) |
| Table 2 sorted? | Rows follow Table 1's design order (not sorted by MSE) |
| Figure 2 shapes, caption, confidence intervals | Five filled shapes in distinct colors; 95% Monte Carlo intervals over ten surfaces (t₉); your caption wording; points to Table 2 for comparisons |
| Figure 3 vs Figure 2; shorter caption | Caption shortened and opens with "Allocation risk"; Figure 2 = average MSE, Figure 3 = assignment-to-assignment variation relative to SRS |
| Figure 4 caption; axis label | "Completed primary application" dropped; caption shortened; axis "Year (observed SUD incidence)" |

### Admin

| Note | Action |
|---|---|
| Funding | Sage standard sentence; confirm with Habib and Simpson after the talk |
| Supplementary material wording | "For the sudden unexpected death application…"; "supplementary information" |
| IRB statement | Moved to "Ethical considerations"; placeholder kept |
| DOIs | Kept |
| Author contributions placeholder | Kept (coauthor approval after the talk) |

## Not changed, and why

- Title, authors, keywords: approved in V5.
- Why community-college service areas are the clusters: reviewers asked; left for
  Habib/Simpson (TODO).

## Grid allocation-risk rerun

All ten surfaces, same corner settings (`docs/findings/2026-10-27-allocation-risk-surfaces10-findings.md`).
The pilot reproduced exactly on surfaces 1–2 and no conclusion changed. New CTJ values:
worst-decile mean 0.1912 vs 0.3462 (control-only) and 0.0583 vs 0.0544 (both arms). The
400-outcome precision step covered one setting only, so the CTJ no longer mentions it. The
SI and Chapter 3 allocation-risk text, `tab:downside` grid columns and the Figure 3 caption
were updated to match.

## Independent review rounds

Two fresh read-only reviewers per round (statistician; clinician/epidemiologist), each
given the manuscript and the author's V5 notes but not this log.

| Round | Statistician (B/S/N) | Clinician (B/S/N) | Main fixes |
|---|---|---|---|
| 1 | 3/14/8 | 2/17/17 | Coverage ranges limited to the six designs and SRS; "detectably" more accurate; "largest treated share" instead of "most treatment"; design MSE as an expectation; Wald, identification, τ units defined; rounding rule stated; plain spillover-regime meaning; design-name capitalization |
| 2 | 0/10/9 | 1/13/9 | Regimes in plain words, with "either could apply" (author); power-study settings stated; "assignment-specific MSE" vs "allocation risk"; why MSE grows with τ; Figure 3 caption settings |
| 3 | 0/3/10 | 0/5/10 | Power gap restated in Continued work; split-half check worded as a check; Table 2 "one plus"; rounding wording |
| 4 | 0/1/3 | 0/4/4 | Grid tail noise share (63%); application power test is two-sided ("interval excluded zero"; code uses abs(t)); small wording fixes |

Reports: `2026-10-27-ctj-v5-{statistician,clinician}-review-r{1..4}.md`. Reviewer
suggestions not taken: removing pooled grid numbers (the author edited that sentence);
cutting the budget paragraph (the author asked to keep budgets); naming a specific
spillover example for the SUD pilot (the author asked to keep it general); the claim that
coverage was checked only under the SAR model (not verified, since the SI has CAR analyses).

## Final counts

Body 3,499 (limit 3,500); abstract 404 (limit 425); six exhibits. `verify_manuscripts.py`
and `code/tests/test_manuscript_exhibits.R` pass; both CTJ layouts and the SI build with no
undefined references.
