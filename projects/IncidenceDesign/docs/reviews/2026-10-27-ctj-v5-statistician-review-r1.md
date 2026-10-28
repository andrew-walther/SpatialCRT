# CTJ v5 statistician review, round 1 (2026-10-27)

Read-only review of `CTJ_Manuscript.tex` and `CTJ_Exhibits.tex` (line numbers refer to `CTJ_Manuscript.tex`). Placeholder values (`\rerun`, `\hl`) were ignored. Numbers checked against `results/srs_benchmark/*.csv`, `results/manuscript_exhibit_revision_20261002/grid_configuration_means.csv`, `results/grid_tdf_2026-10/interval_comparison.csv`, `application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv`, `application/results/methods_ext_2026-10/stored/*coverage*.csv`, the chapter and the SI.

**Verified correct:** all of Table 2 (MSE, coverage and paired ratios with CIs); the pooled 0.0907/0.1285; the <10% configuration effect (HIF's max/min ratio is 1.19, all others are 1.10 or less); the pooled six-design order, identical in all five configurations; 21.1–29.9%; 0.30761/0.07638; 0.07428–0.07682; HIF 0.70–1.40; Checkerboard ~3x; IB 3.74; budgets (27, 15.89, 27–31) and population shares (22–24%, 27%); 20.7–43.4%; 35–37 colleges; 22%→5%; HIF 3.68–4.44 and 1.42–2.52; 23,523 deaths and 25,594,321 person-years; 1,248 settings.

## BLOCKING

**B1. Application coverage ranges include two SI-only variants.** Abstract Results and Results l.222 ("94.1--94.9\%", "92.6--93.6\%"). The lower bounds come from Saturation Regions, a variant reported only in the SI (t + df 0.9412, Wald 0.9263). For the six designs plus SRS, the design-by-regime means are 0.9421–0.9490 (t + df) and 0.9265–0.9362 (Wald), which round to 94.2–94.9% and 92.7–93.6% (`interval_comparison.csv`, App columns). *Fix:* use 94.2–94.9% and 92.7–93.6%, or say "across all eight designs and SRS". In l.222, also say what the range spans: "Across designs and regimes (averaged over years)".

**B2. "no design was more accurate than SRS" under both-arms spillover.** Abstract Results and l.204. Table 2 shows Balanced Quartiles at 0.0421 against 0.0424 for SRS (ratio 0.99). The claim holds only in the detectable sense. *Fix:* replace it with "no design was detectably more accurate than SRS" in both places (+1 word each).

**B3. "incidence guidance directs the most treatment to the highest-burden regions"** (abstract Conclusions; Conclusions l.254). High Incidence Focus directs more treatment to high-incidence clusters than Incidence-Guided Saturation does. That is its whole rule. The intended claim is narrower. *Fix:* "incidence guidance gives the highest-burden region the largest treated share (80%)". This is the same length and also fits the AGENTS framing rule.

## SHOULD

**S1. Design MSE as "the average of $m(a)$ over its possible assignments"** (l.169). This is wrong for non-uniform designs: Isolation Buffer's greedy draw does not give all maximal independent sets equal probability. *Fix:* "...is the expectation of $m(a)$ over the design's assignment distribution."

**S2. q90 is defined (l.169) but never reported** in the main text or in its exhibits. *Fix:* delete "by its 90th percentile (q90) and" (saves about 8 words).

**S3. "Wald interval" is never defined,** yet the abstract, l.222 and l.248 use it. *Fix:* at l.167, after the t interval, add "rather than the Wald interval $\hat\tau\pm1.96\,\widehat{\mathrm{SE}}$". Pay for the words with S2.

**S4. Power paragraph (l.231).** "An effect of one standard deviation" never says which SD (it means τ = 1 = σ, the cluster-level residual SD), and the test of τ = 0 is never stated. *Fix:* "at $\tau=1$, one residual standard deviation". At l.200, write "using the same $t$ distribution to test $\tau=0$".

**S5. "the approximate MCSE-based interval for each setting's difference excluded zero"** (l.222). This interval is not defined in the Methods. *Fix:* "...and each difference exceeded twice its MCSE", or define it once at l.200.

**S6. The budget explanation does not match the rules** (l.194, l.233). The sentence "Halving an odd-sized group requires rounding, so treatment budgets differed" explains only spatial blocking. Checkerboard's 27 comes from the graph rule and Isolation Buffer's 15.89 from the independent set. The rounding rule is never stated, yet l.233 relies on it ("always assigned to control"), and 31 treated in 2019 (more than 58/2) rules out simple rounding down. *Fix:* state the rule at l.194 (for example, "odd groups rounded at random, singletons to control", if that is correct). Then: "Checkerboard, Isolation Buffer and spatial blocking therefore treated different numbers of clusters; we report these budgets..."

**S7. Figure 3 caption: "Settings are described in the Methods."** The Methods give two different allocation-risk settings: grid ρ ∈ {0, 0.5}, γ ∈ {0.5, 0.8}, and NC ρ ∈ {0, 0.5}, γ = 0.8 in four years. *Fix:* "Top row, grid (Methods settings); bottom row, application at $\rho\in\{0,0.5\}$, $\gamma=0.8$, four years." This settles the author's concern about telling Figure 2 and Figure 3 apart.

**S8. A range does not say what it spans** (l.215: "pooled bias $-0.12$ to $-0.40$ and coverage 0.54 to 0.91"). The chapter (l.311) shows this range covers designs and both neighbor types, with rook Checkerboard excluded. *Fix:* "(bias $-0.12$ to $-0.40$ and coverage 0.54 to 0.91 across designs and neighbor types at $\tau=1$)".

**S9. "the low both-arms error of High Incidence Focus"** (l.241). This is true only in the application. On the grid, HIF has 3.82 times SRS's both-arms MSE. *Fix:* "its low both-arms error in the application".

**S10. Conclusions overstate** (l.254 "saturation designs gave the lowest average error and favorable worst-case error"; abstract "saturation designs had the lowest estimation error"). The main text names only one saturation design. Isolation Buffer is lower at τ = 3 (l.213), HIF beats IGS in NC 2018 (ratio 0.697 vs 0.789), and the figures show worst-decile error, not worst-case error. *Fix:* "incidence-guided saturation had the lowest or near-lowest average error and lower worst-decile error than SRS under control-only spillover".

**S11. Design names are inconsistent.** Prose writes "incidence-guided saturation" in lower case (l.55, 204, 217, 222...). The table, the figures and l.213 write "Incidence-Guided Saturation", and the other designs are capitalised in prose. The standing rule requires one name set. *Fix:* capitalise throughout. "Spatial blocking" for the NC version of 2×2 Blocking is acceptable if Table 1 says so once. It currently says only "Treat a random half...".

**S12. Rerun-dependent Methods** (l.184, l.217). The "400 outcomes for selected assignments" selection rule and the worst-decile statement are placeholders. The chapter (l.427, 429, 586, 621, 676) and the SI still describe the two-surface "grid pilot". After the rerun, sync both and define the selection rule in one clause.

**S13. "differences too small to establish either equivalence or superiority"** (l.222). No equivalence margin was set. *Fix:* "differences within simulation error".

**S14. "All 1,248 application settings"** (l.220). The Methods never add up to this count. *Fix:* "Every application setting, including sensitivity analyses, completed..." This is shorter and needs no total.

## NIT

- N1 l.143: "Rather, adjacent populations..." "Rather" follows no negation. Use "In practice, adjacent...".
- N2 l.229: "Allocation risk in the application told the same story" is a pivot cliché. Use "In the application, under control-only spillover at...". Also l.211: "Isolation Buffer behaves differently." can merge into the next sentence ("For Isolation Buffer, no treated cluster...").
- N3 Colons: l.239 ("conditional guidance: consider saturation..."). Use "support considering saturation when..., and retaining SRS...". Section headings, see concern table.
- N4 l.213: 0.1822 at τ = 3 does not appear in the chapter or SI, which give 0.182. Confirm the four-digit value, and confirm IB at τ = 3 is 0.1596.
- N5 Figure 2: in the both-arms panel most points lie at 0.04–0.06, but the lowest axis tick is 0.1. Add a 0.03 tick.
- N6 Table 2 caption: say the interval is a $t_{49}$ interval for the mean paired difference divided by SRS MSE, with no multiplicity adjustment.
- N7 l.147: "keeps clusters of the same arm apart" is not true under queen contiguity, where same-arm clusters share corners. Use "so no two edge-sharing clusters share an arm".
- N8 Title page word counts (3489/411) should be recounted after these edits.

## Reproducibility

The model, estimator, interval, metrics, MCSE units, surfaces, replication counts and application inputs are adequate with the SI. Three items are missing from the main text: the spatial-blocking rounding rule (S6), the 400-outcome selection rule (S12) and the power test (S4). The participation analysis's treatment of non-participating colleges sits in SI S11 only, which is acceptable.

## Author concerns (v5_author_notes.md)

| Concern | Status | Evidence |
|---|---|---|
| Abstract: spillover framing | RESOLVED | l.51 "spill over and contaminate clusters that neighbor treated clusters" |
| Abstract: two criteria | RESOLVED | l.51 average error and allocation risk |
| Abstract: incidence drives allocation, not outcome | RESOLVED | l.53 "outcome that does not depend on incidence" |
| Abstract: "simulated education outcome" unclear | PARTLY | still does not say who is measured (healthcare professionals) |
| Abstract: last Methods sentence | RESOLVED | split into queen primary plus three sensitivity verbs |
| Abstract: τ = 1 | RESOLVED | l.55 |
| Abstract: precision-target sentence | RESOLVED | removed |
| Abstract: conclusions (budget, jargon) | PARTLY | budget/share and nominal coverage added; overstatement B3/S10 |
| Intro: colons | RESOLVED | none left; "Rather," awkward (N1) |
| Intro: "larger study region" | RESOLVED | l.147 |
| Intro: explain alternating allocation | RESOLVED | l.147 (N7 precision) |
| Intro: identifiability | RESOLVED | l.145 exposure "fully determined by its own treatment" |
| Intro: stuffy, long application paragraph | PARTLY | l.149 still four long sentences |
| Methods: model first | RESOLVED | §Outcome model precedes designs and simulation |
| Methods: cite queen/rook; primary/sensitivity | RESOLVED | l.153 |
| Methods: why β = 1 | RESOLVED | l.165 |
| Methods: cite 35/100,000 | RESOLVED | l.180 Mirzaei |
| Methods: ten surfaces | PARTLY | stated as shared, but what a surface is and why ten are not said |
| Methods: Monte Carlo error | RESOLVED | l.182 with Morris |
| Methods: choppy / colons | PARTLY | improved; see N2, N3 |
| Methods: define allocation-risk terms | PARTLY | terms defined; S1, S2 |
| Methods: "pilot" language | RESOLVED | only the planned education pilot remains |
| Methods: no "Because" | RESOLVED | none found |
| Methods: application heading colon | NOT RESOLVED | l.186 "Application: Sudden..."; Results headings l.203, l.219 also |
| Methods: clarify rules | PARTLY | precision rule stated; rounding rule missing (S6) |
| Methods: what replication extension was | PARTLY | l.200 clear; grid l.184 placeholder (S12) |
| Results: choppy, define terms | PARTLY | S3, S5 |
| Results: precision targets named | RESOLVED | l.220 |
| Results: say what ranges span | PARTLY | B1, S8 |
| Results: run totals | RESOLVED | none in Results (S14) |
| Results: "preserves regime distinction" / transfer | RESOLVED | l.227 rewritten |
| Results: "grid pilot" | RESOLVED | absent (chapter still has it, S12) |
| Results: "Budgets differed" fragment | RESOLVED | l.233 full sentence |
| Results: key results only | PARTLY | power, omitted-spillover and τ = 2/3 detail could go to SI |
| Discussion: plain language | PARTLY | "identified", Wald undefined |
| Discussion: education endpoint and resources | RESOLVED | l.251 |
| Exhibits: Figure 2 vs 3 distinct | PARTLY | strips label rows; caption vague (S7) |
| Exhibits: Figure 2 CIs | RESOLVED | intervals drawn; caption says so |
| Exhibits: Figure 3 caption length | RESOLVED | five sentences |
| Exhibits: Table 2 SRS 50 of 100 | RESOLVED | caption |
| Overall: flow, reproducibility | PARTLY | see S6, S12, N2 |
