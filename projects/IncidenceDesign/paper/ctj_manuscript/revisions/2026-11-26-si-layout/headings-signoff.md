# SI Part II headings: narrative-style sign-off list (2026-11-26)

**What.** Proposed visible-text changes to the `\section` and `\subsection` headings of
Part II (Supplementary Results) of
`paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex`, in both builds
(curated and unabridged).

**Why.** You asked for SI headings that, where possible, state the result, as in the DeSurv
SI Appendix ("The DeSurv model shows consistent convergence and stability across
restarts"), and for plainly descriptive headings to stay where a finding headline would be
forced. Part I (Supplementary Methods) stays descriptive, as in DeSurv.

**How to apply.** Approve by ID (for example "accept NH-2, NH-3, NH-5; reject NH-7"). For each
accepted ID, replace the "Current" line with the "Proposed" line. Each "Current" line
matches the source exactly once (`grep -F -c` = 1, checked 2026-11-26). Only the visible text
changes. Every `\label{...}` stays as it is, so internal `\ref`s and the CTJ's
`\ref{SI-<label>}` pointers resolve as before. The CTJ cites SI sections only through
`\ref{SI-...}` (for example `SI-allocation-risk-results`, `SI-igs-sensitivity-results`,
`SI-interval-calibration`, `SI-nc-sensitivity-results`, `SI-app:technical`). Neither the CTJ
sources nor `paper/tools/*.py` quote any Part II heading text, so no section number and no
cross-reference changes. After applying, rebuild SI, CTJ, SI (Build gotchas in AGENTS.md)
and rebuild `Supplementary_Information_Full.tex`.

Section numbers below are from the curated `Supplementary_Information.pdf` (pdftotext,
2026-11-26). The unabridged `Supplementary_Information_Full.pdf` has the same numbers for
Sections 6 to 16, plus 16.5, which appears only in that build.

Evidence follows each section's own text, tables and figures. No heading adds a number the
section does not print. Only NH-7 (\(\beta\), \(\hat\rho\)) and NH-17 (\(t\)) use math,
both wrapped in `\texorpdfstring`. Design names are the canonical nine. No heading
names a single best design, and none says incidence guidance has no benefit.

---

## Recommendation

**Accept (clear, conservative, directly supported):** NH-2, NH-3, NH-4, NH-6, NH-8, NH-9, NH-10,
NH-11, NH-12, NH-14, NH-15, NH-16, NH-18, NH-19.

**Accept with a choice or a look first:** NH-1 (accurate but bland; keep the current heading
if you prefer), NH-5 (headline is the one pair the tests did not separate), NH-7 (long,
and technical for a TOC line), NH-13 (two candidate headlines; A recommended), NH-17
(needs the "main case study" qualifier, because the bootstrap analysis has every method,
including the corrected \(t\), below 94.0% somewhere).

**Keep as they are:** Sections 6, 12 and 16 (umbrella headings) and 16.5 (unabridged only).
Reasons are in the last list.

---

## Proposed changes

### NH-1 · Section 6.1 · both builds
- **Current:** `\subsection{Design performance}\label{design-performance-under-queen-contiguity}`
- **Proposed:** `\subsection{Designs differed widely in MSE but all kept near-nominal coverage}\label{design-performance-under-queen-contiguity}`
- **Reason:** Under queen contiguity, Table S5 (curated; `tab:performance`) gives MSE from 0.091 to 1.087 and coverage from 0.946 to 0.950 for all six designs. The heading avoids the pooled ranking, which hides the regime dependence of Section 6.3. *Optional; the current heading is also fine.*

### NH-2 · Section 6.2 · both builds
- **Current:** `\subsection{Design rankings by incidence configuration}\label{design-rankings-by-incidence-configuration}`
- **Proposed:** `\subsection{The incidence configuration did not change the design ranking}\label{design-rankings-by-incidence-configuration}`
- **Reason:** The text reports the pooled mean-MSE order within every one of the five configurations, with only slight changes in MSE and coverage of 0.941 to 0.952. Article Figure 2 is the curated evidence, and `tab:by-config` (unabridged) also supports it.

### NH-3 · Section 6.3 · both builds
- **Current:** `\subsection{Spatial dependence, spillover magnitude and regime}\label{sensitivity-to-spatial-dependence-spillover-magnitude-and-spillover-regime}`
- **Proposed:** `\subsection{The spillover regime changed the ranking, while spatial dependence and spillover magnitude did not}\label{sensitivity-to-spatial-dependence-spillover-magnitude-and-spillover-regime}`
- **Reason:** The section says that the spillover regime "was the one factor that changed the queen ranking" (article `tab:grid-srs`; `tab:regime`, unabridged). It also says that no design's MSE moved by more than 0.02 across \(\rho\) or (except Checkerboard, with no trend) across \(\psi\) (`tab:rho`, `tab:psi`, unabridged). Figure S3's variance finding is secondary, so the heading leaves it out.

### NH-4 · Section 6.4 · both builds
- **Current:** `\subsection{The eight-design comparison}\label{app:eight-design}`
- **Proposed:** `\subsection{Each related variant had nearly the same MSE as its counterpart}\label{app:eight-design}`
- **Reason:** The text says that each excluded design "is adjacent to its retained counterpart in mean MSE". Table S5 (curated) and `tab:eight-design` (unabridged) give queen MSE 0.090 against 0.091 for Saturation Quadrants and Incidence-Guided Saturation, and 0.130 against 0.131 for Balanced Halves and Balanced Quartiles, with similar rook values. "Related variant" is the term in Section 2.1's heading.

### NH-5 · Section 7 · both builds
- **Current:** `\section{Formal statistical comparisons}\label{formal-statistical-comparison}\label{app:supp-comparisons}`
- **Proposed:** `\section{Formal tests separated most pairs of designs but not Balanced Quartiles from SRS}\label{formal-statistical-comparison}\label{app:supp-comparisons}`
- **Reason:** Nemenyi comparisons separated 17 of 21 pairs and Holm-adjusted Wilcoxon tests 20 of 21, all except Balanced Quartiles against SRS (\(p = 0.15\)), and the adjacent-pair test agreed (text; `tab:adjacent-pairs` and Figure `fig:cd`, unabridged). This wording avoids a "ranked first" headline for one design. *Check:* the section reports queen contiguity only, and the heading does not say so. Under rook contiguity (Section 10, unabridged), the Nemenyi comparisons also failed to separate Incidence-Guided Saturation from those two.

### NH-6 · Section 8 · both builds
- **Current:** `\section{Stability across treatment effect sizes}\label{stability-across-treatment-effect-sizes}`
- **Proposed:** `\section{The design order changed little across treatment effect sizes}\label{stability-across-treatment-effect-sizes}`
- **Reason:** Figure S4 and the text show that the queen orders by average rank and by mean MSE were the same at every value of \(\beta\). Under rook contiguity, only the close pair of \(2\times2\) Blocking and Isolation Buffer traded places, at \(\beta = 3\). "Changed little" rather than "did not change" also allows for Isolation Buffer passing Incidence-Guided Saturation under control-only spillover at \(\beta = 3\).

### NH-7 · Section 9 · both builds
- **Current:** `\section{\texorpdfstring{Estimation error and its dependence on \(\beta\)}{Estimation error and its dependence on beta}}\label{app:technical}\label{app:technical-notes}\label{estimation-error-and-its-dependence-on-tau}`
- **Proposed:** `\section{\texorpdfstring{Estimation error can grow with \(\beta\) through the error in \(\hat\rho\)}{Estimation error can grow with beta through the error in the estimated rho}}\label{app:technical}\label{app:technical-notes}\label{estimation-error-and-its-dependence-on-tau}`
- **Reason:** Equation `eq:error-decomposition` and the text attribute almost all of Incidence-Guided Saturation's rise in bias (0.029 to 0.047) and control-only MSE to the term \(\beta(\rho-\hat\rho)g\). It is absent for Isolation Buffer (\(g \approx 0\)), hence "can grow". This follows the prelim Appendix B heading "Estimation Error Depends on the Treatment Effect". A shorter alternative that reuses that wording: `\section{\texorpdfstring{Estimation error depends on \(\beta\)}{Estimation error depends on beta}}...`.

### NH-8 · Section 10 · both builds
- **Current:** `\section{Rook contiguity}\label{rook-contiguity}\label{app:rook-figures}`
- **Proposed:** `\section{Under rook contiguity, every design except Checkerboard kept nominal coverage}\label{rook-contiguity}\label{app:rook-figures}`
- **Reason:** The text and Table S5 (rook rows) give coverage of 0.948 to 0.952 for every design except Checkerboard (0.162), whose coefficient estimates \(\beta-\psi\). The heading keeps "rook contiguity" so readers can still find the section in the TOC. The rook order differs from the queen order (\(2\times2\) Blocking moves to third), so a "same ranking" headline would be wrong.

### NH-9 · Section 11 · both builds
- **Current:** `\section{Non-oracle sensitivity analysis}\label{non-oracle-sensitivity-analysis}\label{app:oracle-nonoracle}`
- **Proposed:** `\section{Omitting the spillover term biased the estimates of every design}\label{non-oracle-sensitivity-analysis}\label{app:oracle-nonoracle}`
- **Reason:** Table S6 (`tab:oracle-nonoracle`) gives non-oracle bias from \(-0.121\) to \(-0.400\) and coverage from 0.54 to 0.91 outside rook Checkerboard, which is biased under both models (\(-0.652\)). The lower non-oracle MSE for some designs is a secondary point, and the text itself qualifies it.

### NH-10 · Section 12.1 · both builds
- **Current:** `\subsection{Robustness across blocks}\label{robustness-across-blocks-and-the-recommended-design}`
- **Proposed:** `\subsection{No single design was best under both spillover regimes}\label{robustness-across-blocks-and-the-recommended-design}`
- **Reason:** This is the section's own synthesis sentence ("No single design was best under both spillover regimes"). It rests on the block distributions (`tab:robustness`, unabridged; median and worst-case values in the text of both builds) and matches the standing rule against a single recommended design. A heading that names one design's worst case would cut against that rule. The label still says "recommended-design", but it stays unchanged because it is invisible.

### NH-11 · Section 12.2 · both builds
- **Current:** `\subsection{Comparison with SRS}\label{srs-benchmark-and-allocation-specific-risk}`
- **Proposed:** `\subsection{The saturation designs had lower MSE than SRS under control-only spillover and slightly higher under both-arms spillover}\label{srs-benchmark-and-allocation-specific-risk}`
- **Reason:** Figure S5 (`fig:vs-srs`) gives 36.9% lower control-only MSE for Incidence-Guided Saturation and 36.7% for Saturation Quadrants (all 50 units). It gives 8.3% and 4.3% higher both-arms MSE, both significant under queen contiguity. Under rook contiguity the both-arms difference was 3% to 5%. Two TOC lines.

### NH-12 · Section 12.3 · both builds
- **Current:** `\subsection{Allocation-specific risk}\label{allocation-risk-results}`
- **Proposed:** `\subsection{Under control-only spillover, the saturation designs' advantage over SRS held for their worst allocations}\label{allocation-risk-results}`
- **Reason:** For Incidence-Guided Saturation, the estimated 90th percentile (0.173 against 0.296) and the worst-decile mean (0.191 against 0.346) were below SRS in all 200 units. Saturation Quadrants "also performed favorably", and Table S8 (`tab:downside`) gives grid worst-decile ratios of 0.55 and 0.54. The both-arms tail differences are small, which is why the heading names only the control-only regime.

### NH-13 · Section 13 · both builds
- **Current:** `\section{Sensitivity of incidence guidance to regions and saturation spread}\label{igs-sensitivity-results}`
- **Proposed (A, recommended):** `\section{Wider saturation spreads lowered control-only MSE but could raise both-arms MSE}\label{igs-sensitivity-results}`
- **Alternative (B):** `\section{Incidence guidance changed where treatment went more than how accurate the estimate was}\label{igs-sensitivity-results}`
- **Reason:** Figure S6 (`fig:igs-sensitivity`) shows that widening from 60/40 to 90/10 lowered control-only MSE for every number of regions, while both-arms MSE rose for halves and quadrants (up to 1.94 times SRS at 90/10, and 8.13 at 100/0). That trade-off fills most of the section, so A is recommended. B is the section's own topic sentence for its third paragraph (guided and unguided MSE were close, and targeting rose from 52.7% to as much as 60.1%). It fits the "matches Saturation Quadrants' accuracy while targeting" framing, but it covers only one paragraph, and readers could take it as "no benefit" if it were shortened.

### NH-14 · Section 14 · both builds
- **Current:** `\section{Accuracy, inference and reach in the case study}\label{completed-results-accuracy-inference-and-reach}`
- **Proposed:** `\section{In the case study, both saturation designs had lower control-only MSE than SRS in every year}\label{completed-results-accuracy-inference-and-reach}`
- **Reason:** Table S7 (`tab:nc-annual`) gives control-only MSE of 0.223 for Saturation Quadrants and 0.216 to 0.243 for Incidence-Guided Saturation in each year, against 0.308 for SRS. The text adds that both were lower in every primary setting, with intervals that excluded zero. The coverage and reach findings (Figure S7) are secondary. The heading does not touch the targeting result, so the "small gain on the NC map" framing in the text is unaffected.

### NH-15 · Section 15 · both builds
- **Current:** `\section{Allocation downside in the case study}\label{allocation-downside-and-sensitivity}`
- **Proposed:** `\section{In the case study, the saturation designs also had lower worst-decile MSE than SRS under control-only spillover}\label{allocation-downside-and-sensitivity}`
- **Reason:** The NC control-only columns of Table S8 (`tab:downside`) give worst-decile ratios of 0.64 for Saturation Quadrants and 0.68 for Incidence-Guided Saturation, and the text gives 20.7% to 43.4% reductions. The text says the both-arms tail differences "are too small to establish better downside protection", so the heading is limited to control-only spillover. High Incidence Focus (one allocation, 0.64) is also below 1. The heading does not say "only".

### NH-16 · Section 16.1 · both builds
- **Current:** `\subsection{Colleges needed and power}\label{planning-for-realistic-participation}`
- **Proposed:** `\subsection{Under control-only spillover, the saturation designs needed the fewest colleges and detected the smallest effects}\label{planning-for-realistic-participation}`
- **Reason:** Table S10 (`tab:nc-power`) gives 28 colleges for both saturation designs against 29 to 32 or more than 58 for the others. It also gives the smallest control-only detectable effects at 29 colleges (1.46 and 1.47) and at 58 (1.24 and 1.25, against 1.45 or more). Figure S8 shows the colleges needed. The 28-against-29 margin over Balanced Halves is small, but the detectable-effect gap is clear. Under both-arms spillover the designs hardly differed, apart from Checkerboard and Isolation Buffer.

### NH-17 · Section 16.2 · both builds
- **Current:** `\subsection{Interval calibration}\label{interval-calibration}`
- **Proposed:** `\subsection{\texorpdfstring{Only the corrected \(t\) interval met the coverage target in the main case study}{Only the corrected t interval met the coverage target in the main case study}}\label{interval-calibration}`
- **Reason:** In the main case study, the corrected \(t\) interval gave 94.1% to 94.9% coverage, "the only interval that met the target" of 94.0% to 96.0% (Wald 92.6% to 93.6%; Tables S11 `tab:nc-calibration` and S12 `tab:interval-comparison`). The qualifier is needed because in the separate 500-fit analysis every method, the corrected \(t\) included, fell below 94.0% in some design and regime.

### NH-18 · Section 16.3 · both builds
- **Current:** `\subsection{Guidance under partial participation}\label{guidance-partial-participation}`
- **Proposed:** `\subsection{The saturation designs' advantage shrank as fewer colleges participated}\label{guidance-partial-participation}`
- **Reason:** `tab:nc-guidance` gives control-only ratios to SRS of 0.78 with 58 colleges, 0.90 to 0.91 with 40, 0.95 to 0.96 with 29 and 0.99 with 20. The text explains that non-participating neighbors supply exposure variation.

### NH-19 · Section 16.4 · both builds
- **Current:** `\subsection{Sensitivity analyses}\label{nc-sensitivity-results}`
- **Proposed:** `\subsection{The design comparison held under most alternative models and spillover forms}\label{nc-sensitivity-results}`
- **Reason:** This is the section's opening sentence ("The design comparison held ... in most spatial-model comparisons"). `tab:nc-car` keeps the saturation designs at 0.76 to 0.82 times SRS under control-only spillover, and `tab:nc-forms` keeps the count form's ordering. "Most" covers the exceptions the section names: CAR fits to spatial lag data, the Bayesian CAR fit, and binary exposure's assignment-dependent identifiability.

---

## Headings kept as they are

| Section | Build | Current LaTeX line | Why it stays |
|---|---|---|---|
| 6 | both | `\section{Grid benchmark under queen contiguity}\label{detailed-grid-results-and-allocation-risk-pilot}` | Umbrella over four subsections with different findings. A single headline would repeat NH-3 or NH-1, and DeSurv keeps such parent headings descriptive. |
| 12 | both | `\section{Robustness, the SRS benchmark and allocation-specific risk}\label{robustness-srs-allocation-risk}` | Umbrella over 12.1 to 12.3. Its natural headline ("No single design was best under both spillover regimes") is used for 12.1 (NH-10), and the current title keeps "SRS benchmark" visible in the TOC. |
| 16 | both | `\section{Participation, power and sensitivity in the case study}\label{app:nc-ext}` | Umbrella over four subsections (power, calibration, partial participation, model sensitivity) with no shared result. |
| 16.5 | unabridged only | `\subsection{Implications for trial planning}\label{nc-planning-implications}` | A short planning discussion, not a result. A headline would turn it into a recommendation, which the paragraph itself declines ("do not establish ... a final trial protocol"). |

**Not in scope, unchanged:** the run-in `\sihead{...}` lead-ins in Part II ("Colleges needed and
smallest detectable effect", "Spatial model", "Bayesian CAR model", "Bayesian spatial lag
model", "Spillover form", and the unabridged-only "Formal tests on the case study results").
They are bold run-in labels that do not appear in the TOC.

**Part I check:** no Methods heading is now misleading. Two labels carry stale names
(`detailed-grid-results-and-allocation-risk-pilot` on Section 6,
`robustness-across-blocks-and-the-recommended-design` on 12.1). Labels never print, so both
stay unchanged to protect existing `\ref`s.

**Ripple check after applying:** Chapter 3 (`Thesis_Chapter.qmd`) uses its own appendix
headings, which are not shared with the SI, so these edits need not carry over unless you want
the chapter appendix to follow the same convention. `paper/ctj_manuscript/notes/si-reorg-map.md` mentions some old
section titles and could be updated for reference only.

---

## Decisions (author, 2026-11-26)

- NH-1 rejected; the current heading stays.
- Accepted exactly as proposed: NH-2, NH-4, NH-6, NH-8, NH-9, NH-10, NH-17, NH-18, NH-19, and NH-13 option A.
- Accepted with revised wording (replacing the proposal above):
  - NH-3: Only the spillover regime changed the design ranking
  - NH-5: Formal tests separated most pairs of designs
  - NH-7: Estimation error depends on \(\beta\) (with the `\texorpdfstring` plain-text form "Estimation error depends on beta")
  - NH-11: The saturation designs had lower MSE than SRS only under control-only spillover
  - NH-12: The saturation designs' control-only advantage held for their worst allocations
  - NH-14: Both saturation designs had lower control-only MSE than SRS in every case-study year
  - NH-15: In the case study, the saturation designs also reduced worst-decile control-only MSE
  - NH-16: Under control-only spillover, the saturation designs needed the fewest colleges
- Every `\label` is unchanged; only visible heading text changed.
- Applied 2026-11-26 in `Supplementary_Information.tex` (both builds).
