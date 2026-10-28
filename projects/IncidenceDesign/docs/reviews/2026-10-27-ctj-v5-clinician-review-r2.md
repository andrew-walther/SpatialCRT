# CTJ v5 clinician/epidemiologist review, round 2 (2026-10-27)

Reader: an epidemiologist or trialist who reads *Clinical Trials* but is not a spatial statistician. Source: `CTJ_Reading.pdf`, with phrases quoted from `CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`. Highlighted placeholders, the kept colons, "Spatial blocking" and the pending reason for service areas are left alone. Every rewrite is no longer than the text it replaces, or names a cut.

## Findings

### BLOCKING

**B1. Abstract Methods and Results, and the whole paper: the two spillover regimes are never explained in plain words or tied to the pilot.** Every conclusion depends on the regime, but the abstract gives only "(control clusters only or both arms)". The body explains the regimes once, in equation terms. A clinician cannot tell which regime an education pilot is likely to face, so "consider saturation when spillover primarily reaches controls" gives no guidance they can act on.
- Abstract rewrite: "spillover regime (benefit reaching only nearby control clusters, or nearby clusters in both arms)" (same length as the current parenthesis plus 6 words). Pay for it by cutting "Repeated outcomes simulated for each sampled assignment showed how MSE varies across assignments." The criteria sentence already defines allocation risk.
- Discussion, "Plausible spillover regimes should be considered before an allocation is chosen." Replace with: "Investigators should judge whether trained professionals' practices reach only neighboring control colleges or also reinforce neighboring treated colleges." (21 words against 9. Cut the variants sentence in Treatment assignment designs, which repeats the Table 1 caption and saves 25.)

### SHOULD

**S1. Methods, Simulation study, "The parameter grid crossed".** "Grid" also names the 10×10 study region. Use "The simulation crossed".

**S2. "Detectable" is used in Results and the abstract but is defined only in the Table 2 caption.** At the end of the MCSE paragraph, add: "A difference is called detectable when its paired 95% interval excludes no difference." Pay for it by cutting "without comparing numerical results between the two studies" (Introduction).

**S3. Results 1, opening sentence, "...for the five incidence configurations, which changed MSE by less than 10%..."** The referent of "which" is unclear, and I read it twice. Rewrite: "Table 2 reports the queen grid benchmark at τ = 1. Incidence configuration changed MSE by less than 10% for every design except High Incidence Focus (Figure 2)."

**S4. Methods, "so it is not an average total effect."** Clinicians will not parse this. Rewrite: "τ is the effect of a cluster's own treatment before outcomes spread between neighbors, not the total effect."

**S5. Methods, "identification rests mainly on clusters at the grid boundary."** The reason is missing. If the authors confirm it, add: "every interior cluster has exactly half its queen neighbors treated". Pay for it by cutting "(cannot be estimated)", since "not identified" is already glossed by "cannot be separated".

**S6. Application adaptation, the spatial-blocking sentence (about 70 words, ending "so a one-cluster group was never treated").** It is too long, and "half... rounded to the nearest even number" needs two readings. Split it after "log-population". Then: "A random half of each group was treated, rounding an odd group's half to the nearest even number (a group of three treated two, and a single cluster was never treated)."

**S7. The downside criterion has five names: allocation risk, worst-decile mean, "allocation downside", "unfavorable allocation" and "poor realized allocation".** Use "allocation risk" for the concept and "worst-decile MSE" for the summary. Discussion: "Average MSE and allocation risk should be evaluated together."

**S8. "Under ..." opens nine sentences.** Examples are "Under control-only spillover, Incidence-Guided Saturation's mean MSE..." and "Under both-arms spillover, SRS mean MSE...". Move the regime later in some of them: "Incidence-Guided Saturation's control-only mean MSE was 21.1–29.9% lower...", and "With both-arms spillover, every design except...".

**S9. Abstract and Conclusions leave out the power result.** No design reached 80% power under control-only spillover with all 58 colleges. A planner reading only the abstract would think the control-only recommendation is ready to use. Abstract Conclusions: replace "whether confidence intervals reach their nominal coverage" (jargon) with "whether the available colleges give adequate power" (same length). The power result is already in the body.

**S10. Figures 3 and 4 label the NC adaptation "2×2 Blocking", but the text says "spatial blocking".** Add to both captions: "In the application, 2×2 Blocking denotes spatial blocking." Pay for it by cutting "budgets" from the Figure 4 list, since Table 1 gives them.

**S11. Figure 3 caption fragments ("Simulation study, all 50 grid surfaces at...". "Application study, ρ...").** These are not sentences, and "50 grid surfaces" conflicts with "50 configuration-by-surface units". Rewrite: "Grid results use all 50 configuration-by-surface units at ρ ∈ {0, 0.5} and γ ∈ {0.5, 0.8}, and application results use ρ ∈ {0, 0.5} and γ = 0.8 in each year." Also drop the opening "Allocation risk:" colon by starting "Distribution of assignment-specific MSE...".

**S12. Figure 2 caption does not state its message.** As written, Table 2 and Figure 2 look redundant. Open with: "Incidence configuration changes mean MSE little for every design except High Incidence Focus." Pay for it by cutting "Parameter-specific displays are in the SI."

**S13. Methods, "Settings with the same possible assignments, such as SRS in different years, shared one simulation."** The reason is missing, and I read it twice. Rewrite: "SRS ignores incidence, so its settings in different years shared one simulation."

### NIT

- **N1. Capitalization.** "Incidence-guided saturation used four fixed..." (Methods) and the Figure 1 caption, "used by incidence-guided saturation". Capitalize both to match the design name.
- **N2. Treatment assignment designs, "every design treated 50 clusters except Isolation Buffer".** Add "on the grid". The NC budgets differ.
- **N3. Symbols undefined in exhibits.** The Table 1 caption uses "(ρ_X = 0.5)" before the text defines it; write "(incidence dependence ρ_X = 0.5)". Table 2 says "four ρ values and four γ values"; write "four outcome-dependence (ρ) and four spillover-magnitude (γ) values".
- **N4. Table 2 caption order.** "Simple random sampling randomly selects..." comes one sentence before "(SRS)" is defined. Swap the two sentences.
- **N5. Introduction, "which address different objectives and constraints".** Vague filler. Delete it.
- **N6. Treatment assignment designs, "These restrictions change how treatment relates to the spillover covariate."** "Spillover covariate" is undefined. Use "...change how closely a cluster's treatment tracks its share of treated neighbors."
- **N7. Application, "It followed Equation 1 with β = 0, τ = 1".** Gloss τ here as "(one residual standard deviation)" instead of first glossing it in the power paragraph.
- **N8. Limitations, "94–96% in the application at every number of participating colleges".** This is a reassurance placed in Limitations, and it sits beside 94.2–94.9% in Results. Say "in the reduced-college analyses".
- **N9. Abstract length (411 words).** Check it against the journal's limit. I am not sure what the limit is.

## Author's v5 concerns

| Concern | Status |
|---|---|
| Abstract: contamination framing, two criteria, allocation not outcome, education outcome, Methods flow, τ = 1, precision sentence | RESOLVED |
| Abstract Conclusions: budget/population share, no "interval calibration" | PARTLY ("nominal coverage" still jargon; S9) |
| Intro: colons, "study region", alternating allocation, identifiability | RESOLVED |
| Intro: stuffy, long application paragraph | PARTLY (N5) |
| Methods: order, queen/rook citation and roles, β = 1, 35/100k, ten surfaces, MC error, risk-equation terms, replication rule | RESOLVED |
| Methods: choppy or ambiguous sentences | PARTLY (S6, S13) |
| "Pilot" language; no "Because" openers | RESOLVED |
| Application: rules explained for a nontechnical investigator | PARTLY (S6) |
| Results: define terms, name precision targets | PARTLY ("detectable", S2) |
| Results: MSE ranges, run totals, regime/transfer sentences, "grid pilot", key results only | RESOLVED |
| "Budgets differed (...)" fragment | RESOLVED in text; Figure 3 caption fragments (S11) |
| Discussion: plain Limitations/Conclusions | PARTLY (S7, S9) |
| Continued work: endpoint, staffing | RESOLVED |
| Figure 2 vs 3 distinct; Figure 3 caption length | PARTLY (S11, S12) |
| Figure 2 intervals; Table 2 SRS 50 of 100 | RESOLVED |
| Overall flow | PARTLY (S3, S8, B1) |
