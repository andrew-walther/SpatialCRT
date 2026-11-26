---
title: "CTJ Manuscript - Revision Report"
author: "Andrew Walther"
date: "2026-11-26"
geometry: margin=2.2cm
fontsize: 10pt
---

This package holds `1_CTJ_revised.pdf`, the new version in the reading layout, and `2_CTJ_changes_since_annotated_draft.pdf`, which marks every change since the draft you annotated with deletions struck in red and additions underlined in blue. This report shows how we answered each comment and where to check it. The main text had 3,498 words in the annotated draft and 3,853 after your comments were applied. Cuts bring the revised manuscript to 3,433 words, under the 3,500 limit, which excludes the abstract, references, tables and figures. Our count also excludes mathematical expressions. The abstract has 421 of 425 words.

Following our meeting, the treatment effect is now $\beta$, the incidence coefficient $\tau$ and the spillover magnitude $\psi$. Where your notes used the old symbols, the revised manuscript uses the new ones.

## Your comments

| ID | You asked | Revised text | Where |
|------|--------------|----------------------|----------|
| FCL-001 | Citation at sentence end | Moved | Introduction |
| FCL-002 | Who is "its"? | Sentence rewritten without it | Introduction |
| FCL-003 | "checkerboard design" | "in a checkerboard design with neighbors defined by shared edges" | Introduction |
| FCL-004 | Cite the first paper | Earlier paper cited | Introduction |
| FCL-005 | Meaning of the sentence | Now says we judged each design by average accuracy and by the risk that the single allocation drawn gives a poor estimate | Introduction |
| FCL-006 | Cover SRS and Checkerboard in the first sentence | Done, the SRS sentence is gone | Introduction |
| FCL-007 | Sentence belongs last | Now leads into the closing statement of what we compare | Introduction |
| FCL-008 | Cite a paper | Baird et al., Leung and our earlier paper | Introduction |
| FCL-009 | Your wording | Used as written | Introduction |
| FCL-010 | Checkerboard sentence unclear | Two sentences now say that "every neighbor of a treated cluster is a control" and that treatment then "cannot be separated from the spillover effect" | Introduction |
| FCL-011 | Is this checkerboard? | Yes, "We refer to this design here as Checkerboard." | Introduction |
| FCL-012 | "In this article, we aim to plan" | "In this paper, we also use a case study, a simulation based on observed incidence, to recommend a design" | Introduction |
| FCL-013 | "further challenges"? | "Spatial interference violates the assumption" | Introduction |
| FCL-014 | Allocation risk unclear | The abstract defines it as "how large the error can be for the particular assignment a trial happens to draw", and Methods defines it in words before the formula | Abstract and Methods |
| FCL-015 | "pilot trial" | Used in abstract and Introduction | Abstract and Introduction |
| FCL-016 | One-line formula, noise in text | Done | Methods, Outcome model |
| FCL-017 | Queen and rook sentence from nowhere | Moved to the simulation design | Methods, Simulation study |
| FCL-018 to 020 | Remove the parenthesis on noise units | "The noise $\varepsilon$ is independent standard normal. Therefore, the intervention effect $\beta$ is in units of its standard deviation." | Methods, Outcome model |
| FCL-021 | "baseline outcome" confusing | Removed. The coefficient $\tau$ now "sets how strongly incidence shifts the outcome" | Methods, Outcome model |
| FCL-022 | "choose" but "not estimated"? | "The investigator specifies $\tau$ ... It is not estimated at the design stage". A new sentence adds that the same model can estimate $\tau$ after a trial | Methods, Simulation study |
| FCL-023 | Sentence too broken | Split into two | Methods, Outcome model |
| FCL-024 | Reference for the interval | MacKinnon and White for the factor, Leyrat et al. for the t quantile, and now Lee (2004) for the Wald interval | Methods, Outcome model |
| FCL-025 | Spell out MSE | "mean squared error (MSE)" at first use | Introduction |
| FCL-026, 027 | "confidence intervals", "coverage was measured" | Both used | Methods, Outcome model |
| FCL-028 | "identifiable"? | Used throughout the main text and the supplement | Methods and Results |
| FCL-029 | Is one assignment possible? | Sentence removed. Methods says Checkerboard is deterministic and High Incidence Focus has one assignment unless ties occur | Methods, Treatment assignment designs |
| FCL-030 | Why not a coefficient on S(Z)? | Equation 1 has $\psi$ S(Z) | Methods, Equation 1 |
| FCL-031 | Comma | Added | Methods, Simulation study |
| FCL-032 to 034 | Model for X, normal and Poisson, how mapped to (0,1) | Spatial incidence is written as an equation, Poisson counts use log-rates with the same spatial filter and are rank-scaled | Methods, Simulation study |
| FCL-035 | "surface k assignment"? | "Within surface k, the assignments, the simulated outcomes and the fitted model all used the same incidence" | Methods, Simulation study |
| FCL-036 | $\rho_X$? | Appears in the spatial incidence equation in Methods, and the Table 1 caption names it the incidence dependence | Methods, Simulation study and Table 1 |
| FCL-037 | Row is not a design | "Example incidence surface (not a design)" | Table 1 |
| FCL-038 | "region" for "cluster"? | Declined. The cluster is the unit of randomization, and "region" names the four saturation groups | Table 1 |
| FCL-039 | North Carolina Checkerboard rule | Table 1 states the rule in full | Table 1 |
| FCL-040 | k-means on location | k-means forms 14 groups from location, incidence rank and log population, and half of each group is treated | Methods, Case study and Table 1 |
| FCL-041 | Delete "so" | MCSE sentence rewritten | Methods, Simulation study |
| FCL-042 | "configuration" confusing | Defined as the model that generates incidence, and a surface as one draw from it | Methods, Simulation study |
| FCL-043 | Is this a result? | It describes a study design and stays in Methods, reworded as "We also ran a separate allocation-risk study on the grid." Its results are reported with Figure 3 | Methods, Simulation study |
| FCL-044 | Isolation Buffer fair? | Reply below | |
| FCL-045 | What is "simulated"? | "Observed incidence and geography were fixed inputs, and outcomes were simulated" | Methods, Case study |
| FCL-046 | "(1)" | "Equation (1)" | Methods, Case study |
| FCL-047 | Parameter values unclear | Now $\tau=0$ and $\beta=1$, with "The education outcome is not expected to depend on SUD incidence." | Methods, Case study |
| FCL-048 | Place after simulation results? | Renamed "Case study" and kept in Methods because it describes the setup. Its results follow the grid results | Methods and Results |
| FCL-049 | Mention later? | Moved to the end of the paragraph | Results, Simulation study |
| FCL-050 | 0.0907 not in Table 2 | Pooled sentence removed | Results, Simulation study |
| FCL-051 | Meaning? | Deleted | Results, Simulation study |
| FCL-052 | "(SI)"? | Defined at first use, and every pointer names a section, table or figure | Methods, Outcome model |
| FCL-053 | Why 1? | $\tau=1$ shifts the outcome across the incidence range by as much as the primary effect $\beta=1$, the harder case for incidence-based designs | Methods, Simulation study |
| FCL-054 | Regions | k-means on coordinates run 250 times, keeping the partition that best balanced compactness, population, clusters and counties | Methods, Case study |
| FCL-055 | "spatial dependence coefficient" | Used | Table 2 and Figure 2 |
| FCL-056 | True $\gamma$? | "includes the true spillover covariate S(Z), with $\psi$ estimated" | Table 2 |
| FCL-057 | "Reference" | Done | Table 2 |
| FCL-058 | Both arms, no better design? | Reply below | |
| FCL-059 | How sensitive is the winner? | New paragraph, "Incidence-Guided Saturation was not best in every setting." | Discussion |
| FCL-060 | "the best" or "preferred"? | Heading states the comparison, "had lower error than SRS under control-only spillover" | Results, Case study heading |

**FCL-044.** Yes, under control-only spillover. Isolation Buffer does not use incidence, yet its grid MSE was 0.74 times that of SRS, and at $\beta=3$ its MSE was lower than that of Incidence-Guided Saturation. Under both-arms spillover it had 3.80 times the MSE of SRS, and it treats fewer clusters. The Discussion says its protection is useful only when both-arms spillover is implausible and treating fewer clusters is acceptable.

**FCL-058.** Yes. Results says no design was detectably more accurate than SRS under both-arms spillover, and the Conclusions say "SRS remains a reasonable choice when spillover reaches both arms."

**FCL-059.** Under both-arms spillover the design's grid MSE was 8% above that of SRS, and under control-only spillover its advantage narrowed as $\beta$ grew until Isolation Buffer was lower at $\beta=3$.

## Points from our meeting

- **Notation (D-1).** Equation 1 uses $\beta$, $\tau$ and $\psi$ as above.
- **Incidence coefficient (D-2).** Set by the investigator at the design stage and fixed at $\tau=1$ on the grid (Methods, Simulation study).
- **Case study (D-3).** "Application" became "case study" in the headings, abstract and Figure 3.
- **Interval citations (D-4).** All three intervals are cited (Methods, Outcome model).
- **Estimating $\tau$ after a trial (D-5).** Methods now says that after a trial the same model estimates $\tau$, less precisely when treatment follows incidence (Methods, Simulation study). SI Section 3.2 gives the evidence: we repeated the grid benchmark at $\tau = 0$, 0.25, 0.5, 1 and 2.
- **Other meeting items.** Methods now says the North Carolina data enter only the case study and explains the common (0,1) incidence scale (Methods, Simulation study).

## Further revisions

These came from our own further review after your comments were applied.

- **Content.** The fuller definition of $\tau$ at first use now says what $\tau=0$ and $\tau>0$ mean. Methods notes that our earlier study used binary spillover, unlike the proportional spillover here, and points to the supplement comparison of spillover forms. Results now states that incidence configuration did not change which design had the lowest control-only MSE. The Wald interval cites Lee (2004). The Discussion adds that on the grid, widening the quadrants' saturation spread beyond 80/20 lowered control-only MSE but raised both-arms MSE (SI Section 13).
- **Wording.** We removed semicolons, colons in running prose and ", so" joins, took results out of parentheses, and opened sentences with "Therefore," where the logic needed it. We also defined SI and queen and rook contiguity at first use, split several long sentences, dropped an unused abbreviation, and replaced "X, not Y" contrasts with plain statements. For example, the abstract now says the study "uses simulation to evaluate designs for a planned trial," and the Conclusions say the case-study outcomes "are simulated and give no evidence about whether education reduces SUD."
- **Abstract.** The case-study sentence was shortened, and the contiguity sentence was folded into the list of sensitivity analyses.
- **Final checks.** Figure 3 now reads "Case study," two orphaned sentences were reconnected, and a DOI in the references now prints correctly.

## Length cuts

Each cut is recorded in `paper/ctj_manuscript/notes/trimmed-material.md`.

- **Introduction.** The sentence saying we compare no numbers between the two studies was dropped. We still make no such comparison.
- **Methods, estimation.** The `lagsarlm` agreement check moved to SI Section 3, and a repeated SI pointer was dropped.
- **Methods, designs.** The paragraph on how designs handle incidence and exposure, and the SRS sentence, are covered by Table 1.
- **Methods, simulation study.** The (0,1) scale caveat and the two MCSE sentences were each shortened to one sentence.
- **Methods, case study.** The IRB protocol is in Ethical considerations, the heart-failure correction and rounding rule in SI Section 5, and the split-half check in SI Section 5.4. The Checkerboard rule is in Table 1.
- **Results.** Precision targets appear only in Methods, and worst-decile detail moved to SI Section 15.
- **Discussion and Continued work.** Power, reach and omitted-spillover sentences that repeated Results were removed, and two paragraphs on High Incidence Focus and the spillover regime were merged.

## Still open

The supplement has been revised to match this version. The coauthor approval statement and the IRB determination remain placeholders.
