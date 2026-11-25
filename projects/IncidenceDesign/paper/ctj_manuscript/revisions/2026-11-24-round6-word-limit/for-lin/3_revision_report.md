---
title: "CTJ Manuscript - Revision Report"
author: "Andrew Walther"
date: "2026-11-25"
geometry: margin=2.2cm
fontsize: 10pt
---

This package holds `1_CTJ_revised.pdf`, the new version in the reading layout, and `2_CTJ_changes_since_annotated_draft.pdf`, which marks every change since the draft you annotated with deletions struck in red and additions underlined in blue. This report shows how we answered each comment and where to check it. The main text had 3,498 words in the annotated draft and 3,853 after your comments were applied. Cuts bring the revised manuscript to 3,435 words, under the 3,500 limit, which excludes the abstract, references, tables and figures. Our count also excludes mathematical expressions. The abstract has 422 of 425 words.

Following our meeting, the treatment effect is now $\beta$, the incidence coefficient $\tau$ and the spillover magnitude $\psi$. Where your notes used the old symbols, the revised manuscript uses the new ones. Page numbers refer to `1_CTJ_revised.pdf`.

## Your comments

| ID | You asked | Revised text | Where |
|------|--------------|----------------------|----------|
| FCL-001 | Citation at sentence end | Moved | Introduction, p. 1 |
| FCL-002 | Who is "its"? | Sentence rewritten without it | Introduction, p. 1 |
| FCL-003 | "checkerboard design" | "in a checkerboard design with neighbors defined by shared edges" | Introduction, p. 1 |
| FCL-004 | Cite the first paper | Earlier paper cited | Introduction, p. 1 |
| FCL-005 | Meaning of the sentence | Now says we judged each design by average accuracy and by the risk that the single allocation drawn gives a poor estimate | Introduction, p. 1 |
| FCL-006 | Cover SRS and Checkerboard in the first sentence | Done, the SRS sentence is gone | Introduction, p. 1 |
| FCL-007 | Sentence belongs last | Now leads into the closing statement of what we compare | Introduction, p. 1 |
| FCL-008 | Cite a paper | Baird et al., Leung and our earlier paper | Introduction, p. 1 |
| FCL-009 | Your wording | Used as written | Introduction, p. 1 |
| FCL-010 | Checkerboard sentence unclear | Two sentences now say that "every neighbor of a treated cluster is a control" and that treatment then "cannot be separated from the spillover effect" | Introduction, p. 1 |
| FCL-011 | Is this checkerboard? | Yes, "We refer to this design here as Checkerboard." | Introduction, p. 1 |
| FCL-012 | "In this article, we aim to plan" | "In this paper, we also use a case study, a simulation based on observed incidence, to recommend a design" | Introduction, p. 1 |
| FCL-013 | "further challenges"? | "Spatial interference violates the assumption" | Introduction, p. 1 |
| FCL-014 | Allocation risk unclear | The abstract defines it as "how large the error can be for the particular assignment a trial happens to draw", and Methods defines it in words before the formula | Abstract, p. 1 and Methods, p. 2 |
| FCL-015 | "pilot trial" | Used in abstract and Introduction | p. 1 |
| FCL-016 | One-line formula, noise in text | Done | Methods, Outcome model, p. 2 |
| FCL-017 | Queen and rook sentence from nowhere | Moved to the simulation design | Methods, Simulation study, p. 2 |
| FCL-018 to 020 | Remove the parenthesis on noise units | "The noise $\varepsilon$ is independent standard normal. Therefore, the intervention effect $\beta$ is in units of its standard deviation." | Methods, Outcome model, p. 2 |
| FCL-021 | "baseline outcome" confusing | Removed. The coefficient $\tau$ now "sets how strongly incidence shifts the outcome" | Methods, Outcome model, p. 2 |
| FCL-022 | "choose" but "not estimated"? | "The investigator specifies $\tau$ ... It is not estimated at the design stage" | Methods, Simulation study, p. 2 |
| FCL-023 | Sentence too broken | Split into two | Methods, Outcome model, p. 2 |
| FCL-024 | Reference for the interval | MacKinnon and White for the factor, Leyrat et al. for the t quantile, and now Lee (2004) for the Wald interval | Methods, Outcome model, p. 2 |
| FCL-025 | Spell out MSE | "mean squared error (MSE)" at first use | Introduction, p. 1 |
| FCL-026, 027 | "confidence intervals", "coverage was measured" | Both used | Methods, Outcome model, p. 2 |
| FCL-028 | "identifiable"? | Used throughout the main text. The supplement still says "identified" and changes in its next revision | Methods, p. 2 and Results, pp. 4–5 |
| FCL-029 | Is one assignment possible? | Sentence removed. Methods says Checkerboard is deterministic and High Incidence Focus has one assignment unless ties occur | Methods, Designs, p. 2 |
| FCL-030 | Why not a coefficient on S(Z)? | Equation 1 has $\psi$ S(Z) | Methods, p. 2 |
| FCL-031 | Comma | Added | Methods, Simulation study, p. 2 |
| FCL-032 to 034 | Model for X, normal and Poisson, how mapped to (0,1) | Spatial incidence is written as an equation, Poisson counts use log-rates with the same spatial filter and are rank-scaled | Methods, Simulation study, pp. 2–3 |
| FCL-035 | "surface k assignment"? | "Within surface k, the assignments, the simulated outcomes and the fitted model all used the same incidence" | Methods, Simulation study, p. 3 |
| FCL-036 | $\rho_X$? | Defined in the model and in the Table 1 caption | p. 2 and Table 1, p. 3 |
| FCL-037 | Row is not a design | "Example incidence surface (not a design)" | Table 1, p. 3 |
| FCL-038 | "region" for "cluster"? | Declined. The cluster is the unit of randomization, and "region" names the four saturation groups | Table 1, p. 3 |
| FCL-039 | North Carolina Checkerboard rule | Table 1 states the rule in full | Table 1, p. 3 |
| FCL-040 | k-means on location | k-means forms 14 groups from location, incidence rank and log population, and half of each group is treated | Methods, Case study, p. 4 and Table 1 |
| FCL-041 | Delete "so" | MCSE sentence rewritten | Methods, Simulation study, p. 4 |
| FCL-042 | "configuration" confusing | Defined as the model that generates incidence, and a surface as one draw from it | Methods, Simulation study, p. 2 |
| FCL-043 | Is this a result? | It describes a study design and stays in Methods, reworded as "We also ran a separate allocation-risk study on the grid." Results follow on p. 5 | Methods, p. 4 |
| FCL-044 | Isolation Buffer fair? | Reply below | |
| FCL-045 | What is "simulated"? | "Observed incidence and geography were fixed inputs, and outcomes were simulated" | Methods, Case study, p. 4 |
| FCL-046 | "(1)" | "Equation (1)" | Methods, Case study, p. 4 |
| FCL-047 | Parameter values unclear | Now $\tau=0$ and $\beta=1$, with "The education outcome is not expected to depend on SUD incidence." | Methods, Case study, p. 4 |
| FCL-048 | Place after simulation results? | Renamed "Case study" and kept in Methods because it describes the setup. Its results follow the grid results | Methods, p. 4 |
| FCL-049 | Mention later? | Moved to the end of the paragraph | Results, Simulation study, p. 4 |
| FCL-050 | 0.0907 not in Table 2 | Pooled sentence removed | Results, Simulation study, p. 4 |
| FCL-051 | Meaning? | Deleted | Results, p. 4 |
| FCL-052 | "(SI)"? | Defined at first use, and every pointer names a section, table or figure | Methods, p. 2 |
| FCL-053 | Why 1? | $\tau=1$ shifts the outcome across the incidence range by as much as the primary effect $\beta=1$, the harder case for incidence-based designs | Methods, Simulation study, p. 2 |
| FCL-054 | Regions | k-means on coordinates run 250 times, keeping the partition that best balanced compactness, population, clusters and counties | Methods, Case study, p. 4 |
| FCL-055 | "spatial dependence coefficient" | Used | Table 2, p. 5 and Figure 2, p. 6 |
| FCL-056 | True $\gamma$? | "includes the true spillover covariate S(Z), with $\psi$ estimated" | Table 2, p. 5 |
| FCL-057 | "Reference" | Done | Table 2, p. 5 |
| FCL-058 | Both arms, no better design? | Reply below | |
| FCL-059 | How sensitive is the winner? | New paragraph, "Incidence-Guided Saturation was not best in every setting." | Discussion, p. 8 |
| FCL-060 | "the best" or "preferred"? | Heading states the comparison, "had lower error than SRS under control-only spillover" | Results, p. 5 |

**FCL-044.** Yes, under control-only spillover. Isolation Buffer does not use incidence, yet its grid MSE was 0.74 times that of SRS, and at $\beta=3$ its MSE was lower than that of Incidence-Guided Saturation. Under both-arms spillover it had 3.80 times the MSE of SRS, and it treats fewer clusters. The Discussion says its protection is useful only when both-arms spillover is implausible and treating fewer clusters is acceptable (p. 7).

**FCL-058.** Yes. Results says no design was detectably more accurate than SRS under both-arms spillover (p. 4), and the Conclusions say "SRS remains a reasonable choice when spillover reaches both arms" (p. 8).

**FCL-059.** Under both-arms spillover the design's grid MSE was 8% above that of SRS, and under control-only spillover its advantage narrowed as $\beta$ grew until Isolation Buffer was lower at $\beta=3$.

## Points from our meeting

- **Notation (D-1).** Equation 1 uses $\beta$, $\tau$ and $\psi$ as above (Methods, p. 2).
- **Incidence coefficient (D-2).** Set by the investigator, not estimated, and fixed at $\tau=1$ on the grid (Methods, p. 2).
- **Case study (D-3).** "Application" became "case study" in the headings, abstract and Figure 3 (pp. 1, 4, 5, 7).
- **Interval citations (D-4).** All three intervals are cited (Methods, p. 2).
- **Estimating $\tau$ after a trial (D-5).** This is future work in the dissertation chapter, and the paper has no sentence on it yet. A sensitivity analysis over a range of $\tau$ values is being added to the supplement.
- **Other meeting items.** Methods now says the North Carolina data enter only the case study and explains the common (0,1) incidence scale (p. 3).

## Further revisions

These came from our own further review after your comments were applied.

- **Content.** The fuller definition of $\tau$ at first use now says what $\tau=0$ and $\tau>0$ mean (p. 2). Methods notes that our earlier study used binary spillover, unlike the proportional spillover here, and points to the supplement comparison of spillover forms (p. 2). Results now states that incidence configuration did not change which design had the lowest control-only MSE (p. 4). The Wald interval cites Lee (2004).
- **Wording.** We removed semicolons, colons in running prose and ", so" joins, took results out of parentheses, and opened sentences with "Therefore," where the logic needed it.
- **Abstract.** The case-study sentence was shortened, and the contiguity sentence was folded into the list of sensitivity analyses.
- **Final checks.** Figure 3 now reads "Case study," two orphaned sentences were reconnected (pp. 2, 6), and a DOI in the references now prints correctly.

## Length cuts

Each cut is recorded in `paper/ctj_manuscript/notes/trimmed-material.md`.

- **Introduction.** The sentence saying we compare no numbers between the two studies was dropped. We still make no such comparison.
- **Methods, estimation.** The `lagsarlm` agreement check moved to SI Section 3, and a repeated SI pointer was dropped.
- **Methods, designs.** The paragraph on how designs handle incidence and exposure, and the SRS sentence, are covered by Table 1.
- **Methods, simulation study.** The (0,1) scale caveat and the two MCSE sentences were each shortened to one sentence.
- **Methods, case study.** The IRB protocol is in Ethical considerations, the heart-failure correction and rounding rule in SI Section 5, and the split-half check in SI Section 5.4. The Checkerboard rule is in Table 1.
- **Results.** Precision targets appear only in Methods, and worst-decile detail moved to SI Section 14.
- **Discussion and Continued work.** Power, reach and omitted-spillover sentences that repeated Results were removed, and two paragraphs on High Incidence Focus and the spillover regime were merged.

## Still open

The supplement revision comes next. The coauthor approval statement and the IRB determination remain placeholders.
