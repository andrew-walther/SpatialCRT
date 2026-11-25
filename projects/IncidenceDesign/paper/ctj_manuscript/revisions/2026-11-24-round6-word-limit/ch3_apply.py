#!/usr/bin/env python3
r"""
Script: ch3_apply.py
Purpose: Apply the author-approved Phase 1e Chapter 3 carry-over (ch3-signoff.md) to a
         Chapter 3 Quarto source, by ID; each old text must occur exactly once.
Author: Andrew Walther
Created: 2026-11-25
Dependencies: Python 3 standard library

Usage (project root):
  python3 .../ch3_apply.py paper/thesis_chapter/Thesis_Chapter.qmd
  python3 .../ch3_apply.py <bios-dissertation>/thesis/chapter3-project2-incidence/project2-incidence.qmd
"""
import sys
from pathlib import Path

EDITS = [
 ('CH-1', r"The 95\% Wald interval uses $q = z_{0.975} = 1.96$", r"The 95\% Wald interval [@lee_asymptotic_2004] uses $q = z_{0.975} = 1.96$"),
 ('CH-2', r"The incidence coefficient $\tau$ sets how strongly cluster-level incidence shifts the outcome. Its values",
          r"The incidence coefficient $\tau$ sets how strongly cluster-level incidence shifts the outcome, holding treatment and spillover fixed. With $\tau = 0$, incidence guides only the allocation, as in the case study. With $\tau > 0$, clusters with higher incidence have higher outcomes, as on the grid. Its values"),
 ('CH-3', r"The incidence coefficient $\tau$ is set by the investigator according to whether incidence is expected to shift the outcome, and it is not estimated at the design stage. On the grid we set $\tau = 1$. Incidence lies between 0 and 1, and this value therefore raises the outcome, before spatial spread, by one noise standard deviation from the lowest to the highest incidence, the same size as the primary intervention effect $\beta = 1$.",
          r"The investigator specifies $\tau$ according to whether incidence is expected to shift the outcome. It is not estimated at the design stage, and we set $\tau = 1$ in the grid simulation. Incidence ranges from 0 to 1, and this value raises the outcome by one noise standard deviation from the lowest to the highest incidence before spatial spread. That shift equals the primary intervention effect $\beta = 1$."),
 ('CH-4', r"Half of each group, rounded, is treated at random. Treatment budgets, the numbers of treated clusters, therefore differ across designs.",
          r"Half of each group, rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), is treated at random. Therefore, treatment budgets, the numbers of treated clusters, differ across designs."),
 ('TS-01', r"The literature therefore offers", r"Therefore, the literature offers"),
 ('TS-02', r"The intervention effect $\beta$ is therefore in units of the noise standard deviation.", r"Therefore, the intervention effect $\beta$ is in units of the noise standard deviation."),
 ('TS-03', r"regardless of its own assignment, and neighboring treatment therefore adds", r"regardless of its own assignment. Therefore, neighboring treatment adds"),
 ('TS-04', r"on the $(0, 1)$ scale, and one value of $\tau$ therefore sets", r"on the $(0, 1)$ scale. Therefore, one value of $\tau$ sets"),
 ('TS-05', r"Its results are therefore a best case for each design.", r"Therefore, its results are a best case for each design."),
 ('TS-06', r"A trial realizes only one allocation, and we therefore judge each design on two criteria: the error of the estimate averaged over the allocations the design can produce, and how that error",
           r"A trial draws only one allocation. Therefore, we judge each design on two criteria, the error of the estimate averaged over the allocations the design can produce and how that error"),
 ('TS-08', r"simulated datasets and therefore carries Monte Carlo error.", r"simulated datasets. Therefore, it carries Monte Carlo error."),
 ('TS-09', r"The two regimes in Equation \ref{eq:spillover} therefore define", r"Therefore, the two regimes in Equation \ref{eq:spillover} define"),
 ('TS-10', r"spillover covariate and therefore has no omitted-variable bias.", r"spillover covariate. Therefore, it has no omitted-variable bias."),
 ('TS-11', r"share their assignments and noise, and these comparisons are therefore not independent", r"share their assignments and noise. Therefore, these comparisons are not independent"),
 ('TS-12', r"did not use split-half selection, and its worst-decile means are therefore biased upward", r"did not use split-half selection. Therefore, its worst-decile means are biased upward"),
 ('TS-13', r"We therefore use the NC case study in @sec-application as a setting-specific assessment instead of extrapolating the size of the tail protection from these settings.",
           r"Instead of extrapolating the size of the tail protection from these settings, we use the NC case study in @sec-application as a setting-specific assessment."),
 ('TS-14', r"Both average estimation error and the risk of a poor allocation therefore matter.", r"Therefore, both average estimation error and the risk of a poor allocation matter."),
 ('TS-16', r"We therefore repeated the case study for subsets", r"Therefore, we repeated the case study for subsets"),
 ('TS-17', r"highest-incidence clusters, and adjusting for incidence therefore leaves", r"highest-incidence clusters. Therefore, adjusting for incidence leaves"),
 ('TS-18', r"possible assignments and therefore the same results.", r"possible assignments and the same results."),
 ('TS-19', r"includes the spatial multiplier, and the CAR coefficient on treatment therefore estimates", r"includes the spatial multiplier. Therefore, the CAR coefficient on treatment estimates"),
 ('TS-20', r"The colleges-needed numbers above therefore come", r"Therefore, the colleges-needed numbers above come"),
 ('TS-21', r"was identifiable, and the designs therefore differed in precision", r"was identifiable. Therefore, the designs differed in precision"),
 ('TS-23', r"SRS and the saturation designs are therefore both reasonable choices.", r"Therefore, SRS and the saturation designs are both reasonable choices."),
 ('TS-24', r"Their accuracy advantage over SRS under control-only spillover therefore appears", r"Therefore, their accuracy advantage over SRS under control-only spillover appears"),
 ('TS-25', r"exposure term, and its results are therefore a best case for each design.", r"exposure term. Therefore, its results are a best case for each design."),
 ('TS-26', r"Their observed rates therefore differ in precision", r"Therefore, their observed rates differ in precision"),
 ('TS-27', r"Uncertainty in their measurement and changes between planning and implementation are therefore not propagated.", r"Therefore, uncertainty in their measurement and changes between planning and implementation are not propagated."),
 ('TS-28', r"Results from different code or settings are therefore never combined.", r"Therefore, results from different code or settings are never combined."),
 ('TS-29', r"The college mean therefore has standard deviation", r"Therefore, the college mean has standard deviation"),
 ('TS-30', r"standard deviations is therefore an individual-level standardized effect", r"standard deviations then corresponds to an individual-level standardized effect"),
 ('TS-31', r"The variance due to outcome noise therefore involves", r"Therefore, the variance due to outcome noise involves"),
 ('TS-32', r"Retaining Incidence-Guided Saturation as the representative of the saturation mechanism therefore does not rest", r"Therefore, retaining Incidence-Guided Saturation as the representative of the saturation mechanism does not rest"),
 ('TS-33', r"We therefore also ran paired Wilcoxon", r"Therefore, we also ran paired Wilcoxon"),
 ('TS-34', r"Their rows therefore repeat and are not independent", r"Therefore, their rows repeat and are not independent"),
 ('TS-35', r"They therefore describe the consistency of the ordering", r"Therefore, they describe the consistency of the ordering"),
]
# TS-22 with the author's softer wording of the line 599 recommendation (2026-11-25) ----
PENDING = [
 ('TS-22', r"""(Table \ref{tab:downside}), and they are therefore the designs we recommend over SRS. This recommendation assumes""",
           r"""(Table \ref{tab:downside}). Therefore, they are the designs to consider first in that setting. This guidance assumes"""),
]


def main():
    """Apply every edit to the given file, or stop without writing if any old text is not unique."""
    path = Path(sys.argv[1])
    text = path.read_text()
    for k, old, new in EDITS + PENDING:
        n = text.count(old)
        if n != 1:
            sys.exit(f'{k}: old text found {n} times in {path}; nothing written')
        text = text.replace(old, new)
    path.write_text(text)
    print(f'Applied {len(EDITS + PENDING)} edits to {path}')


if __name__ == '__main__':
    main()
