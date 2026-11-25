#!/usr/bin/env python3
r"""
Script: round6_apply.py
Purpose: Apply the author-approved round 6 CTJ edits (decisions.md) by ID: accepted cuts
         from cuts.py (with the author's edited wording where it differs), new items
         (C-16, TAU-1, A1) and, once decided, accepted human-writer edits from hw_edits.py.
Author: Andrew Walther
Created: 2026-11-25
Dependencies: Python 3 standard library; cuts.py and hw_edits.py beside this script;
              paper/tools/ctj_wordcount.py (texcount)

Every replacement's old text must occur exactly once in the text it is applied to, in
order, or the script stops without writing anything.

Usage (project root):
  python3 .../round6_apply.py --out <file.tex>   -> writes the edited manuscript elsewhere
  python3 .../round6_apply.py --in-place         -> edits CTJ_Manuscript.tex and CTJ_Exhibits.tex
"""
import argparse
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
MAIN = ROOT / 'paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex'
EXHIBITS = ROOT / 'paper/ctj_manuscript/manuscript_ctj/CTJ_Exhibits.tex'


def load(name):
    """Import a sibling script as a module."""
    spec = importlib.util.spec_from_file_location(name, HERE / f'{name}.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cuts = load('cuts')

# Accepted cut IDs, in the order they are applied (decisions.md, 1b) ----
ACCEPTED_CUTS = ['T14', 'T12', 'C10', 'C3', 'N2', 'C4', 'T4', 'T11', 'N3', 'T5', 'C1', 'T6',
                 'T6b', 'N7', 'N6', 'N12b', 'T2', 'T10', 'T2b', 'T7', 'C5', 'C8', 'N12', 'T13',
                 'T3', 'C6', 'N13', 'C7', 'R2', 'R3', 'N5']

# Author's edited wording replaces the proposed text for these IDs ----
EDITED = {
 'T11': r""" to put them on the same $(0,1)$ scale as the other configurations, so that one value of $\tau$ shifts the outcome by the same amount from the lowest to the highest incidence in every configuration.""",
 'N3': r"""Each performance metric is an average over a finite number of simulated trials, and its Monte Carlo standard error (MCSE) measures the simulation error in that average.\cite{morris_using_2019}""",
 'T10': r"""For Isolation Buffer, the two regimes define the same model because no treated cluster has a treated neighbor, and its small differences between regimes are Monte Carlo error.""",
 'T7': r"""However, each worst decile held only ten assignments. Under both-arms spillover, most of the variance in Incidence-Guided Saturation's estimated risk was Monte Carlo error (SI Section 14). Therefore, small both-arms differences do not show better downside protection.""",
 'T3': r"""The low both-arms error of High Incidence Focus in the case study (Figure~\ref{fig:nc-mse}) came from a fixed assignment and was sensitive to incidence adjustment. Therefore, it does not support a general recommendation. Which regime applies depends on how an intervention spreads between neighboring areas, and investigators need to judge this when planning. Saturation keeps both arms in every region while varying the local share treated, which helps separate the direct effect from neighboring exposure.""",
 'N5': r"""The investigator specifies $\tau$ according to whether incidence is expected to shift the outcome. It is not estimated at the design stage, and we set $\tau=1$ in the grid simulation. Incidence ranges from 0 to 1, and this value raises the outcome by one noise standard deviation from the lowest to the highest incidence before spatial spread. That shift equals the primary intervention effect $\beta=1$. Treatment that follows incidence is then partly confounded with it, the harder case for designs that use incidence.""",
}

# New items not in cuts.py, applied after the cuts (old, new) ----
NEW = [
 ('C-16', r"""with $N-p$, as the HC1 covariance estimator does,\cite{mackinnon_heteroskedasticity_1985}""",
          r"""with $N-p$, a standard small-sample correction,\cite{mackinnon_heteroskedasticity_1985}"""),
 ('TAU-1', r"""The incidence coefficient $\tau$ sets how strongly incidence shifts the outcome.""",
           r"""The incidence coefficient $\tau$ sets how strongly incidence shifts the outcome, holding treatment and spillover fixed. With $\tau=0$, incidence guides only the allocation, as in the case study. With $\tau>0$, clusters with higher incidence have higher outcomes, as on the grid."""),
 ('A1', r"""more reliably than the Wald interval, which uses""",
        r"""more reliably than the Wald interval,\cite{lee_asymptotic_2004} which uses"""),
]

# Accepted human-writer IDs (decisions.md, 1c); filled in after the 1c decisions ----
ACCEPTED_HW = ['HW-01', 'HW-02', 'HW-07', 'HW-09', 'HW-10', 'HW-12', 'HW-13',
               'HW-16', 'HW-17', 'HW-21', 'HW-22', 'HW-23', 'HW-25', 'HW-26', 'HW-27',
               'HW-30', 'HW-31', 'HW-34', 'HW-37', 'HW-38', 'HW-39',
               'HW-41', 'HW-42', 'HW-43', 'HW-44', 'HW-45', 'HW-46', 'HW-47', 'HW-49',
               'HW-50', 'HW-51']
HW_EDITED = {}

# Human-writer items accepted with the author's own wording, as full (old, new) pairs ----
HW_CUSTOM = [
 ('HW-03', r"""Incidence-Guided Saturation, which treats larger shares of clusters in higher-incidence regions, had the lowest or near-lowest estimation error and lower allocation risk than SRS, while giving the highest-burden region the largest treated share.""",
           r"""Incidence-Guided Saturation had the lowest or near-lowest estimation error and lower allocation risk than SRS. It also gave the highest-burden region the largest treated share."""),
 ('HW-05', r"""in the same arm; we call this design Checkerboard here.\cite{walther_spatialcrt_2026}""",
           r"""in the same arm.\cite{walther_spatialcrt_2026} We refer to this design here as Checkerboard."""),
 ('HW-08', r"""A trial realizes only one assignment, so we evaluated both""",
           r"""A trial draws only one assignment. Therefore, we evaluated both"""),
 ('HW-11', r"""shows an example assignment of each on the grid and summarizes its adaptation to the case study, and the supplementary information (SI) gives full specifications.""",
           r"""shows an example assignment of each on the grid, and summarizes its adaptation to the case study. The supplementary information (SI) gives full specifications."""),
 ('HW-03b', r"""On the grid with treatment effect $\beta=1$, Incidence-Guided Saturation had a control-only mean MSE""",
            r"""On the grid with treatment effect $\beta=1$, Incidence-Guided Saturation, which treats larger shares of clusters in higher-incidence regions, had a control-only mean MSE"""),
 ('HW-19', r"""These are general simulation settings; the North Carolina data only enter the case study.""",
           r"""These are general simulation settings. The North Carolina data enter only the case study."""),
 ('HW-24', r"""Treatment budgets therefore differed between designs (Table~\ref{tab:allocation-rules}), and we report them rather than equalizing them after the fact.

We report each design's number of treated clusters""",
           r"""Therefore, treatment budgets differed between designs (Table~\ref{tab:allocation-rules}). Rather than equalizing budgets after the fact, we report each design's number of treated clusters"""),
 ('HW-28', r"""Isolation Buffer and Incidence-Guided Saturation were nearly level at $\beta=2$ (0.160 and 0.160), and at $\beta=3$ Isolation Buffer was lower (0.160 versus 0.182), because for most designs an error in estimating $\rho$ carries into $\hat\beta$ in proportion to $\beta$, which Isolation Buffer largely avoids (SI Section 6).""",
           r"""At $\beta=2$, Isolation Buffer and Incidence-Guided Saturation had nearly equal MSE (0.160 and 0.160). At $\beta=3$, Isolation Buffer was lower (0.160 versus 0.182), because it largely avoids the error in estimating $\rho$ that for most designs passes into $\hat\beta$ in proportion to $\beta$ (SI Section 6)."""),
 ('HW-29', r"""underestimated $\beta$ (bias $-0.121$ to $-0.400$ and $t$-based coverage 54--91\% across designs and contiguity types at $\beta=1$), yet some designs had lower MSE than under the correct model because their reduced variance outweighed the added bias. A low MSE alone therefore does not show""",
           r"""underestimated $\beta$. At $\beta=1$, bias ranged from $-0.121$ to $-0.400$ and $t$-based coverage was 54--91\% across designs and contiguity types. Some designs still had lower MSE than under the correct model, because their reduced variance outweighed the added bias. Therefore, a low MSE alone does not show"""),
 ('HW-33', r"""The two changes were made together, and their separate contributions therefore cannot be distinguished.""",
           r"""The two changes were made together. Therefore, their separate contributions cannot be distinguished."""),
 ('HW-36', r"""investigators need to judge this when planning.""",
           r"""investigators need to judge which is more plausible when planning."""),
 ('TH-1', r"""independent standard normal, and the intervention effect $\beta$ is therefore in units of its standard deviation.""",
          r"""independent standard normal. Therefore, the intervention effect $\beta$ is in units of its standard deviation."""),
 ('TH-2', r"""A design's ranking on the grid therefore need not carry over""",
          r"""Therefore, a design's ranking on the grid need not carry over"""),
]

# Caption and Table 1 cell edits with the author's own wording (CTJ_Exhibits.tex) ----
CAP_CUSTOM = [
 ('HW-40', r"""(blue, treated; white, control)""", r"""(blue = treated, white = control)"""),
 ('HW-48', r"""Designs share surfaces and outcome noise, so Table~\ref{tab:grid-srs} gives""",
           r"""Designs share surfaces and outcome noise. Therefore, Table~\ref{tab:grid-srs} gives"""),
 ('TC-1', r"""Treat alternating clusters in a checkerboard pattern, so no two clusters""",
          r"""Treat alternating clusters in a checkerboard pattern, in which no two clusters"""),
 ('TC-2', r"""until none remain, so no two treated clusters are neighbors (about 21""",
          r"""until none remain. No two treated clusters are then neighbors (about 21"""),
]


def hw_rows(part):
    """Human-writer rows keyed by the IDs of human-writer-review.md (HW-01 to HW-51).

    hw_edits.py numbers its tuples differently after HW-29 (it has an HW-29b), so the
    review's IDs are assigned here by position, which both files share.
    """
    hw = load('hw_edits')
    keyed = {f'HW-{i + 1:02d}': r for i, r in enumerate(hw.E + hw.CAP)}
    keep = hw.E if part == 'E' else hw.CAP
    return {k: r for k, r in keyed.items() if any(r is x for x in keep)}


def replace_once(text, old, new, label):
    """Replace old by new, requiring exactly one occurrence."""
    n = text.count(old)
    if n != 1:
        sys.exit(f'{label}: old text found {n} times; nothing written')
    return text.replace(old, new)


def edit_main(tex):
    """Apply cuts (with edits), new items and accepted HW edits to the manuscript."""
    by_id = {it['id']: it for it in cuts.ITEMS}
    for k in ACCEPTED_CUTS:
        tex = replace_once(tex, by_id[k]['old'], EDITED.get(k, by_id[k]['new']), k)
    for k, old, new in NEW:
        tex = replace_once(tex, old, new, k)
    rows = hw_rows('E')
    for k in ACCEPTED_HW:
        if k in rows:
            tex = replace_once(tex, rows[k][2], HW_EDITED.get(k, rows[k][3]), k)
    for k, old, new in HW_CUSTOM:
        tex = replace_once(tex, old, new, k)
    return tex


def edit_exhibits(tex):
    """Apply accepted caption HW edits to the exhibits file."""
    rows = hw_rows('CAP')
    for k in ACCEPTED_HW:
        if k in rows:
            tex = replace_once(tex, rows[k][2], HW_EDITED.get(k, rows[k][3]), k)
    for k, old, new in CAP_CUSTOM:
        tex = replace_once(tex, old, new, k)
    return tex


def main():
    """Write the edited sources to --out or in place, and report the body count."""
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path)
    ap.add_argument('--in-place', action='store_true')
    args = ap.parse_args()
    main_tex = edit_main(MAIN.read_text())
    exhibits_tex = edit_exhibits(EXHIBITS.read_text())
    if args.in_place:
        MAIN.write_text(main_tex)
        EXHIBITS.write_text(exhibits_tex)
    elif args.out:
        args.out.write_text(main_tex)
        args.out.with_name(args.out.stem + '_exhibits.tex').write_text(exhibits_tex)
    else:
        sys.exit('Give --out <file> or --in-place')
    print(f'Body words: {cuts.count(main_tex)}')


if __name__ == '__main__':
    main()
