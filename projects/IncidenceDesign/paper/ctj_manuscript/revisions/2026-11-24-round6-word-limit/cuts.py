#!/usr/bin/env python3
r"""
Script: cuts.py
Purpose: Proposed CTJ word-limit cuts and small revisions (round 6) as exact
         current -> proposed replacements, checked against CTJ_Manuscript.tex and
         measured with the scripted word count.
Author: Andrew Walther
Created: 2026-11-24
Dependencies: Python 3 standard library; texcount; paper/tools/ctj_wordcount.py

Every item's `old` text must occur exactly once in the manuscript; the script stops
otherwise, so an applied or stale item cannot be applied twice. Nothing is changed
unless --apply is given with explicit IDs.

Usage (project root):
  python3 paper/ctj_manuscript/revisions/2026-11-24-round6-word-limit/cuts.py
      -> checks every item, prints words saved per item and the total if all applied,
         and writes signoff.md beside this script
  python3 .../cuts.py --apply T2,T4      -> applies only those IDs to the manuscript
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'paper/tools'))
import ctj_wordcount as wc  # noqa: E402

MAIN = ROOT / 'paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex'

# id, section, source, note, old, new ----
ITEMS = [
 dict(id='T14', where='Introduction (L145)', source='trim plan T14',
      note='Delete. The no-cross-study-comparison rule (AGENTS 2026-10-02) is kept by not '
           'making such comparisons; it does not need a sentence.',
      old=r""" We do not compare numerical results between the two studies.""",
      new=r""""""),
 dict(id='T12', where='Methods, estimation (L158)', source='trim plan T12',
      note='Move to the SI (already stated there).',
      old=r""" Our own ML routine reproduced the estimates of the \texttt{lagsarlm} function in the R package \texttt{spatialreg} to within $10^{-6}$.""",
      new=r""""""),
 dict(id='N2', where='Methods, allocation risk (L162)', source='new',
      note='Shorten the reason for averaging over outcomes (reviewer: avoid implying noise is removed; see L210, L222).',
      old=r"""Outcome noise alone can make one trial's error large even when the assignment's expected MSE is modest, so we estimated each $m(a)$ by averaging squared errors over many simulated outcomes before ranking assignments.""",
      new=r"""We estimated each $m(a)$ by averaging squared errors over many simulated outcomes, to limit the influence of outcome noise on the ranking of assignments."""),
 dict(id='T4', where='Methods, designs (L170)', source='trim plan T4',
      note='Table 1 already gives each design and its category; keep the point about '
           'exposure.',
      old=r"""The designs differ in how they handle incidence and spatial exposure. Balanced Quartiles spreads treatment evenly across ranked incidence groups, and High Incidence Focus aligns treatment with incidence. $2\times2$ Blocking imposes the same local proportion of treated clusters everywhere, while saturation deliberately varies that proportion between regions. Isolation Buffer prevents treated neighbors, and Checkerboard imposes strong interspersion. These restrictions change how closely a cluster's treatment tracks its share of treated neighbors.""",
      new=r"""The designs' restrictions on incidence and geography (Table~\ref{tab:allocation-rules}) change how closely a cluster's treatment tracks its share of treated neighbors."""),
 dict(id='T1', where='Methods, simulation study (L173)', source='trim plan T1 (compromise)',
      note='NOT ADVISED; cannot be combined with T11 (shared sentence). CONFLICTS WITH FCL-032 (Dr. Lin asked for the model for X). This compromise keeps '
           'the spatial formula and names the Poisson mechanism; the full move to SI '
           'Section 1.2 would save about 60 more words.',
      old=r"""The five configurations were independent Uniform$(0,1)$ incidence; spatial incidence $\mathbf{X}=\Phi\{(\mathbf{I}-\rho_X\mathbf{W})^{-1}\mathbf{u}\}$, with $\mathbf{u}$ independent standard normal, $\Phi$ the standard normal distribution function and $\rho_X=0.2$ or $0.5$; and Poisson incidence, in which death counts were drawn from Poisson distributions with log-rates $\log(35/100{,}000)+(\mathbf{I}-\rho_X\mathbf{W})^{-1}\mathbf{u}$ and 100,000 people per cluster, at the same two values of $\rho_X$. The base rate of 35 per 100,000 lies between published United States SUD rates for working-age men and women.\cite{mirzaei_years_2019} The counts were then rank-scaled to $(r_i-0.5)/100$, where $r_i$ is cluster $i$'s rank with average ranks for ties, to put them on the same $(0,1)$ scale as the other configurations.""",
      new=r"""The five configurations were independent Uniform$(0,1)$ incidence; spatial incidence $\mathbf{X}=\Phi\{(\mathbf{I}-\rho_X\mathbf{W})^{-1}\mathbf{u}\}$, with $\mathbf{u}$ independent standard normal, $\Phi$ the standard normal distribution function and $\rho_X=0.2$ or $0.5$; and Poisson death counts with spatially correlated log-rates around 35 per 100,000, between published United States SUD rates for working-age men and women,\cite{mirzaei_years_2019} at the same two values of $\rho_X$ and rank-scaled to $(0,1)$ (SI Section 1.2)."""),
 dict(id='T11', where='Methods, simulation study (L173)', source='trim plan T11',
      note='Keeps the S-04 reason (one tau means the same shift in every configuration) in the rank-scaling sentence; drops the rest. Overlaps T1, which is not advised.',
      old=r""" to put them on the same $(0,1)$ scale as the other configurations. All configurations put incidence on the $(0,1)$ scale, and one value of $\tau$ therefore sets the same shift between the lowest and highest possible incidence in every configuration, although the spread between those limits differs somewhat by mechanism.""",
      new=r""" to put them on the same $(0,1)$ scale as the other configurations, on which one value of $\tau$ gives the same shift."""),
 dict(id='N4', where='Methods, simulation study (L173)', source='new',
      note='The next subsection says where the NC data enter. Check whether a reviewer asked '
           'for this sentence before cutting.',
      old=r""" These are general simulation settings; the North Carolina data only enter the case study.""",
      new=r""""""),
 dict(id='N3', where='Methods, simulation study (L175)', source='new',
      note='Shorter MCSE definition; still defined at first use (FCL feedback, 71d0ce3).',
      old=r"""Every performance metric is estimated from a finite number of simulated datasets and therefore carries simulation (Monte Carlo) error. Its Monte Carlo standard error (MCSE) measures this uncertainty and shrinks as more datasets are simulated.\cite{morris_using_2019}""",
      new=r"""Each performance metric carries simulation error, measured by its Monte Carlo standard error (MCSE).\cite{morris_using_2019}"""),
 dict(id='T5', where='Methods, case study data (L180)', source='trim plan T5',
      note='The IRB protocol is in Ethical considerations; the heart-failure correction is '
           'already in SI Section 5.1. Keep it here if the correction must be visible in the '
           'main text.',
      old=r"""SUD deaths were identified from deidentified death certificates, obtained under University of North Carolina at Chapel Hill Institutional Review Board protocol 14-2036 (see Ethical considerations), with the case-adjudication algorithm of Gan et al.\cite{gan_factors_2019} The corrected counts include heart-failure deaths that had previously been excluded as presumed non-sudden, giving 23,523 deaths""",
      new=r"""SUD deaths were identified from deidentified death certificates (see Ethical considerations) with the case-adjudication algorithm of Gan et al.,\cite{gan_factors_2019} giving 23,523 deaths"""),
 dict(id='T6', where='Methods, case-study adaptations (L187)', source='trim plan T6',
      note='Table 1 gives the Checkerboard rule; the rounding rule moves to the SI.',
      old=r""" Checkerboard became a deterministic rule that visits clusters from most to fewest neighbors and gives each cluster the arm held by fewer of its already-assigned neighbors.""",
      new=r""""""),
 dict(id='T6b', where='Methods, case-study adaptations (L187)', source='trim plan T6',
      note='Rounding rule to the SI: SI 5.2 says only "rounded", so add the round-half-to-even rule (1.5 to 2, 0.5 to 0) there when applying.',
      old=r"""and a random half of each group was treated, with fractional halves rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), so a one-cluster group was never treated.""",
      new=r"""and a random half of each group was treated (rounding in SI Section 5.2)."""),
 dict(id='N7', where='Methods, budgets (L189)', source='new',
      note='Drop the lead-in clause.',
      old=r"""Treating half the clusters need not reach half the population, so we report each design's number of treated clusters and its share of the working-age population, a population distinct from""",
      new=r"""We report each design's number of treated clusters and its share of the working-age population, a population distinct from"""),
 dict(id='T2', where='Methods, case-study replication (L193)', source='trim plan T2',
      note='Keep one clause; the split-half check and 400-outcome detail go to SI Section 5.',
      old=r"""For allocation risk at $\rho\in\{0,0.5\}$ and $\psi=0.8$, we simulated 400 outcomes for each of the same assignments, so that each assignment's risk was estimated precisely enough to rank. As a check against overstating worst-decile error by selecting on noise, worst-decile assignments were also chosen with one half of these outcomes and evaluated with the other half.""",
      new=r"""For allocation risk at $\rho\in\{0,0.5\}$ and $\psi=0.8$, we simulated 400 outcomes per assignment and checked the worst decile by selecting it with half of these outcomes and evaluating it with the other half."""),
 dict(id='T10', where='Results, grid (L204)', source='trim plan T10',
      note='Shorter; keeps the Monte Carlo clause so Table 2 values (0.161 vs 0.159) are not read as a real difference.',
      old=r"""For Isolation Buffer, no treated cluster has a treated neighbor, so only control clusters are ever exposed, the two regimes define the same outcome model, and its small differences between regimes come from separate Monte Carlo draws.""",
      new=r"""For Isolation Buffer the two regimes define the same model, because no treated cluster has a treated neighbor, and its small differences between regimes are Monte Carlo error."""),
 dict(id='T2b', where='Results, case study (L213)', source='trim plan T2',
      note='The targets are stated in Methods (L193).',
      old=r""" and met the precision targets, an MCSE of mean MSE at most 5\% of its value and an MCSE of coverage at most 0.01.""",
      new=r""" and met the precision targets."""),
 dict(id='T7', where='Results, case-study allocation risk (L222)', source='trim plan T7',
      note='Three sentences on worst-decile uncertainty become one.',
      old=r"""Only ten sampled assignments formed each worst decile, however, so which assignments belonged to it was uncertain in every setting. Under both-arms spillover, the two outcome halves selected largely different worst deciles for Incidence-Guided Saturation, and most of the variance in Incidence-Guided Saturation's estimated risk came from Monte Carlo error (SI Section 14), so small both-arms differences cannot support a claim of better downside protection.""",
      new=r"""Each worst decile held only ten assignments, however, and under both-arms spillover most of the variance in Incidence-Guided Saturation's estimated risk was Monte Carlo error (SI Section 14), and small both-arms differences therefore do not show better downside protection."""),
 dict(id='T13', where='Discussion (L234)', source='trim plan T13 (swapped by review)',
      note='Keep the 22% to 5% numbers in Results (L224), where L232 relies on them, and delete '
           'the Discussion repeat. Reject if R-09 (FCL-059) should stay word for word.',
      old=r""" In the case study, its advantage over SRS fell from 22\% with 58 colleges to 5\% with 29.""",
      new=r""""""),
 dict(id='T3', where='Discussion (L236, L238)', source='trim plan T3',
      note='Merge the High Incidence Focus and regime paragraphs. Drops "is consistent with its control-only advantage". The FCL-058 reply should then point to this merged paragraph.',
      old=r"""Incidence can also target treatment directly, but the low both-arms error of High Incidence Focus in the case study (Figure~\ref{fig:nc-mse}) came from a fixed assignment rather than from randomization and was sensitive to incidence adjustment, so it does not support a general recommendation.

Either regime could apply, depending on how a proposed intervention spreads between neighboring areas, so investigators need to judge which is plausible when planning. Saturation varies the local proportion of treated clusters while keeping both arms in every region, which helps separate the direct effect of treatment from neighboring exposure and is consistent with its control-only advantage.""",
      new=r"""The low both-arms error of High Incidence Focus in the case study (Figure~\ref{fig:nc-mse}) came from a fixed assignment and was sensitive to incidence adjustment, and therefore does not support a general recommendation. Which regime applies depends on how an intervention spreads between neighboring areas, and investigators need to judge this when planning. Saturation keeps both arms in every region while varying the local share treated, which helps separate the direct effect from neighboring exposure."""),
 dict(id='N13', where='Discussion, continued work (L244)', source='new (replaces trim plan T9)',
      note='Limitations already says the rankings are not comparisons at equal cost, and Results '
           'gives the reach figures.',
      old=r""" Isolation Buffer treats fewer clusters and High Incidence Focus reaches only 22--24\% of the working-age population, so these results do not identify the most efficient design for a fixed implementation cost or number of participants. Before using them""",
      new=r""" Before using these results"""),
 dict(id='N6', where='Methods, case-study settings (L191)', source='new',
      note='Simulation bookkeeping; SI Section 5.3 already says it.',
      old=r""" Designs whose assignments ignore incidence, such as SRS, have the same possible assignments every year, so their settings for different years shared one simulation.""",
      new=r""""""),
 dict(id='N9', where='Methods, simulation study (L175)', source='new',
      note='The unit of analysis is explained in SI Section 4; keep the count.',
      old=r"""MCSEs and paired comparisons between designs used the 50 independent configuration-by-surface units (five configurations by ten surfaces), rather than treating settings that share a surface as independent replications.""",
      new=r"""MCSEs and paired comparisons between designs used the 50 independent configuration-by-surface units (five configurations by ten surfaces)."""),
 dict(id='N12', where='Results, case-study sensitivity (L228)', source='new',
      note='Already in SI Section 14 (regional orderings in SI 5.2); N12b adds the pointer in Methods.',
      old=r""" Alternative incidence summaries reordered some saturation regions, and the prespecified mean-rank rule was kept.""",
      new=r""""""),
 dict(id='N5', where='Methods, simulation study (L173)', source='new; compresses M-11 (FCL-047, FCL-053, D-2)',
      note='OPTIONAL. Keeps every point Dr. Lin asked for (pre-specified, not estimated, size '
           'relative to beta, the harder case) in fewer words; drops "may need adjusting".',
      old=r"""The incidence coefficient $\tau$ is set by the investigator according to whether incidence is expected to shift the outcome, and it is not estimated at the design stage. On the grid we set $\tau=1$. Incidence lies between 0 and 1, and this value therefore raises the outcome, before spatial spread, by one noise standard deviation from the lowest to the highest incidence, the same size as the primary intervention effect $\beta=1$. With incidence shifting the outcome, treatment that follows incidence is partly confounded with it, which is the harder case for designs that use incidence. The value may need adjusting for a particular outcome.""",
      new=r"""The investigator sets $\tau$ by whether incidence is expected to shift the outcome; it is not estimated at the design stage. On the grid we set $\tau=1$, and because incidence lies between 0 and 1, this raises the outcome, before spatial spread, by one noise standard deviation from the lowest to the highest incidence, the size of the primary effect $\beta=1$. Treatment that follows incidence is then partly confounded with it, the harder case for designs that use incidence."""),
 dict(id='N8', where='Introduction (L143)', source='new; compresses I-04 (FCL-002 to FCL-010)',
      note='Same example in one sentence; Methods (L160) repeats it formally. Keeps FCL-002, -003, -010 content.',
      old=r"""For example, in a checkerboard design with neighbors defined by shared edges, every neighbor of a treated cluster is a control and every neighbor of a control is treated. A cluster's share of treated neighbors is then fixed by its own arm, and the treatment effect cannot be separated from the spillover effect.""",
      new=r"""For example, in a checkerboard design with neighbors defined by shared edges, every neighbor of a treated cluster is a control and vice versa, so a cluster's share of treated neighbors is fixed by its own arm and the treatment effect cannot be separated from the spillover effect."""),
 dict(id='R2', where='Results, grid (L197)', source='audit S-03 (partial)', kind='add',
      note='Answers why five configurations are reported. Checked against the grid t + df '
           'results: under control-only spillover Incidence-Guided Saturation had the lowest '
           'MSE in all five configurations (0.129 to 0.139); under both-arms spillover SRS and '
           'Balanced Quartiles traded first place within 0.002.',
      old=r"""Incidence configuration changed each design's MSE by less than 10\%, except for High Incidence Focus (Figure~\ref{fig:grid}).""",
      new=r"""Incidence configuration changed each design's MSE by less than 10\%, except for High Incidence Focus, and did not change which design had the lowest MSE under control-only spillover (Figure~\ref{fig:grid})."""),
 dict(id='R3', where='Methods, outcome model (L156)', source='audit Q-01 (partial; settled 2026-11-14)', kind='add',
      note='Your 2026-11-14 decision: state the Project 1 / Project 2 difference in Methods and '
           'point to the existing case-study comparison of spillover forms.',
      old=r"""The incidence coefficient $\tau$ sets how strongly incidence shifts the outcome.""",
      new=r"""Spillover here is proportional to the share of treated neighbors, whereas our earlier study used binary (on-or-off) spillover;\cite{walther_spatialcrt_2026} SI Section 15.4 compares proportional, count and binary forms in the case study. The incidence coefficient $\tau$ sets how strongly incidence shifts the outcome."""),
 dict(id='N12b', where='Methods, case-study sensitivity list (L191)', source='review of N12',
      note='Results no longer reports this analysis after N12, so point to it here.', kind='add',
      old=r"""instead of mean rank.""",
      new=r"""instead of mean rank (SI Section 14)."""),
 dict(id='C1', where='Methods, case-study outcome (L185)', source='cut-list review',
      note='L158 already says incidence is fit only in the grid study; the tau = 1 sensitivity '
           '(L191) still says incidence was added.',
      old=r""" The fit included an intercept, treatment and the true spillover term but not incidence.""",
      new=r""""""),
 dict(id='C2', where='Methods, outcome model (L156)', source='cut-list review',
      note='OPTIONAL. Restates "elementwise product"; it was added as a plain-language gloss '
           '(meeting item W-04), so keep it if clinical readers need it.',
      old=r""", which multiplies the two vectors cluster by cluster""",
      new=r""""""),
 dict(id='C3', where='Methods, outcome model (L156)', source='cut-list review',
      note='Removes a colon in prose.',
      old=r"""The intervention effect $\beta$ is a direct effect: the effect of a cluster's own treatment""",
      new=r"""It is a direct effect, the effect of a cluster's own treatment"""),
 dict(id='C4', where='Methods, designs (L165)', source='cut-list review',
      note="SRS is defined in the Introduction, Table 1 and the Table 2 caption; Table 1 gives "
           "Isolation Buffer's count.",
      old=r"""SRS randomly selects 50 of the 100 clusters for treatment, and on the grid every design treated 50 clusters except Isolation Buffer, whose number treated depended on the neighbor graph.""",
      new=r"""On the grid every design treated 50 clusters except Isolation Buffer."""),
 dict(id='C5', where='Results, case-study allocation risk (L222)', source='cut-list review',
      note="Figure 3's caption and Discussion L232 already make the point; Dr. Lin asked to delete "
           "the Methods version of this sentence (FCL-029, M-09).",
      old=r""" Checkerboard and High Incidence Focus each have a single possible assignment for these inputs, so they show no variation across assignments but can still have substantial error.""",
      new=r""""""),
 dict(id='C6', where='Limitations (L241)', source='cut-list review',
      note='The omitted-spillover result is in Results (L208).',
      old=r"""The primary model used the true form of spillover, and omitting the spillover term produced biased estimates and poor coverage.""",
      new=r"""The primary model used the true form of spillover."""),
 dict(id='C7', where='Discussion, continued work (L244)', source='cut-list review',
      note='Repeats L224; removes a ", so" join.',
      old=r"""Under control-only spillover even all 58 colleges gave less than 80\% power at $\beta=1$, so the expected effect size and outcome variability must be established first.""",
      new=r"""The expected effect size and outcome variability must also be established first, because power was low under control-only spillover."""),
 dict(id='C8', where='Results, budgets (L226)', source='cut-list review',
      note='Methods L187 already says budgets differed (Table 1).',
      old=r"""Treatment budgets differed by design and year (Table~\ref{tab:allocation-rules} and SI Table S5). High Incidence Focus treated 29 clusters""",
      new=r"""High Incidence Focus treated 29 clusters (SI Table S5)"""),
 dict(id='C9', where='Introduction (L147)', source='cut-list review',
      note='OPTIONAL. Touches your I-06 edit.',
      old=r"""In this paper, we also use a case study""",
      new=r"""We also use a case study"""),
 dict(id='C10', where='Methods, estimation (L158)', source='cut-list review',
      note='The next sentence but one points to SI 15.2 again.',
      old=r"""approximate and was chosen by calibration (SI Section 15.2).""",
      new=r"""approximate and was chosen by calibration."""),
]


RECOMMENDED = ['T14', 'N8', 'T12', 'C10', 'C3', 'N2', 'C4', 'T4', 'T11', 'N9', 'N3', 'T5', 'C1',
               'T6', 'T6b', 'N7', 'N6', 'N12b', 'T2', 'T10', 'T2b', 'T7', 'C5', 'C8', 'N12', 'T3',
               'T13', 'C6', 'N13', 'C7', 'R2', 'R3']
OPTIONAL = ['N5', 'C2', 'C9', 'N4', 'T1']   # touch text added for reviewer items; N4, T1 not advised


def check(tex):
    """Confirm every item's old text occurs exactly once; return problems."""
    return [f"{it['id']}: old text found {tex.count(it['old'])} times"
            for it in ITEMS if tex.count(it['old']) != 1]


def count(tex):
    """Words (text + headers, math excluded) of a manuscript text."""
    start = tex.index('\\section{Introduction}')
    body = tex[start:tex.index('\\section*{', start)]
    return sum(t + h for _, t, h, _ in wc.texcount_sections(body))


def main():
    """Check items, report savings, write signoff.md, or apply chosen IDs."""
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', default='')
    args = ap.parse_args()
    tex = MAIN.read_text()
    problems = check(tex)
    if args.apply:
        ids = args.apply.split(',')
        chosen = [it for it in ITEMS if it['id'] in ids]
        if len(chosen) != len(ids) or any(p.split(':')[0] in ids for p in problems):
            sys.exit(f'Unknown or unmatched IDs: {ids}; problems: {problems}')
        for it in chosen:
            tex = tex.replace(it['old'], it['new'])
        MAIN.write_text(tex)
        print(f'Applied {ids}; body words now {count(tex)}')
        return
    if problems:
        sys.exit('\n'.join(problems))
    assert sorted(RECOMMENDED + OPTIONAL) == sorted(it['id'] for it in ITEMS), 'unclassified item'
    base = count(tex)
    by_id = {it['id']: it for it in ITEMS}

    def apply_ids(ids):
        t = tex
        for k in ids:
            t = t.replace(by_id[k]['old'], by_id[k]['new'])
        return t

    saved = {it['id']: base - count(tex.replace(it['old'], it['new'])) for it in ITEMS}
    rec_total = count(apply_ids(RECOMMENDED))
    opt_total = count(apply_ids(RECOMMENDED + ['N5', 'C2', 'C9']))
    lines = ['# CTJ round 6: word-limit cuts and small revisions (sign-off list)', '',
             'Author: Andrew Walther. Generated by `cuts.py` on 2026-11-24 from '
             '`CTJ_Manuscript.tex`; nothing has been applied.', '',
             f'Body now **{base}** words (text and headers, math excluded; '
             '`paper/tools/ctj_wordcount.py`). Limit 3,500.', '',
             f'- Recommended set ({len(RECOMMENDED)} items, including the two small additions R2 '
             f'and R3): **{rec_total}** words.',
             f'- Recommended plus the optional N5, C2 and C9: **{opt_total}** words.',
             '- N4 and T1 are listed but not advised: N4 is your own M-14 sentence, and T1 '
             'would cut most of the incidence model Dr. Lin asked for (FCL-032).', '',
             'Reply by ID: accept, reject, or edit. Accepted items are applied with '
             '`cuts.py --apply <IDs>`. Moved text goes to the named SI section or to '
             '`paper/thesis_chapter/notes/trimmed-material.md`. Then both layouts are rebuilt, '
             'the checkers run, and the Chapter 3 carry-over is proposed as a second list. The '
             'remaining partial audit items are handled outside the word budget: a reference '
             'for the Wald interval (FCL-024, your choice of source), "identified" in the SI '
             '(FCL-028, SI pass), "application" in the cover letter and checklist (D-3), and the '
             '", so" joins (W-03, reading-flow pass after the cut).', '',
             '| ID | Set | Where | Source | Words saved |', '|---|---|---|---|---|']
    for k in RECOMMENDED + OPTIONAL:
        it = by_id[k]
        tag = 'recommended' if k in RECOMMENDED else 'optional'
        lines.append(f"| {k} | {tag} | {it['where']} | {it['source']} | {saved[k]} |")
    lines += ['', 'Negative numbers are additions.', '', '## Items', '']
    for k in RECOMMENDED + OPTIONAL:
        it = by_id[k]
        lines += [f"### {k} ({it['where']}; saves {saved[k]})", '', it['note'], '',
                  '**Current**', '', '> ' + it['old'].strip().replace('\n\n', '\n>\n> '), '',
                  '**Proposed**', '',
                  '> ' + (it['new'].strip().replace('\n\n', '\n>\n> ') or '*(deleted)*'), '']
    (HERE / 'signoff.md').write_text('\n'.join(lines) + '\n')
    for k in RECOMMENDED + OPTIONAL:
        print(f"{k:<5} {saved[k]:>4}  {'rec' if k in RECOMMENDED else 'opt'}  {by_id[k]['where']}")
    print(f'Body {base}; recommended -> {rec_total}; plus N5, C2, C9 -> {opt_total}')

if __name__ == '__main__':
    main()
