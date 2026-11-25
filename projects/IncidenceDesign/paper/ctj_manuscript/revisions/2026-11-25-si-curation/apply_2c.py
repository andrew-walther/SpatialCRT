#!/usr/bin/env python3
"""Apply approved Phase 2c SI edits by ID (sign-off: 2c-signoff.md).

Author: Andrew Walther. Created: 2026-11-25.
Usage: python3 apply_2c.py <ID> [<ID> ...]   (run from anywhere)
Each edit is an exact (old, new) replacement that must match exactly once.
"""
import sys
from pathlib import Path

SI = Path(__file__).resolve().parents[2] / 'supplementary_info/Supplementary_Information.tex'

EDITS = {
    'R1': [('\\subsection{Robustness across blocks and the recommended design}',
            '\\subsection{Robustness across blocks}')],
    'R2': [('These results do not change the\nrecommendation.',
            'These results do not change the\ncomparison of designs.')],
    'R3': [('Among the six designs in the detailed comparison, Incidence-Guided\n'
            'Saturation combines low pooled MSE with comparatively favorable\n'
            'control-only performance. The SRS comparison in Section\n'
            '\\ref{srs-benchmark-and-allocation-specific-risk} qualifies this. A',
            'No single design was best under both spillover regimes. Under\n'
            'control-only spillover, Incidence-Guided Saturation had the lowest MSE\n'
            'of the six designs, and it had the lowest worst-case block MSE (Table\n'
            '\\ref{tab:robustness}). Under both-arms spillover, SRS and Balanced\n'
            'Quartiles were slightly more accurate (Section\n'
            '\\ref{srs-benchmark-and-allocation-specific-risk}). A')],
    'R4': [('High Incidence Focus is year-sensitive\n'
            'under control-only spillover and cannot be recommended uniformly.',
            'Under control-only spillover, the accuracy of High Incidence Focus\n'
            'depended on the planning year.')],
    'S1': [('Both variants remain in the main-text\n'
            'nine-design benchmark and in the case study, and Section\n'
            '\\ref{app:eight-design} reports the comparison of all eight designs at the\n'
            'primary scenario',
            'Both variants remain in the case study, and Section\n'
            '\\ref{app:eight-design} reports the comparison of all eight designs and SRS at the\n'
            'primary scenario')],
    'S2': [('retained designs are those of Table \\ref{tab:performance}, given here to\n'
            'four and three decimal places.',
            'retained designs are those of Table \\ref{tab:performance}.')],
    'S3': [('\\caption{Parameter grid of the original eight-design simulation. The crossed parameters form '
            '$5 \\times 5 \\times 2 \\times 4 \\times 4 \\times 2 \\times 8 = 12{,}800$ scenarios, 9,600 of them '
            'for the six designs of the detailed comparison; the settings in the lower panel are the same in every scenario.}',
            '\\caption{Parameter grid of the grid simulation. The crossed parameters form '
            '$5 \\times 5 \\times 2 \\times 4 \\times 4 \\times 2 \\times 8 = 12{,}800$ scenarios for the eight designs, 9,600 of them '
            'for the six designs of the detailed comparison. SRS, run afterward on the same grid, adds 1,600 scenarios '
            '(14,400 in all). The settings in the lower panel are the same in every scenario.}')],
    'T1': [('These tests compared each pair of designs adjacent in the mean-MSE order, and each design with SRS',
            'These tests compared each pair adjacent in the mean-MSE order, with SRS in its place in that order, and each design with SRS')],
    'T2': [('We also compared each pair of designs adjacent in the mean-MSE\n'
            'order over the 50 independent configuration-by-surface units. Every\n'
            'adjacent pair differed significantly. The mean paired differences in\n'
            'MSE, with \\(t_{49}\\) 95\\% intervals, were \\(-0.040\\) (\\(-0.043\\),\n'
            '\\(-0.037\\)) for Incidence-Guided Saturation against Balanced\n'
            'Quartiles, \\(-0.029\\)',
            'We also compared each pair adjacent in the mean-MSE order, with SRS in\n'
            'its place, over the 50 independent configuration-by-surface units. Every\n'
            'adjacent pair differed significantly except SRS and Balanced Quartiles.\n'
            'The mean paired differences in MSE, with \\(t_{49}\\) 95\\% intervals, were\n'
            '\\(-0.038\\) (\\(-0.040\\), \\(-0.035\\)) for Incidence-Guided Saturation\n'
            'against SRS, \\(-0.002\\) (\\(-0.005\\), \\(0.001\\)) for SRS against\n'
            'Balanced Quartiles (Wilcoxon \\(p = 0.15\\); SRS lower in 29 of the 50\n'
            'units), \\(-0.029\\)'),
           ('Checkerboard, with Wilcoxon \\(p \\le 4.9\\times10^{-9}\\) in every case.\n'
            'The better design of each pair had',
            'Checkerboard, with Wilcoxon \\(p \\le 4.9\\times10^{-9}\\) in every other case.\n'
            'The better design of each of those pairs had')],
    'T3': [('Wilcoxon \\(p = 9.9\\times10^{-10}\\)).',
            'Wilcoxon \\(p = 9.9\\times10^{-10}\\)). Balanced Quartiles and SRS were not\n'
            'separated. Balanced Quartiles had the lower MSE in 32 of the 50 units\n'
            '(mean difference \\(-0.001\\); \\(-0.003\\), \\(0.000\\); Wilcoxon \\(p = 0.11\\)).'),
           ('Every\nother adjacent pair differed in all 50 units.',
            'Every\nother adjacent pair, including SRS against \\(2\\times2\\) Blocking, differed in all 50 units.')],
}

def main(ids):
    text = SI.read_text()
    for i in ids:
        for old, new in EDITS[i]:
            n = text.count(old)
            if n != 1:
                raise SystemExit(f'{i}: expected 1 match, found {n}: {old[:70]!r}')
            text = text.replace(old, new)
    SI.write_text(text)
    print('applied', ' '.join(ids))

if __name__ == '__main__':
    main(sys.argv[1:])
