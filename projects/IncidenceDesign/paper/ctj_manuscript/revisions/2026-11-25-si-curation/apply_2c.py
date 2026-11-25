#!/usr/bin/env python3
"""Apply approved Phase 2c SI edits by ID (sign-off: 2c-signoff.md).

Author: Andrew Walther. Created: 2026-11-25.
Usage: python3 apply_2c.py <ID> [<ID> ...]   (run from anywhere)
Each edit is an exact (old, new[, count[, 'regex']]) replacement that must match exactly
count times (default once).
"""
import re
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
    # Batch 2: "identified" -> "identifiable" in the identifiability sense only (FCL-028).
    'I1': [('treatment and spillover cannot be separately identified;',
            'treatment and spillover are not separately identifiable;')],
    'I2': [('treatment is identified but control-only spillover',
            'treatment is identifiable but control-only spillover')],
    'I3': [('\\(\\beta\\) is not identified: the fitting routine',
            '\\(\\beta\\) is not identifiable: the fitting routine')],
    'I4': [('\\(\\beta\\) is\nidentified, but imprecisely.',
            '\\(\\beta\\) is\nidentifiable but imprecisely estimated.')],
    'I5': [('Under rook contiguity Checkerboard is not identified ($\\mathbf{W}',
            'Under rook contiguity, $\\beta$ is not identifiable for Checkerboard ($\\mathbf{W}')],
    # I6 (seven captions) and I8 (Table S14 caption) share one phrase: 8 matches in all.
    'I6+I8': [('$\\beta$ is not identified for Checkerboard',
               '$\\beta$ is not identifiable for Checkerboard', 8)],
    'I7': [('($\\beta$ not identified)', '($\\beta$ not identifiable)', 3)],
    # I9: five prose places, some broken across source lines (apply after I3).
    'I9': [(r'(\\\(\\beta\\\) is not\s+)identified\b', r'\1identifiable', 5, 'regex')],
    'I10': [('treatment effect is identified on this irregular graph',
             'treatment effect is identifiable on this irregular graph'),
            ('exact non-identification under rook', 'exact non-identifiability under rook')],
    'I11': [('Isolation Buffer is never identified,', 'Isolation Buffer is never identifiable,'),
            ('no strategy was identified in every fit', 'no strategy was identifiable in every fit'),
            ('The saturation designs were identified most often',
             'The saturation designs were identifiable most often'),
            ('among identified fits', 'among identifiable fits')],
    'I12': [('n.i.: not identified in every fit', 'n.i.: not identifiable in every fit')],
    # Batch 3, F6: the band note moves from the figure's axis label into the caption.
    'F6': [('The dashed red line marks the nominal 0.95 level, and the gray dashed line marks the mean coverage of simple random sampling.}',
            'The dashed red line marks the nominal 0.95 level, and the gray dashed line marks the mean coverage of simple random sampling. Shaded bands show $\\pm 1$ times the average scenario-level Monte Carlo standard error of coverage.}')],
}

def main(ids):
    text = SI.read_text()
    for i in ids:
        for edit in EDITS[i]:
            old, new = edit[0], edit[1]
            want = edit[2] if len(edit) > 2 else 1
            if len(edit) > 3 and edit[3] == 'regex':
                n = len(re.findall(old, text))
                if n != want:
                    raise SystemExit(f'{i}: expected {want} matches, found {n}: {old!r}')
                text = re.sub(old, new, text)
                continue
            n = text.count(old)
            if n != want:
                raise SystemExit(f'{i}: expected {want} matches, found {n}: {old[:70]!r}')
            text = text.replace(old, new)
    SI.write_text(text)
    print('applied', ' '.join(ids))

if __name__ == '__main__':
    main(sys.argv[1:])
