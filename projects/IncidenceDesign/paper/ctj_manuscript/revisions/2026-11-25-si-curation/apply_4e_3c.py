#!/usr/bin/env python3
"""Apply the approved Phase 4e and 3c edits (sign-off: 4e-3c-signoff.md).

Author: Andrew Walther. Created: 2026-11-25.
Usage: python3 apply_4e_3c.py   (run from anywhere; applies N1, N3, N4, N5, C1, C2)

The three new LaTeX blocks (N1, N4, N5) are read from the ```latex fences of the sign-off
file in that order, so the applied text is the signed-off text. Each anchor must match
exactly once, or the script stops before writing anything.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SI = PAPER / 'supplementary_info/Supplementary_Information.tex'
CTJ = PAPER / 'manuscript_ctj/CTJ_Manuscript.tex'

blocks = re.findall(r'```latex\n(.*?)```', (HERE / '4e-3c-signoff.md').read_text(), re.S)
if len(blocks) != 3:
    sys.exit(f'expected 3 latex blocks in the sign-off, found {len(blocks)}')
N1, N4, N5 = (b.rstrip('\n') + '\n' for b in blocks)
# D3 knock-on (author 2026-11-25: ratios to two decimals, as already in the SI)
N5 = N5.replace('Quadrants changed little (0.59 to 0.59)',
                'Quadrants changed little (0.59 at both spreads)')

SI_EDITS = {
    'N1': ('\\FloatBarrier\n\\section{Simulation design, resampling and metrics}',
           N1 + '\n\\FloatBarrier\n\\section{Simulation design, resampling and metrics}'),
    'N3': ('coefficient \\(\\tau = 1\\), and \\(\\sigma = 1\\). The spatial multiplier',
           'coefficient \\(\\tau = 1\\), and \\(\\sigma = 1\\). Section~\\ref{incidence-coefficient} '
           'explains how \\(\\tau\\) can be estimated after a trial. The spatial multiplier'),
    'N4': ('\\FloatBarrier\n\\section{Estimation and identification}',
           N4 + '\n\\FloatBarrier\n\\section{Estimation and identification}'),
    'N5': ('\\FloatBarrier\n\\section{Case study: accuracy, inference and reach}',
           '\\FloatBarrier\n' + N5 + '\n\\FloatBarrier\n\\section{Case study: accuracy, inference and reach}'),
}
CTJ_EDITS = {
    'C1': ('It is not estimated at the design stage, and we set $\\tau=1$ in the grid simulation.',
           'It is not estimated at the design stage, and we set $\\tau=1$ in the grid simulation. '
           'After a trial, the same model estimates it, less precisely when treatment follows '
           'incidence (SI Section~\\ref{SI-incidence-coefficient}).'),
    'C2': ('which helps separate the direct effect from neighboring exposure.',
           'which helps separate the direct effect from neighboring exposure. On the grid, '
           'widening the quadrants\' saturation spread beyond 80/20 lowered control-only MSE but '
           'raised both-arms MSE (SI Section~\\ref{SI-igs-sensitivity-results}).'),
}


def apply(path, edits):
    text = path.read_text()
    for key, (old, new) in edits.items():
        n = text.count(old)
        if n != 1:
            sys.exit(f'{key}: anchor matched {n} times in {path.name}; nothing written')
        text = text.replace(old, new)
    return text


si_text = apply(SI, SI_EDITS)
ctj_text = apply(CTJ, CTJ_EDITS)
SI.write_text(si_text)
CTJ.write_text(ctj_text)
print('applied', ', '.join(list(SI_EDITS) + list(CTJ_EDITS)))
