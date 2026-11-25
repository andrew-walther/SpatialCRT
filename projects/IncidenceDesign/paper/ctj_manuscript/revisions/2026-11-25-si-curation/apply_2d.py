#!/usr/bin/env python3
"""Apply the approved Phase 2d SI table edits (sign-off: 2d-signoff.md).

Author: Andrew Walther. Created: 2026-11-25.
TD1-TD3 copy the Chapter 3 rows of tab:designs, tab:parameter-grid and tab:excluded-designs
(approved earlier there) into the SI; TD4 types the bias -0.007 as $-0.007$ in tab:nc-annual;
TD5 renames tab:gamma to tab:psi. Each row is taken from Thesis_Chapter.qmd, not retyped, and
must replace exactly one SI row whose first cell matches.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'paper/tools'))
import verify_manuscripts as v  # table() extracts a labeled tabular body

SI = v.SUPPLEMENT
ROWS = {'tab:designs': ['Incidence-Guided Saturation', 'Isolation Buffer', '2$\\times$2 Blocking', 'Checkerboard'],
        'tab:parameter-grid': ['Incidence configuration', 'Poisson mechanism'],
        'tab:excluded-designs': ['Saturation Quadrants', 'Balanced Halves']}


def row(body, first_cell):
    """Return the single table row of body that starts with first_cell."""
    hits = [l for l in body.split('\n') if l.strip().startswith(first_cell + ' &')]
    if len(hits) != 1:
        raise SystemExit(f'{first_cell}: expected 1 row, found {len(hits)}')
    return hits[0].strip()


def main():
    chapter = v.CHAPTER.read_text()
    si = SI.read_text()
    for label, cells in ROWS.items():                       # TD1-TD3
        for cell in cells:
            old, new = row(v.table(si, label), cell), row(v.table(chapter, label), cell)
            if si.count(old) != 1:
                raise SystemExit(f'{label}/{cell}: SI row not unique')
            si = si.replace(old, new)
    body = v.table(si, 'tab:nc-annual')                     # TD4
    fixed = re.sub(r'(?<![$\w])-0\.007\b', '$-0.007$', body)
    if fixed.count('$-0.007$') - body.count('$-0.007$') != 4:
        raise SystemExit('TD4: expected four -0.007 cells')
    si = si.replace(body, fixed)
    for old, new in [('\\label{tab:gamma}', '\\label{tab:psi}'), ('\\ref{tab:gamma}', '\\ref{tab:psi}')]:
        if si.count(old) != 1:                               # TD5
            raise SystemExit(f'TD5: {old} found {si.count(old)} times')
        si = si.replace(old, new)
    SI.write_text(si)
    print('applied TD1-TD5')


if __name__ == '__main__':
    main()
