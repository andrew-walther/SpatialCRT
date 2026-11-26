#!/usr/bin/env python3
"""Apply the approved Phase 2i text edits to Chapter 3 (sign-off: 2i-signoff.md).

Author: Andrew Walther. Created: 2026-11-26.
Usage: python3 apply_2i.py [--dry-run] [--out FILE]

Every edit is read from the sign-off's ```qmd (or ````qmd) fences: each "Current" block
is replaced by the "Proposed" block of the same label within its item, so the applied
text is the signed-off text. Every item in the file was accepted (see "Author decisions").
Each Current block must match exactly once in Thesis_Chapter.qmd; otherwise the script
stops and writes nothing. Items are applied in file order against the original text's
non-overlapping anchors.

C3N-10 runs last: it renumbers the appendix sections after A6 (headings and \\applabel
numbers) as its table gives. Each old heading and \\applabel must match exactly once.

Figure, table-in-SI and verifier items (C3F-1 to C3F-7) are not text edits and are
applied separately.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHAPTER = HERE.parents[2] / 'thesis_chapter/Thesis_Chapter.qmd'
SIGNOFF = HERE / '2i-signoff.md'

# C3N-10: label key -> (old number, new number, heading title)
RENUMBER = [
    ('app:supp-comparisons', 'A7', 'A8', 'Pairwise Tests and Critical-Difference Diagrams'),
    ('app:sensitivity-tables', 'A8', 'A9',
     'Results by Incidence Configuration, Spatial Dependence and Spillover Magnitude'),
    ('app:rook-figures', 'A9', 'A10', 'Rook Contiguity and the Bias--Variance Decomposition'),
    ('app:oracle-nonoracle', 'A10', 'A11', 'Sensitivity Analysis Without the Spillover Term'),
    ('app:nc-results', 'A11', 'A13', 'Additional Case Study Results'),
    ('app:nc-ext-results', 'A12', 'A14', 'Case Study Participation, Power and Sensitivity Results'),
]


def read_edits():
    """Return a list of dicts (item, label, old, new) in file order."""
    lines = SIGNOFF.read_text().split('\n')
    edits, item, pending, i = [], None, {}, 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('## '):
            item = None  # batch headings and closing sections
        if line.startswith('### C3'):
            item, pending = line[4:].split('.')[0], {}
        fence = line.strip()
        if fence in ('```qmd', '````qmd') and item:
            close = fence[:-3]
            j = i + 1
            while lines[j].strip() != close:
                j += 1
            label = next(lines[k].strip() for k in range(i - 1, -1, -1) if lines[k].strip())
            block = '\n'.join(lines[i + 1:j])
            m = re.match(r'\*\*(Current|Proposed)(.*)\*\*$', label)
            if m and m.group(1) == 'Current':
                pending[m.group(2)] = block
            elif m and m.group(2) in pending:
                edits.append(dict(item=item, label=m.group(2).strip() or '(single)',
                                  old=pending.pop(m.group(2)), new=block))
            i = j
        i += 1
    return edits


def apply_edits(text, edits):
    """Locate every anchor in the original text, check overlaps, then splice."""
    spans = []
    for e in edits:
        n = text.count(e['old'])
        if n != 1:
            sys.exit(f"{e['item']} {e['label']}: anchor matches {n} times; nothing written")
        k = text.index(e['old'])
        spans.append((k, k + len(e['old']), e))
    spans.sort(key=lambda s: s[0])
    for (a1, b1, e1), (a2, b2, e2) in zip(spans, spans[1:]):
        if a2 < b1:
            sys.exit(f"{e1['item']} and {e2['item']}: anchors overlap; nothing written")
    out, pos = [], 0
    for a, b, e in spans:
        out.append(text[pos:a])
        out.append(e['new'])
        pos = b
    out.append(text[pos:])
    return ''.join(out)


def renumber(text):
    for key, old, new, title in RENUMBER:
        for a, b in [(f'## {old}. {title} ', f'## {new}. {title} '),
                     (f'\\applabel{{{old}}}{{{key}}}', f'\\applabel{{{new}}}{{{key}}}')]:
            if text.count(a) != 1:
                sys.exit(f'C3N-10: "{a}" matches {text.count(a)} times; nothing written')
            text = text.replace(a, b)
    return text


def main():
    edits = read_edits()
    text = renumber(apply_edits(CHAPTER.read_text(), edits))
    items = sorted({e['item'] for e in edits})
    print(f'{len(edits)} edits from {len(items)} items, then C3N-10 on {len(RENUMBER)} sections')
    if '--dry-run' in sys.argv:
        print('dry run: nothing written')
        return
    dest = Path(sys.argv[sys.argv.index('--out') + 1]) if '--out' in sys.argv else CHAPTER
    dest.write_text(text)
    print(f'wrote {dest}')


if __name__ == '__main__':
    main()
