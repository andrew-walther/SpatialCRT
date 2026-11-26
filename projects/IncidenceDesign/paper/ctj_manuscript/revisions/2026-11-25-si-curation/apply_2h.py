#!/usr/bin/env python3
"""Apply the approved Phase 2h edits (sign-off: 2h-signoff.md).

Author: Andrew Walther. Created: 2026-11-26.
Usage: python3 apply_2h.py [--dry-run] [--out DIR]

Everything applied is read from the sign-off, so the applied text is the signed-off text.
- Batch NY (SI) and Batch E/NYC (CTJ): each ```latex "Current" block is replaced by the
  "Proposed" block of the same label within its item.
- Batch PC (SI): pointer edits as above, and each "Moves to the unabridged build" line
  wraps the span from its first words to its last words in \\iffullsi ... \\fi. When the
  last words run into the next block (a following \\iffullsi figure or subsection, or the
  sentence "To check that"), that tail stays outside the wrap.
- Batch SC (SI): SC-1 adds the short title from its table to each caption
  (\\caption[short]{full}); a label listed twice gets its titles in source order. The four
  subfigure captions keep their own text. SC-2 was declined and is not applied.

Order: NY, then PC in sign-off order, then SC-1; NYC on the CTJ. Each anchor must match
exactly once (with any run of whitespace allowed between its words when there is no exact
match); otherwise the script stops and writes nothing. The title-page word count is set
after this script, from ctj_wordcount.py.

--out DIR writes the results to DIR instead of the manuscript files (for checking).
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SI = PAPER / 'supplementary_info/Supplementary_Information.tex'
CTJ = PAPER / 'manuscript_ctj/CTJ_Manuscript.tex'
SIGNOFF = HERE / '2h-signoff.md'

SKIP_ITEMS = {'SC-2'}  # declined 2026-11-26: each list keeps its own page
SUBFIGURES = {'fig:cd-queen', 'fig:cd-rook', 'fig:tau-coverage-queen', 'fig:tau-coverage-rook'}
WRAP_TAILS = ['\\iffullsi\n\\begin{figure}', 'To check that',
              '\\iffullsi\n\\subsection{Implications for trial planning}']


# Read the sign-off ----

def read_signoff():
    """Return (ops, short_titles).

    ops: list of dicts in file order, kind 'replace' (item, label, target, old, new) or
    'wrap' (item, start, end). short_titles: {label: [title, ...]} from SC-1's table.
    """
    lines = SIGNOFF.read_text().split('\n')
    ops, titles, item, target, pending, i = [], {}, None, 'si', {}, 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('## Batch E'):
            target = 'ctj'
        elif line.startswith('## ') and not line.startswith('## Batch'):
            target = 'si'
            item = None  # summary, follow-up and decision sections hold no edits
        if line.startswith('### '):
            item, pending = line[4:].split('.')[0], {}
        if item == 'SC-1':
            m = re.match(r'\| `([^`]+)` \| (.+) \|$', line)
            if m:
                titles.setdefault(m.group(1), []).append(m.group(2))
        m = re.match(r'- \*\*Moves to the unabridged build:\*\* from "(.*)" to "(.*)"\.$', line)
        if m and item and item not in SKIP_ITEMS:
            ops.append(dict(kind='wrap', item=item, start=m.group(1), end=m.group(2)))
        if line.strip() == '```latex':
            j = i + 1
            while lines[j].strip() != '```':
                j += 1
            label = next(lines[k].strip() for k in range(i - 1, -1, -1) if lines[k].strip())
            block = '\n'.join(lines[i + 1:j])
            m = re.match(r'\*\*(Current|Proposed)(.*)\*\*$', label)
            if m and m.group(1) == 'Current':
                pending[m.group(2)] = block
            elif m and m.group(2) in pending and item and item not in SKIP_ITEMS:
                ops.append(dict(kind='replace', item=item, label=m.group(2).strip() or '(single)',
                                target=target, old=pending.pop(m.group(2)), new=block))
            i = j
        i += 1
    return ops, titles


# Matching ----

def flexible(s):
    """Regex for s with any whitespace run allowed between its tokens."""
    return re.compile(r'\s+'.join(re.escape(tok) for tok in s.split()))


def find_once(text, s, what, pos=0):
    """(start, end) of the single match of s in text[pos:]; exact first, then flexible."""
    n = text.count(s, pos)
    if n == 1:
        k = text.find(s, pos)
        return k, k + len(s)
    if n == 0:
        hits = list(flexible(s).finditer(text, pos))
        if len(hits) == 1:
            return hits[0].start(), hits[0].end()
        n = len(hits)
    sys.exit(f'{what}: anchor matches {n} times; nothing written')


def find_next(text, s, what, pos):
    """(start, end) of the first match of s at or after pos (exact, then flexible)."""
    k = text.find(s, pos)
    if k >= 0:
        return k, k + len(s)
    m = flexible(s).search(text, pos)
    if m:
        return m.start(), m.end()
    sys.exit(f'{what}: end anchor not found; nothing written')


# Apply ----

def apply_op(text, op):
    name = f"{op['item']} {op.get('label', 'wrap')}"
    if op['kind'] == 'replace':
        a, b = find_once(text, op['old'], name)
        return text[:a] + op['new'] + text[b:]
    i, _ = find_once(text, op['start'], name + ' start')
    _, k = find_next(text, op['end'], name, i)
    span, keep = text[i:k], ''
    for tail in WRAP_TAILS:
        if span.endswith(tail):
            span, keep = span[:-len(tail)], tail
    return text[:i] + '\\iffullsi\n' + span.rstrip('\n') + '\n\\fi\n\n' + keep + text[k:]


def caption_end(t, j):
    """Index of the brace that closes the \\caption{ starting at j."""
    d = 0
    for k in range(j + 8, len(t)):
        if t[k] == '{' and t[k - 1] != '\\':
            d += 1
        elif t[k] == '}' and t[k - 1] != '\\':
            d -= 1
            if d == 0:
                return k
    sys.exit('SC-1: unbalanced caption; nothing written')


def apply_short_titles(t, titles):
    out, i, seen, done = [], 0, {}, 0
    while True:
        j = t.find('\\caption{', i)
        if j < 0:
            out.append(t[i:])
            break
        k = caption_end(t, j)
        m = re.search(r'\\label\{([^}]*)\}', t[k:k + 1500])
        # The regional-ordering table's label sits away from its caption
        lab = m.group(1) if m else 'tab:nc-order'
        if lab in SUBFIGURES:
            out.append(t[i:k + 1])
            i = k + 1
            continue
        if lab not in titles:
            sys.exit(f'SC-1: no short title for {lab}; nothing written')
        v = titles[lab]
        n = seen.get(lab, 0)
        seen[lab] = n + 1
        out.append(t[i:j] + '\\caption[' + v[min(n, len(v) - 1)] + ']' + t[j + 8:k + 1])
        i = k + 1
        done += 1
    unused = set(titles) - set(seen)
    if unused:
        sys.exit(f'SC-1: titles never used: {sorted(unused)}; nothing written')
    return ''.join(out), done


def main():
    dry = '--dry-run' in sys.argv
    out = Path(sys.argv[sys.argv.index('--out') + 1]) if '--out' in sys.argv else None
    ops, titles = read_signoff()
    texts = {'si': SI.read_text(), 'ctj': CTJ.read_text()}
    order = ([o for o in ops if o['item'].startswith('NY-')]
             + [o for o in ops if o['item'].startswith('PC-')]
             + [o for o in ops if o['item'].startswith('NYC-')])
    if len(order) != len(ops):
        sys.exit('unexpected item prefix in sign-off; nothing written')
    for op in order:
        tgt = op.get('target', 'si')
        texts[tgt] = apply_op(texts[tgt], op)
    texts['si'], n_sc = apply_short_titles(texts['si'], titles)
    count = lambda p, k: sum(o['item'].startswith(p) and o['kind'] == k for o in ops)
    print(f"NY {count('NY-', 'replace')} edits, PC {count('PC-', 'replace')} edits and "
          f"{count('PC-', 'wrap')} wraps, NYC {count('NYC-', 'replace')} edits, "
          f'SC-1 {n_sc} short titles')
    if dry:
        print('dry run: nothing written')
        return
    si_dest, ctj_dest = (out / SI.name, out / CTJ.name) if out else (SI, CTJ)
    si_dest.write_text(texts['si'])
    ctj_dest.write_text(texts['ctj'])
    print(f'wrote {si_dest} and {ctj_dest}')


if __name__ == '__main__':
    main()
