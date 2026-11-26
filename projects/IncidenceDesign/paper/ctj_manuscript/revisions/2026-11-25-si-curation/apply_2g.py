#!/usr/bin/env python3
"""Apply the approved Phase 2g edits (sign-off: 2g-merged-signoff.md).

Author: Andrew Walther. Created: 2026-11-25.
Usage: python3 apply_2g.py [--dry-run]

Every edit is read from the sign-off's ```latex fences, so the applied text is the
signed-off text. Each "Current" block is paired with the "Proposed" block of the same
label within its item. Items under "## Batch E" edit CTJ_Manuscript.tex; all others edit
Supplementary_Information.tex.

Not applied:
- HC-7 to HC-10 (declined by the author);
- HW-6 block 2 (rejected; it has no proposed block);
- HC-20 (the title-page word count is set after this script, from ctj_wordcount.py).

Order. Some current blocks quote the final wording of another item (for example SR-6
builds on HW-10 and HW-27), so the script applies edits in passes: each pass applies
every edit whose anchor matches exactly once in the text as it stands. An anchor that
matches more than once stops the script. Build-on anchors may differ from the applied
text only in line breaks, so an anchor with no exact match is tried with any run of
whitespace allowed between its words. Anything left unapplied after the passes stops the
script, and nothing is written.

SH-1 runs last: it adds the \\sihead command to the SI preamble and turns each bold
run-in head (\\textbf{Head.} at the start of a line) into \\sihead{Head} on its own line.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SI = PAPER / 'supplementary_info/Supplementary_Information.tex'
CTJ = PAPER / 'manuscript_ctj/CTJ_Manuscript.tex'
SIGNOFF = HERE / '2g-merged-signoff.md'

SKIP_ITEMS = {'HC-7', 'HC-8', 'HC-9', 'HC-10', 'HC-20'}
SH1_MACRO = ('% Bold paragraph head on its own line; the paragraph starts on the next line\n'
             '\\newcommand{\\sihead}[1]{\\par\\medskip\\noindent\\textbf{#1}\\par'
             '\\nopagebreak\\noindent\\ignorespaces}\n')
SH1_HEADS = 14  # 8 printed in the curated SI and 6 in the unabridged build only


# Read the sign-off ----

def read_edits():
    """Return a list of dicts (item, label, target, old, new) in file order."""
    lines = SIGNOFF.read_text().split('\n')
    edits, item, target, pending, i = [], None, 'si', {}, 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('## Batch E'):
            target = 'ctj'
        elif line.startswith('## Batch F') or (line.startswith('## ') and not line.startswith('## Batch')):
            target = 'si'
        if line.startswith('### '):
            item, pending = line[4:].split('.')[0], {}
        if line.strip() == '```latex':
            j = i + 1
            while lines[j].strip() != '```':
                j += 1
            label = next(lines[k].strip() for k in range(i - 1, -1, -1) if lines[k].strip())
            block = '\n'.join(lines[i + 1:j])
            m = re.match(r'\*\*(Current|Proposed)(.*)\*\*$', label)
            if m and m.group(1) == 'Current':
                pending[m.group(2)] = block
            elif m and m.group(2) in pending and item not in SKIP_ITEMS:
                edits.append(dict(item=item, label=m.group(2).strip() or '(single)',
                                  target=target, old=pending.pop(m.group(2)), new=block))
            i = j
        i += 1
    return edits


# Apply ----

def flexible(old):
    """Regex for old with any whitespace run allowed between its tokens."""
    return re.compile(r'\s+'.join(re.escape(tok) for tok in old.split()))


def try_apply(text, e):
    """Return (new text, status). Status: 'ok', 'none' or 'multi'."""
    n = text.count(e['old'])
    if n == 1:
        return text.replace(e['old'], e['new']), 'ok'
    if n > 1:
        return text, 'multi'
    hits = list(flexible(e['old']).finditer(text))
    if len(hits) == 1:
        h = hits[0]
        return text[:h.start()] + e['new'] + text[h.end():], 'ok'
    return text, 'multi' if hits else 'none'


def apply_all(texts, edits):
    left, log = list(edits), []
    while left:
        progress, nxt = False, []
        for e in left:
            texts[e['target']], status = try_apply(texts[e['target']], e)
            if status == 'ok':
                progress = True
                log.append(f"{e['item']} {e['label']}")
            elif status == 'multi':
                sys.exit(f"{e['item']} {e['label']}: anchor matches more than once; nothing written")
            else:
                nxt.append(e)
        if not progress:
            names = ', '.join(f"{e['item']} {e['label']}" for e in nxt)
            sys.exit(f'unapplied after all passes: {names}; nothing written')
        left = nxt
    return log


def apply_sh1(si):
    anchor = '\\makeatletter\\newcommand{\\InputRows}[1]{\\@@input #1 }\\makeatother\n'
    if si.count(anchor) != 1:
        sys.exit('SH-1: preamble anchor not found once; nothing written')
    si = si.replace(anchor, anchor + SH1_MACRO)
    si, n = re.subn(r'^\\textbf\{([^{}]+)\.\}[ \t]*', r'\\sihead{\1}\n', si, flags=re.M)
    if n != SH1_HEADS:
        sys.exit(f'SH-1: expected {SH1_HEADS} run-in heads, found {n}; nothing written')
    return si


def main():
    dry = '--dry-run' in sys.argv
    edits = read_edits()
    texts = {'si': SI.read_text(), 'ctj': CTJ.read_text()}
    log = apply_all(texts, edits)
    texts['si'] = apply_sh1(texts['si'])
    print(f'{len(log)} edits applied ({sum(e["target"] == "si" for e in edits)} SI, '
          f'{sum(e["target"] == "ctj" for e in edits)} CTJ), plus SH-1 on {SH1_HEADS} heads')
    if dry:
        print('dry run: nothing written')
        return
    SI.write_text(texts['si'])
    CTJ.write_text(texts['ctj'])
    print(f'wrote {SI.name} and {CTJ.name}')


if __name__ == '__main__':
    main()
