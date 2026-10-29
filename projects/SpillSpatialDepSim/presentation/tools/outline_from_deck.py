"""Rebuild the slide entries of presentation/OUTLINE.md from the rendered deck source.

Script: outline_from_deck.py
Purpose: Keep OUTLINE.md's slide-by-slide section in step with spillover_talk.qmd, so
         the outline describes the deck as built, with the computed numbers filled in.
Author: Andrew Walther
Created: 2026-10-28
Dependencies: Python 3 standard library only

Run from the project root, after rendering the deck:
  python3 presentation/tools/outline_from_deck.py
  python3 presentation/tools/build_outline_html.py
Adapted from IncidenceDesign/presentation/tools/outline_from_deck.py. It reads
spillover_talk.revealjs.md, the knitted copy Quarto keeps (keep-md: true), so inline R
values appear as numbers. Everything in OUTLINE.md above "<!-- SLIDES (generated) -->"
is kept as written. Below it, one "### n. title" entry per slide is written with Figure,
Slide text, Notes, Time, Carry-over and any Prelim or Conference note from the dicts
below. A slide missing from CARRY is reported, not silently tagged.
"""
import html
import re
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
MD = HERE / 'spillover_talk.revealjs.md'
OUT = HERE / 'OUTLINE.md'
MARK = '<!-- SLIDES (generated) -->'

# Carry-over tags (keyed by data-menu-title): both = same slide in the prelim and the
# conference talk; prelim = prelim only; conference = standalone/conference only (the
# prelim drops it).
CARRY = {
    'Title': 'conference', 'Outline': 'conference', 'Trials of places': 'both', 'Spillover': 'both',
    'Spatial dependence': 'both', 'The design gap': 'both', 'The question': 'both', 'Neighbors': 'both',
    'The two designs': 'both', 'Two spillover types': 'both', 'Model': 'both', 'Simulation design': 'both',
    'Where BSS falls': 'both', 'Control clusters only': 'both', 'Control & intervention clusters': 'both',
    'Why: collinearity': 'both', 'Larger grids': 'both', 'NC probation districts': 'both',
    'NC MCO areas': 'both', 'Recommendation': 'both', 'Choosing a design': 'both', 'Take-home': 'both',
    'Limitations and next steps': 'both', 'Thank you & questions': 'conference',
}
PRELIM_NOTE = {
    'Title': 'Dropped; the prelim has its own title slide and opens Part I with its divider.',
    'Outline': 'Dropped; the prelim has one agenda for all three parts.',
    'Trials of places': 'Serves as the prelim\'s shared motivation for cluster trials of places. Part II builds on it.',
    'Spillover': 'Shared motivation; Part II refers back to it rather than repeating the examples.',
    'Spatial dependence': 'Shared motivation; Part II refers back to it.',
    'The design gap': 'Part I\'s own prior work and gap. Part II adds heterogeneous incidence to it.',
    'Model': 'Part I defines its own model and notation in full; Part II defines its own (no harmonizing).',
    'Limitations and next steps': 'Keeps the paper reference. Add one closing sentence pointing to Part II (larger grids, heterogeneous incidence, queen neighbors, more designs).',
    'Thank you & questions': 'Dropped; the acknowledgements and funding move to the single prelim closing slide after Part III, grouped with Parts II and III.',
}
CONF_NOTE = {
    'Trials of places': 'Cut first for a statistics audience in a 10–15 minute slot (keep one sentence on the next slide).',
    'Spatial dependence': 'In a 10–15 minute slot, fold into one sentence on "Spillover".',
    'The design gap': 'In a 15–20 minute slot, keep; in a 10–15 minute slot, reduce to the bottom card.',
    'Neighbors': 'In a 10–15 minute slot, cut and say "rook neighbors" on the designs slide.',
    'Model': 'For a statistics audience, can add the reduced form from backup.',
    'Control & intervention clusters': 'In a 10–15 minute slot, cut; "Where BSS falls" and "Why: collinearity" carry the point.',
    'Larger grids': 'In a 10–15 minute slot, cut to one sentence on "Choosing a design".',
    'NC MCO areas': 'In a 15–20 minute slot, can merge into "NC probation districts"; cut in a 10–15 minute slot.',
    'Take-home': 'Keep in every length.',
}


def visible_text(block):
    """Strip notes, HTML tags, fenced-div markers and markdown emphasis from a slide body."""
    block = re.sub(r'::: \{\.notes\}.*?\n:::', '', block, flags=re.S)
    block = re.sub(r'<svg.*?</svg>', ' ', block, flags=re.S)
    block = re.sub(r'<!--.*?-->', ' ', block, flags=re.S)
    block = re.sub(r'<table.*?</table>', ' [table] ', block, flags=re.S)
    block = re.sub(r'<[^>]+>', ' ', block)
    block = re.sub(r'^:::.*$', ' ', block, flags=re.M)
    block = re.sub(r'^\|[-:| ]+\|$', ' ', block, flags=re.M)
    block = re.sub(r'^#+\s*', '', block, flags=re.M)
    block = block.replace('**', '').replace('&middot;', '·').replace('&nbsp;', ' ')
    block = html.unescape(re.sub(r'^\s*(- |\d+\. )', '', block, flags=re.M))
    return re.sub(r'\s+', ' ', block).strip()


def main():
    """Parse the knitted deck into slides and rewrite the generated part of OUTLINE.md."""
    src = MD.read_text()
    slides = re.split(r'(?m)^(?=## )', src)[1:]
    entries, section, n_main, n_back, missing, total = [], None, 0, 0, [], 0.0
    for sl in slides:
        head = sl.split('\n', 1)[0]
        title = re.sub(r'\s*\{.*\}\s*$', '', head[3:]).replace('Backup &middot; ', '').strip() or 'Title'
        menu = html.unescape(re.search(r'data-menu-title="([^"]+)"', head).group(1))
        part = re.search(r'data-footer-part="([^"]*)"', head).group(1)
        backup = menu.startswith('B:')
        part = part or ('Background' if menu in ('Title', 'Outline') else 'Discussion')
        if part != section:
            section = part
            n_sec = len([e for e in entries if e.startswith('\n## ')]) + 1
            entries.append('\n## ' + ('Backup' if backup else f'{n_sec}. {part}') + '\n')
        figs = re.findall(r'src="figures/([^"]+)"', sl)
        notes = re.search(r'::: \{\.notes\}\n(.*?)\n:::', sl, flags=re.S)
        notes = notes.group(1) if notes else ''
        tm = re.match(r'\[([\d.]+) min', notes.strip())
        body = re.sub(r'^\[[^\]]*\]\s*', '', notes.strip())
        sents = [x[2:] if x.startswith('- ') else x for x in body.split('\n') if x.strip()]
        note = ' '.join(sents[:2]) if not backup else ' '.join(sents)
        text = visible_text(sl.split('\n', 1)[1] if '\n' in sl else '')
        if len(text) > 420:
            text = text[:417].rsplit(' ', 1)[0] + ' …'
        if backup:
            n_back += 1
            label = f'B{n_back}'
        else:
            n_main += 1
            label = str(n_main)
            if tm and CARRY.get(menu) != 'conference':
                total += float(tm.group(1))
        tag = 'both' if backup else CARRY.get(menu)
        if tag is None:
            missing.append(menu)
        lines = [f'### {label}. {title}',
                 f'- Figure: {", ".join(figs) if figs else "none"}',
                 f'- Slide text: {text or "(title or image only)"}',
                 f'- Notes: {note}']
        if tm:
            lines.append(f'- Time: {tm.group(1)}')
        lines.append(f'- Carry-over: {tag or "UNTAGGED"}')
        if menu in PRELIM_NOTE:
            lines.append(f'- Prelim note: {PRELIM_NOTE[menu]}')
        if menu in CONF_NOTE:
            lines.append(f'- Conference note: {CONF_NOTE[menu]}')
        entries.append('\n'.join(lines) + '\n')
    pre = OUT.read_text().split(MARK)[0].rstrip() + '\n\n' + MARK + '\n'
    OUT.write_text(pre + '\n'.join(entries))
    print(f'Wrote {OUT}: {n_main} main, {n_back} backup slides; '
          f'{total:.2f} min on slides kept in the prelim')
    for m in missing:
        print('WARNING: no carry-over tag for', m)


if __name__ == '__main__':
    main()
