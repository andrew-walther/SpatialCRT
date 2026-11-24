#!/usr/bin/env python3
"""
Script: ctj_wordcount.py
Purpose: Count CTJ main-text words per section under two rules, so the 3,500-word
         limit can be tracked after every edit.
Author: Andrew Walther
Created: 2026-11-24
Dependencies: Python 3 standard library; texcount (TeX Live)

The body runs from \\section{Introduction} up to the first unnumbered back-matter
section (Acknowledgements); the abstract, declarations, references, captions and
exhibits are excluded. Two counts are reported because the repo has used both:

  text_headers        texcount words in text + headers, math excluded (the
                      submission checklist's rule, 3,498 on 2026-10-27)
  text_headers_math1  the same plus one word per inline or displayed math
                      expression (close to the 2026-11-15 trim plan's rule; that hand
                      count, 3,838, also dropped words texcount keeps, such as \ref
                      targets, so the two do not match exactly)

Neither is confirmed as the journal's own rule; check Clinical Trials' instructions
before submission. Usage (from anywhere):
    python3 paper/tools/ctj_wordcount.py            # table by section and totals
"""
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAIN = ROOT / 'paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex'
LIMIT = 3500


def body_text():
    """Return the main-text LaTeX from Introduction to the first back-matter section."""
    tex = MAIN.read_text()
    start = tex.index('\\section{Introduction}')
    end = tex.index('\\section*{', start)
    return tex[start:end]


def texcount_sections(fragment):
    """Run texcount on a fragment; return [(title, text, headers, math)] per section.

    @param fragment LaTeX source of the body
    @return list of tuples with words in text, words in headers and the number of
            inline plus displayed math expressions for each (sub)section
    """
    with tempfile.NamedTemporaryFile('w', suffix='.tex', delete=False) as handle:
        handle.write(fragment)
        path = handle.name
    out = subprocess.run(['texcount', '-sub=subsection', '-nocol', path],
                         capture_output=True, text=True, check=True).stdout
    Path(path).unlink()
    rows = []
    # Subcount lines look like "  528+6+0 (1/0/53/1) Subsection: Title"
    for m in re.finditer(r'^\s*(\d+)\+(\d+)\+\d+ \(\d+/\d+/(\d+)/(\d+)\) (?:Sub)?[Ss]ection: (.*)$',
                         out, re.M):
        text, headers, inl, disp, title = m.groups()
        rows.append((title.strip(), int(text), int(headers), int(inl) + int(disp)))
    if not rows:
        raise RuntimeError('texcount output had no section subcounts; format changed?')
    return rows


def main():
    """Print words by section under both rules and the distance to the limit."""
    rows = texcount_sections(body_text())
    print(f"{'Section':<70} {'text+hdr':>9} {'+math':>7}")
    tot_a = tot_b = 0
    for title, text, headers, math in rows:
        a, b = text + headers, text + headers + math
        tot_a += a; tot_b += b
        print(f"{title[:70]:<70} {a:>9} {b:>7}")
    print(f"{'TOTAL':<70} {tot_a:>9} {tot_b:>7}")
    print(f"Limit {LIMIT}: over by {tot_a - LIMIT} (text+headers) / {tot_b - LIMIT} (math as one word)")


if __name__ == '__main__':
    main()
