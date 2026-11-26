#!/usr/bin/env bash
# ============================================================
# Script: make_diff.sh
# Purpose: Build the marked-up diff (reading layout) from the draft Dr. Lin annotated
#          (tag ctj-pre-fcl) to the current CTJ manuscript.
# Author: Andrew Walther
# Created: 2026-11-25
# Dependencies: TinyTeX (pdflatex, bibtex, latexdiff), python3, git
# ============================================================
# latexdiff cannot mark changes inside the exhibit macro definitions or the abstract
# macro, so inline_exhibits.py expands them at their call sites first. latexdiff also
# leaves the reading-layout title block unmarked, so the marked abstract paragraphs of
# the submission branch are copied into it. Math is marked whole (--math-markup=whole).
# Run from the project root:  bash paper/ctj_manuscript/revisions/2026-11-24-round6-word-limit/for-lin/tools/make_diff.sh <out_dir>
set -euo pipefail
OUT=${1:?output directory}
ROOT=$(pwd); M=$ROOT/paper/ctj_manuscript/manuscript_ctj; T=$(dirname "$0")
mkdir -p "$OUT/manuscript_ctj"; ln -sfn "$ROOT/paper/ctj_manuscript/figures" "$OUT/figures"
# The CTJ prints SI numbers through xr-hyper (../supplementary_info/.build/*.aux); build the SI first
ln -sfn "$ROOT/paper/ctj_manuscript/supplementary_info" "$OUT/supplementary_info"
for f in CTJ_Manuscript.tex CTJ_Exhibits.tex; do
  git show ctj-pre-fcl:projects/IncidenceDesign/paper/ctj_manuscript/manuscript_ctj/$f > "$OUT/old_$f"
done
python3 "$T/inline_exhibits.py" "$OUT/old_CTJ_Manuscript.tex" "$OUT/old_CTJ_Exhibits.tex" "$OUT/old_inl.tex"
python3 "$T/inline_exhibits.py" "$M/CTJ_Manuscript.tex" "$M/CTJ_Exhibits.tex" "$OUT/new_inl.tex"
D=$OUT/manuscript_ctj; cd "$D"
ln -sf "$M/sagej.cls" .; ln -sf "$M/SageV.bst" .; ln -sf "$ROOT/paper/SpatialCRT_IncidenceDesign.bib" .
cp "$M/CTJ_Reading.tex" .
latexdiff --type=UNDERLINE --math-markup=whole "$OUT/old_inl.tex" "$OUT/new_inl.tex" > CTJ_Manuscript.tex
python3 - CTJ_Manuscript.tex <<'PY'
import sys
p=sys.argv[1]; L=open(p).read().split('\n')
a=L.index('\\begin{abstract}'); b=L.index('\\end{abstract}')
marked=[l for l in L[a+1:b] if l.startswith('\\textbf{')]
i=[k for k,l in enumerate(L) if l.startswith('\\textbf{Background/Aims:}')][0]
for h,m in zip(['\\textbf{Background/Aims:}','\\textbf{Methods:}','\\textbf{Results:}','\\textbf{Conclusions:}'],marked):
    j=next(k for k in range(i,a) if L[k].startswith(h)); L[j]=m
open(p,'w').write('\n'.join(L))
PY
for i in 1 2; do pdflatex -interaction=nonstopmode CTJ_Reading.tex >/dev/null; bibtex CTJ_Reading >/dev/null; done
pdflatex -interaction=nonstopmode CTJ_Reading.tex >/dev/null
echo "Diff PDF: $D/CTJ_Reading.pdf"
