# review/manuscript/ — corrected Project 1 manuscript and correction memo

Built 2026-10-25 from the accepted source
(`paper/Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/`).
Nothing has been sent to the journal.

| File | What |
|---|---|
| `Revisions_V2c_original.tex` | Unchanged copy of the accepted `Revisions_V2c.tex` (input to the diff) |
| `Revisions_V2c_corrected.tex` / `.pdf` | Corrected manuscript: 3×3 and 3×4 tables and figures replaced, text edits T1–T13 |
| `Revisions_V2c_corrections_markup.tex` / `.pdf` | latexdiff of original vs corrected; tables shown with the new values only |
| `figures/` | Copy of the accepted figures, with `BetaMSEresults/MSE_Estimates_3x3.png` and `_3x4.png` replaced by the corrected ones |
| `correction_memo.md` / `.pdf` | Letter to the journal plus proposed Correction-note text |
| `memo_header.tex` | Render-time header that lets the memo keep the characters ×, ψ, ρ |

## Rebuild

From `review/` (TeX is TinyTeX: `export PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$PATH"`):

```bash
Rscript code/build_tables_figures.R
python3 code/text_corrections.py tex manuscript/Revisions_V2c_original.tex /tmp/v2c_text.tex
python3 code/make_latex_tables.py /tmp/v2c_text.tex manuscript/Revisions_V2c_corrected.tex corrected
cp figures/paper_size/MSE_Estimates_3x3_corrected.png manuscript/figures/BetaMSEresults/MSE_Estimates_3x3.png
cp figures/paper_size/MSE_Estimates_3x4_corrected.png manuscript/figures/BetaMSEresults/MSE_Estimates_3x4.png
cd manuscript
latexmk -pdf Revisions_V2c_corrected.tex
latexdiff --type=UNDERLINE --config="PICTUREENV=(?:picture|DIFnomarkup|tabular)[\w\d*@]*" \
  Revisions_V2c_original.tex Revisions_V2c_corrected.tex > Revisions_V2c_corrections_markup.tex
latexmk -pdf Revisions_V2c_corrections_markup.tex
latexmk -c
pandoc correction_memo.md -o correction_memo.pdf --pdf-engine=pdflatex -H memo_header.tex \
  -V geometry:margin=1in -V fontsize=11pt -V colorlinks=false
```

(`pandoc` here is RStudio's bundled one:
`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc`.)

## Checks

- `make_latex_tables.py` in `published` mode reproduces the accepted 3×3 and
  3×4 tables except three cells that were typed slightly off in the paper
  (3×3 Int. & Control BSS: MSE 0.053 printed 0.054, SD×10 0.010 printed
  0.007; 3×4 Control Only BSS MSE×10 0.434 printed 0.435). All three are
  replaced by the correction.
- Every number in the new text matches a cell of
  `review/tables/beta_mse_table_<grid>_corrected.csv` or the optimal-allocation
  counts (`review/tables/optimal_combos_<grid>_corrected.csv`).
- The corrected PDF compiles with no undefined references, 46 pages like the
  accepted PDF, and contains no "paradox".
