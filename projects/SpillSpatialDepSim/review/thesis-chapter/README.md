# review/thesis-chapter/ — corrected Chapter 2 (not yet integrated)

The corrected Chapter 2, ready to copy into bios-dissertation once Andrew
approves.

| Path | What | Copy back? |
|---|---|---|
| `draft/project1-spillover-draft.qmd` / `.tex` / `.pdf` | Corrected chapter | yes |
| `draft/figures/` | Chapter figures; `BetaMSEresults/MSE_Estimates_3x3.png` and `_3x4.png` replaced | the two replaced PNGs |
| `manuscript-declarations.md` | Verbatim copy of the paper's abstract and Declarations, with the abstract edits (T13) and the funder name fixed (T14) | yes |
| `references.bib` | Copy of bios-dissertation `prelim/references.bib` (2026-10-25), used only to render here | no |
| `draft/bios-prelim.cls` | Symlink to bios-dissertation's class file, used only to render here | no |

Built 2026-10-25 from the bios-dissertation copy of the chapter. Changes: the
same text edits as the corrected manuscript (`../code/text_corrections.py`,
edits T1–T12), the corrected 3×3 and 3×4 tables
(`../code/make_latex_tables.py`) and figures, and two sentences (edit N1) at
the end of the paragraph that introduces the results, saying the 3×3 and 3×4
results were recomputed after acceptance. It renders to 52 pages, as before,
with no unresolved cross-references and the same reference list.

The `.qmd` keeps its bios-dissertation paths (`bibliography:
../../../references.bib` and the `apa.csl` path), so it can be copied back
unchanged. To render it here, override those two paths:

```bash
cd draft
PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$PATH" \
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render project1-spillover-draft.qmd \
  -M bibliography:../references.bib \
  -M csl:/Users/ajwalther/GithubProjects/bios-dissertation/templates/source/dissertation-resources-amber/apa.csl
```

## Integrating into bios-dissertation (only after Andrew approves)

```bash
SRC=~/GithubProjects/SpatialCRT/projects/SpillSpatialDepSim/review/thesis-chapter
DST=~/GithubProjects/bios-dissertation/prelim/project-proposals/project1-spillover
cp $SRC/draft/project1-spillover-draft.{qmd,tex,pdf} $DST/draft/
cp $SRC/draft/figures/BetaMSEresults/MSE_Estimates_3x{3,4}.png $DST/draft/figures/BetaMSEresults/
cp $SRC/manuscript-declarations.md $DST/
```

Then, in bios-dissertation:

1. Re-render the chapter there (its own paths now apply) and check it is
   still 52 pages.
2. Rebuild the full prelim (`python3 prelim/full-draft/tools/build_prelim.py`).
3. `prelim/project-proposals/notes/check-prose-fidelity.py` will now report
   the intended differences from `Revisions_V2c.tex`. The exception is
   recorded in `DECISIONS.md` (2026-10-25).
4. Update `ROADMAP.md` ("Chapter 2" status) and
   `prelim/project-proposals/project1-spillover/README.md` (add a line on the
   correction).
