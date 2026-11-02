# AGENTS.md — Manuscript Revisions (Project 1, BMC Medical Research Methodology)

**Status (2026-11-02): final.** The paper was accepted after Revision 2c. This
folder is the archive of the peer-review rounds; nothing here is being edited.
Treat it as read-only unless Andrew asks otherwise.

- **Sent to the journal (accepted version):**
  `Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2_Submission/Revisions_V2c.{tex,pdf}`.
  This is the clean package that was uploaded, and it keeps its own real `figures/`
  folder so it can be zipped or uploaded to Overleaf as is.
- **Working copy of the same version:** `Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/`.
  Same `.tex`, `.pdf` and `.bib`, plus the annotated diff and the V2b → V2c diff.
- **Post-publication correction (2026-10):** the corrected manuscript, markup and
  memo to BMC are in `../../review/manuscript/`, not here. See `../../review/README.md`.

For the project as a whole see `../../AGENTS.md`; for the `paper/` layout see
`../README.md`.

---

## Rounds

Each round folder holds the same kinds of things: the input manuscript
(`Original & Annotated Manuscripts/`), the revised LaTeX package, the journal's
comments, and the response letter.

| Folder | Responds to | Manuscript | Response letter | Notes |
|---|---|---|---|---|
| `Revision 1/` | Round 1 reviews (editor, two reviewers; `.docx`) | `Revisions_V1.tex` (22 changes from `main.tex`, the initial submission) | `response_letter_V3` (earlier drafts in `revision response/archive/`) | Planning documents in `revision response/planning/` |
| `Revision 2/` | Round 2 reviews (`BMC MRM review spatial 2.pdf`, editor comments) | `Revisions_V2.tex` | `response_letter_V4` | Working notes in `notes/` (`.qmd`/`.html`) |
| `Revision 2a/` | Round 2 (same comments) | `Revisions_V2.tex` plus a shortened `Revisions_short_V2.tex` | `response_letter_V4` | Length-reduction attempt |
| `Revision 2b/` | Round 2 (same comments) | `Revisions_V2b.tex` and `_V2b_short/` | `response_letter_V5` and `_V5_short` | `Notes/` holds six session handoffs and the letter-build scripts; `Diffs/` holds the change summaries |
| `Revision 2c/` | 6/19/2026 meeting with Dr. Lin (notation, ICC statement) | `Revisions_V2c.tex` (accepted) | `response_letter_V6` | Plan in `Notes/SpatialCRT_manuscript_updates_plan_6_21_26.md` |
| `Send to Dr. Lin/` | — | PDFs sent to Dr. Lin before the V2b and V2c submissions | — | See its `README.md` |

The annotated diffs (`Revisions_Annotated_V*.tex`) are all baselined against
`Revisions_V1.tex`.

## Figures

The figures did not change across the revision rounds. Every round's `figures/`
folder is a relative symlink to `../figures/revision/` (in `paper/`), except
`_V2_Submission/figures/`, which is a real copy. The initial-submission packages
in `../Walther_Submission/` link to `../figures/submission/`. Each `.tex` still
uses `figures/...` paths, so every round compiles without edits.

## Rebuilding a round

From inside a round's LaTeX folder:

```bash
pdflatex -interaction=nonstopmode Revisions_V2c.tex
bibtex Revisions_V2c
pdflatex -interaction=nonstopmode Revisions_V2c.tex
pdflatex -interaction=nonstopmode Revisions_V2c.tex
```

Build by-products (`.aux`, `.log`, `.out`, ...) are git-ignored by `../.gitignore`;
`.bbl` files are tracked because journal LaTeX submissions use them. To avoid
overwriting an archived PDF, build into another folder with
`-output-directory=<dir>` (run `bibtex` there with `BIBINPUTS` set to the round's
folder).

## Gotchas

- Two `LaTeX Error: Command \float@end already defined.` messages appear in every
  compile (a `float`/`algorithm` package clash). They are harmless.
- The initial-submission `main.tex` (`../Walther_Submission/`) needs `linegoal.sty`,
  which the local TinyTeX does not have (`tlmgr install linegoal`).
- "Checkerboard paradox", which appears in the Revision 1–2c text (R2.M4), was
  retired on 2026-10-22. The pattern came from wrong 3×3/3×4 BSS allocations; the
  real effect is collinearity (z = 1 − x when every control is exposed, so
  β̂ ≈ β − ψ). See `../../review/COMPARISON.md`.
- BSS allocation counts: 2 / 6 / 2 on 2×4 / 3×3 / 3×4 (3×3 uses the relaxed rule,
  no two intervention districts adjacent). SRS: 70 / 126 / 924.
