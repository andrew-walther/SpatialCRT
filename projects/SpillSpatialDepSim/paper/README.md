# Manuscript — Spatial Cluster Randomized Trials: Sampling Design with Spillover Effects & Spatial Dependence

**Authors:** Andrew Walther, Tonya Van Deinse, Feng-Chang Lin (University of North Carolina at Chapel Hill)

**Journal:** *BMC Medical Research Methodology*

> ## Final accepted manuscript
> **[`Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2_Submission/Revisions_V2c.pdf`](https://github.com/andrew-walther/SpatialCRT/blob/main/projects/SpillSpatialDepSim/paper/Manuscript%20Revisions/Revision%202c/Walther_SpatialCRT_LaTeX_Revisions_V2_Submission/Revisions_V2c.pdf)**
>
> `_V2_Submission/` is the clean package that was sent to the journal. The
> `Walther_SpatialCRT_LaTeX_Revisions_V2c/` folder next to it is the working copy:
> byte-identical `.pdf`/`.tex`/`.bib`, plus the annotated diff and the V2b → V2c diff.
> Either is the accepted manuscript.
>
> A post-publication correction of the 3×3 and 3×4 results (2026-10) lives in
> `../review/manuscript/`, not here.

## Directory map

| Folder | Contents |
|--------|----------|
| `Cover Letter/` | Submission cover letter (+ generic templates used to draft it) |
| `Manuscript Drafts/` | Pre-submission draft PDFs, V0–V11 |
| `Manuscript Edits/` | Annotated/tracked-changes passes on early drafts |
| `figures/` | The two figure sets the LaTeX packages use: `submission/` (initial submission) and `revision/` (unchanged from Revision 1 through 2c). Each package's `figures/` is a relative symlink to one of these, except `_V2_Submission/figures/`, which is a real copy |
| `Manuscript Figures/` | Source figure images used across drafts |
| `Manuscript Notes/` | Working notes, meeting photos, discussion drafts |
| `Manuscript Revisions/` | All peer-review rounds (Revision 1 → 2 → 2a → 2b → 2c), each with its LaTeX source, compiled PDF, `journal comments/`, and `revision response/`; its `AGENTS.md` maps the rounds |
| `Publication/` | Research Square preprint (posted during initial review, Feb 5, 2026; reflects the manuscript as initially submitted, not the accepted revisions) |
| `Walther_Submission/` | Initial submission package: LaTeX source + compiled PDF (two versions), cover letter, and `Submission_Figures/` (the numbered figure set uploaded with it) |

This folder holds copies, not a build pipeline — everything above is a static PDF/DOCX
snapshot from when it was produced. A full copy of the folder as it was before cleanup is
kept in the author's OneDrive.

**Cleanup (2026-11-02).** Duplicate figure folders became symlinks to `figures/`;
`Manuscript Revisions/Revision 2c/_archive_2b/` (every file also in `Revision 2b/`), the
top-level `Submission_Figures/` (a renamed copy of `Walther_Submission/Submission_Figures/`
plus an `Obsolete/` set that is in `Manuscript Figures/`) and byte-identical copies in
`Send to Dr. Lin/` were removed; LaTeX/R build files are no longer tracked (see
`.gitignore`; `.bbl` files stay tracked).

## Rebuilding the final manuscript

If you edit `Revisions_V2c.tex`, rebuild the PDF from inside
`Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2_Submission/`
(or the `_V2c/` copy — see note above):

```bash
pdflatex -interaction=nonstopmode Revisions_V2c.tex
bibtex Revisions_V2c
pdflatex -interaction=nonstopmode Revisions_V2c.tex
pdflatex -interaction=nonstopmode Revisions_V2c.tex
```

`references.bib` and the `figures/` folder are already alongside `Revisions_V2c.tex`,
so no other paths need to be set up (in `_V2c/` and older rounds `figures/` is a symlink,
which works locally but not in a zip uploaded to Overleaf; use `_V2_Submission/` for that). Requires a standard TeX distribution (TinyTeX,
MacTeX, TeX Live) — no packages beyond `article`-class defaults (`amsmath`, `booktabs`,
`graphicx`, `hyperref`, `tikz`, etc., all standard).
