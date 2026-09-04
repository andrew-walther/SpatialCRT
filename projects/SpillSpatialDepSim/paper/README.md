# Manuscript — Spatial Cluster Randomized Trials: Sampling Design with Spillover Effects & Spatial Dependence

**Authors:** Andrew Walther, Tonya Van Deinse, Feng-Chang Lin (University of North Carolina at Chapel Hill)

**Journal:** *BMC Medical Research Methodology*

> ## 📄 Final accepted manuscript
> **[`Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/Revisions_V2c.pdf`](https://github.com/andrew-walther/SpatialCRT/blob/main/projects/SpillSpatialDepSim/paper/Manuscript%20Revisions/Revision%202c/Walther_SpatialCRT_LaTeX_Revisions_V2c/Revisions_V2c.pdf)**

## Directory map

| Folder | Contents |
|--------|----------|
| `Cover Letter/` | Submission cover letter (+ generic templates used to draft it) |
| `Manuscript Drafts/` | Pre-submission draft PDFs, V0–V11 |
| `Manuscript Edits/` | Annotated/tracked-changes passes on early drafts |
| `Manuscript Figures/` | Source figure images used across drafts |
| `Manuscript Notes/` | Working notes, meeting photos, discussion drafts |
| `Manuscript Revisions/` | All peer-review rounds (Revision 1 → 2 → 2a → 2b → 2c), each with its LaTeX source, compiled PDF, `journal comments/`, and `revision response/` |
| `Publication/` | Research Square preprint (posted during initial review, Feb 5, 2026; reflects the manuscript as initially submitted, not the accepted revisions) |
| `Submission_Figures/` | Numbered figure set matching the initial submission's figure references |
| `Walther_Submission/` | Initial submission package: LaTeX source + compiled PDF (two versions) + figures |

This folder holds copies, not a build pipeline — everything above is a static PDF/DOCX
snapshot from when it was produced. Some content is redundant (e.g., duplicate figure
sets across submission rounds); a cleanup pass is planned but not yet done.

## Rebuilding the final manuscript

If you edit `Revisions_V2c.tex`, rebuild the PDF from inside
`Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/`:

```bash
pdflatex -interaction=nonstopmode Revisions_V2c.tex
bibtex Revisions_V2c
pdflatex -interaction=nonstopmode Revisions_V2c.tex
pdflatex -interaction=nonstopmode Revisions_V2c.tex
```

`references.bib` and the `figures/` folder are already alongside `Revisions_V2c.tex`,
so no other paths need to be set up. Requires a standard TeX distribution (TinyTeX,
MacTeX, TeX Live) — no packages beyond `article`-class defaults (`amsmath`, `booktabs`,
`graphicx`, `hyperref`, `tikz`, etc., all standard).
