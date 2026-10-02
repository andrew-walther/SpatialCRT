# CTJ paired layouts — initial checkpoint, 2026-10-02

**Later author-approved revision:** original title and six revised exhibits are
now implemented; reading 8 pp, submission 21 pp, chapter 66 pp, SI 44 pp.
Author order is Walther, Habib, Simpson, Lin. The following describes the earlier
layout-only checkpoint. Current details: [exhibit findings](ctj_exhibit_revision_findings_2026-10-02.md).

The author requested both versions to judge the eventual article appearance.
Completed: `CTJ_Manuscript.pdf` (18-page submission draft) and `CTJ_Reading.pdf`
(8-page two-column author preview), built from the same master and exhibits.
This approval covers paired layouts. The proposed restoration of the original
title and six-exhibit revision remains pending; current title and scientific
figure content are retained in this comparison.

## Why the layouts differ

The [live Clinical Trials instructions](https://journals.sagepub.com/author-instructions/ctj)
request double-spaced, minimally formatted review text; tables on separate sheets
at the end; separate figure legends; no more than six main exhibits; and no more
than 3,500 body words/425 structured-abstract words. They do not require authors
to supply two-column publisher proofs. The reading preview uses the existing
Sage `Afour` layout; the submission version uses `Review` with A4 review margins.
The supplied `sagej.cls` is unchanged. No DOI, volume or publication status was
invented in the reading header.

The earlier two-column CTJ source is preserved at
`paper/archive_manuscript/pre_completion_20261002/CTJ_Manuscript.tex`; its PDF
is recoverable from Git before `a47b0ae`. Its title was *Sampling Design for
Spatial Cluster Randomized Trials Under Heterogeneous Incidence*. Its colorful
grid figures and other assets are still present. The older numerical claims and
synthetic application map must not be restored with that visual style.

## Source ownership and builds

- `CTJ_Manuscript.tex`: canonical CTJ prose, abstract and metadata; default
  submission build, conditional reading layout and placement calls.
- `CTJ_Reading.tex`: flag and input only; contains no second manuscript.
- `CTJ_Exhibits.tex`: two table/four figure definitions shared by both layouts.
  Edit cells/captions here. Each is called once in the selected layout.
- `Supplementary_Information.tex`/PDF: unchanged in this layout step.

From `paper/ctj_manuscript/`, using the existing TinyTeX installation:

```sh
PATH=/Users/ajwalther/Library/TinyTeX/bin/universal-darwin:$PATH \
TEXINPUTS=.:../SAGE_Journal_Template: BSTINPUTS=.:../SAGE_Journal_Template: \
latexmk -pdf -interaction=nonstopmode -halt-on-error CTJ_Manuscript.tex

PATH=/Users/ajwalther/Library/TinyTeX/bin/universal-darwin:$PATH \
TEXINPUTS=.:../SAGE_Journal_Template: BSTINPUTS=.:../SAGE_Journal_Template: \
latexmk -pdf -interaction=nonstopmode -halt-on-error CTJ_Reading.tex
```

This is an existing multi-file Sage/bibliography project, compiled with its
existing local workflow; no TeX installation or dependency was added. The
reading wrapper was queued for display in Codex's source editor.

## Verification

- Both latexmk builds succeed. Reading has no final compilation warning;
  neither layout has overfull boxes or undefined references/citations.
  Submission retains existing geometry-override warnings and six underfull
  vertical boxes associated with intentional separate exhibit pages.
- All eight reading pages and 18 submission pages inspected as rendered PNGs,
  with full-size checks of the title/abstract, tables and complex figures.
  No clipping, lost exhibits or unresolved references. Complex graphics span
  both columns; references begin after all reading exhibits.
- Original scientific body and abstract compared with the pre-layout master:
  unchanged after selecting the submission branch and ignoring whitespace.
  All six exhibit labels, figure assets, captions and table values retained.
- `python3 paper/tools/verify_manuscripts.py` passes all existing numerical,
  supplement, figure, reference and citation checks, including the deliberately
  wrong-number fixture. Added checks reject a content-duplicating reading wrapper,
  missing shared definitions, duplicate placement calls or calls on the wrong
  side of the bibliography.
- Fresh texcount: 2,163 prose words plus 21 heading words = **2,184**.
  Earlier 2,217 additionally counted 32 inline and one displayed mathematical
  expression. Abstract **283**. No prose was cut to change this count.
- Existing application outputs, Chapter 3 and its sync artifacts, restricted
  inputs and Project 1 remain unchanged. No simulation or Git push occurred.

## Code/layout teach-back

### Master selection and shared front matter

**Purpose:** produce two useful presentations without two accounts of the study.
**Logic:** the wrapper flag selects `Afour`; absence of the flag selects `Review`.
Shared macros hold title, author line, affiliations, correspondence, keywords
and the verbatim abstract. Review retains its separate title/abstract pages and
line numbers. Reading puts shared front matter above two columns, breaks the
same equation across two lines to fit, and uses a clearly marked preview header.
**Choice:** a small conditional master avoids copied manuscripts and extra
packages. **Inputs/outputs:** the same TeX/bib/figure sources, two named PDFs.
**Gotchas:** build from the CTJ directory with the Sage search paths; do not edit
the class or maintain article prose in the wrapper. **R notes:** no R change.

### Six shared exhibit macros

**Purpose:** keep displayed numbers/captions/assets identical in both versions.
**Logic:** `CTJGridTable` gives the full nine-design queen benchmark;
`CTJAnnualTable` gives year-specific SRS/saturation means; `CTJGridFigure` gives
the detailed six-design grid panels; `CTJMeanFigure`, `CTJRiskFigure` and
`CTJBudgetFigure` give primary NC means, refined risk and geographic reach.
The macros emit starred floats: full-width reading exhibits, ordinary floats
in the one-column version. Reading calls follow relevant results paragraphs;
submission calls follow references with separate-page breaks. A reading flush
before references prevents queued wide figures drifting into the bibliography.
**Choice:** existing graphics remain full-width because reducing these dense
panels to a single column would harm readability. **Inputs/outputs:** existing
assets/cells/captions, six numbered exhibits per build. **Gotchas:** call each
once per selected branch; labels and table/figure counters must stay consistent.
Wide floats can move to the next page. **R notes:** no plots or results rerun.

### Verification block

**Purpose:** extend the current manuscript reconciliation to shared exhibits.
**Logic:** read the master and shared definitions together for existing table,
asset and reference checks; require the wrapper to select/input the master;
reject manuscript commands in it; require one exhibit definition and two source
calls (one per branch), on opposite sides of the bibliography. The cap accepts
starred as well as ordinary floats. **Choice:** extend the existing checker
rather than introduce a new framework. **Inputs/outputs:** document sources and
authorized aggregate CSVs; explicit PASS or raised assertion. **Gotchas:** source
placement checks complement compiled/visual review; they do not prove every
scientific prose claim. **R notes:** Python standard library only, no dependency.

## Pending editorial proposal

The preview confirms the author's concern: the three tall NC graphics are
dense and repetitive, and the current main-exhibit balance favors application
detail over simulation explanation. These are pending editorial decisions, not
new analysis requirements. Proposed six-exhibit set:

1. Concise all-nine-design rules table.
2. Current grid MSE/coverage table by regime.
3. Updated colorful grid performance figure.
4. Compact, explicitly scoped grid bias–variance figure.
5. Observed NC incidence/region map; never the old synthetic map.
6. Compact NC annual performance figure by regime.

Keep allocation-risk findings discussed in the main text while placing detailed
tail/population plots in SI. Restore the original method-focused title if the
author approves. Future exhibit/title edits must reach both PDFs through the
shared sources and undergo numerical and displayed-size review again.
