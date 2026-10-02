# Full prelim and CTJ reference-flow checkpoint

Author approved implementation October 2. The initial full-length working prelim was
`bios-dissertation/prelim/prelim.pdf`: 177 pages using the existing bios-prelim
class and identity file. Bodies are 24 pages for Chapter 1, 42 for Chapter 2,
42 for Chapter 3 and 18 for Chapter 4. Appendices A/B are 4/25 pages; references
12 and front matter 10. The author retained Chapter 3's full appendix after
reviewing Amber Young's 143-page A–N appendix as context; no appendix cuts made.

The assembly reads current working drafts, including uncommitted author-feedback
literature-review and master-bibliography changes. Those files were not edited
or staged. All chapter prose/math/captions are preserved by mechanical assembly;
the Project 1 analysis and accepted manuscript remain read-only. The new
299-word front abstract is an initial draft and explicitly limits current
Chapter 4 supervision claims. Its longer standalone abstract was not rewritten.

Build from bios-dissertation's root:

```bash
PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$PATH" python3 prelim/tools/build_prelim.py
python3 prelim/tools/test_prelim_assembly.py
```

Five tests pass, all 177 PDF pages reviewed as contact sheets, and important
boundaries/exhibits reviewed at larger size. Four chapters, two existing
appendices, 122 bibliography keys and 37 figure assets checked. There are no
unresolved `??` placeholders. Quarto's 25 raw-table notices are documented;
these labels deliberately use LaTeX references, which resolve. Source hashes
and output PDF hash are in `prelim/assembly_manifest.json`; complete build
teach-back and page inventory are in `prelim/notes/full-assembly-2026-10-02.md`.

## CTJ layout purpose and logic

The user's desired reading layout starts references after declarations and
uses both columns on the final page. In `CTJ_Manuscript.tex`, `flushend` is
loaded only for the reading flag; both forced bibliography page breaks in
that branch were removed. The already-installed package balances the final
page automatically, avoiding a manually placed balancing command that can
fall in the wrong column. The `lastpage` package is review-only because Sage
already defines the reading page label and two competing labels can disagree
when the bibliography's final page is balanced.

Inputs: shared CTJ master, reading wrapper, existing exhibits and bibliography.
Outputs: updated eight-page reading PDF and rebuilt 21-page submission PDF.
No scientific prose, author ordering, numbers, exhibits or word count changed.
References begin on page 7 and continue across both columns on page 8, leaving
white space below both columns. Column bottoms need not align to the exact
same line because reference entries have hanging indentation and line-break
constraints. The reading render has no warnings, overfull boxes or unresolved
references. Submission extracted text matches its previous version exactly;
all numerical/exhibit manuscript checks pass. This is an author preview, not
publisher proof formatting. No R dependencies or computation were changed.

## Remaining review

- Review full prelim and its initial front abstract; pending literature-review
  revisions remain working content, not implied author approval.
- Chapter 4 source timeline mentions Summer 2027 defense while bios-dissertation
  records an April 2027 target. Flag for author/source-owner review, without
  rewriting Project 3 here. No scientific Project 3 claims were strengthened.
- Full prelim assembly must be rerun after chapter syncs or literature-review
  updates; existing Project 2/3 sync hooks still generate their individual
  chapters. The hooks were not expanded to auto-render the full document.
- Any condensed committee version is a separate later deliverable. No final
  department-submission/ETD compliance claim is made. Author review, CTJ
  declarations and submission approval remain pending. No push performed.

## Author-requested tables and figures lists

The updated prelim is 187 pages. Its master calls `\listoftables` and
`\listoffigures` after the contents and before main matter, overriding the
class default for prelim mode without editing the class or chapter sources.
The List of Tables has 29 entries on printed pages xi–xv; the List of Figures
has 35 entries on xvi–xx. Both include appendix exhibits and appear in the
contents. Existing full captions remain the single source for each entry.

All ten new pages were visually reviewed, and all 64 entry page references
match the body captions. Five assembly tests and the PDF manifest hash check
pass. All 167 chapter/appendix/reference pages have identical extracted text
and unchanged printed numbering; only their PDF positions move by ten pages.
Pending literature-review/bib edits remain unstaged. No scientific changes
or changes to CTJ were needed. No push.
