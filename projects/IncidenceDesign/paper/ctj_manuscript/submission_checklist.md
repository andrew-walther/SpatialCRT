# Clinical Trials submission checklist — current full draft

Checked 2026-10-02 against [journal instructions](https://journals.sagepub.com/author-instructions/ctj).

- [x] Full original-research draft derived from the canonical thesis chapter.
- [x] Body 3,275 prose/heading words; 3,402 including end declarations
  (3,440 with mathematical-expression units as well);
  structured abstract 285 words; six keywords. Recounted 2026-10-12 with
  texcount on the body, declaration and abstract fragments, the same method
  that gave the earlier 3,215/3,381/283.
  The body count excludes the abstract, declarations, references and exhibits.
- [x] Six main exhibits: two tables and four figures; extra detail in supplement.
- [x] Standalone introduction and citation to verified preceding BMC publication.
- [x] Author/title page, affiliations, corresponding email and short running head.
- [x] Double-spaced single-column review text; end-positioned exhibits and separate legends.
- [x] Confirmed current author list, corresponding author, no funding/conflicts.
- [x] Data-availability statement preserves restricted sources; aggregate outputs linked.
- [x] Generative-assistance acknowledgement removed from the draft (author, 2026-10-12).
- [x] Current numerical/source verification and fresh rendered PDFs.
- [x] Separate shared-source two-column reading preview; submission layout retained.
- [x] Restore approved original title and revise six exhibits for two-column readability.
- [x] Author-approved order: Walther, Habib, Simpson, Lin; Lin corresponding.
- [ ] Author review of the revised full drafts and chosen exhibits.
- [ ] Ethics/data-use statement: the in-text placeholder was removed (author, 2026-10-12).
  Before submission, check whether the journal requires an ethics statement and, if so,
  add the IRB determination for the death-certificate data (Habib prepared the IRB
  documentation; Methods credits him).
- [ ] Corresponding-author telephone removed; email only. Check whether the submission
  system requires a telephone number.
- [ ] ORCIDs: Walther 0000-0002-1932-9653 is on the title page; collect the coauthors'
  ORCIDs in the submission system if required.
- [ ] Acknowledgements omitted for now. Check whether the journal requires an author
  contributions statement or a generative-AI disclosure; if so, confirm each author's
  role (Habib: data preparation and IRB documentation).
- [ ] Obtain every coauthor's approval of the final text.
- [ ] Confirm originality, exclusivity and permissions; disclose dissertation availability
  as applicable. No essential material from the preceding accepted study was copied.
- [ ] Complete the author-side cover letter, date/signature and any requested reviewer fields.
- [ ] Submit after author approval; no submission or public release was performed.
  (Repository commits are pushed to `origin/main`, verified 2026-10-11.)

This is a completed simulation-based design study using historical geographic
inputs. It reports no implemented education trial or clinical efficacy results;
trial registration/CONSORT results must not be fabricated. Interval undercoverage
and finite allocation-tail uncertainty remain scientific limitations in the text.

Files: `CTJ_Manuscript.{tex,pdf}`, `Supplementary_Information.{tex,pdf}`,
`Figure_Legends.md`, the four referenced vector PDFs and the shared bibliography.
Rebuild from this directory with TinyTeX on PATH, `TEXINPUTS` and `BSTINPUTS`
including `../SAGE_Journal_Template`; run `latexmk -pdf` on each TeX source.

`CTJ_Reading.tex`/PDF is an author preview, not a submission file or publisher
proof. It selects a layout branch in `CTJ_Manuscript.tex`; both layouts load
`CTJ_Exhibits.tex`. Rebuild both after changing the master or exhibit definitions.
