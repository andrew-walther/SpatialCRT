# Clinical Trials submission checklist — current full draft

Checked 2026-10-02 against [journal instructions](https://journals.sagepub.com/author-instructions/ctj).

- [x] Full original-research draft derived from the canonical thesis chapter.
- [x] Body 3,280 prose/heading words; 3,425 including end declarations
  (about 3,480 with mathematical-expression units); limit 3,500.
  Structured abstract 317 words (limit 425); six keywords. Recounted 2026-10-12
  after the comprehensive-review edits (abbreviations defined in the body, plain
  language), using texcount on the body, declaration and abstract fragments.
  The body count excludes the abstract, declarations, references and exhibits.
- [x] Six main exhibits: two tables and four figures; extra detail in supplement. Main text
  shows the six detailed designs and SRS; the two related variants are in the supplement
  (author 2026-10-12).
- [x] Standalone introduction and citation to verified preceding BMC publication.
- [x] Author/title page, affiliations, corresponding email and short running head.
- [x] Double-spaced single-column review text; end-positioned exhibits and separate legends.
- [x] Confirmed current author list, corresponding author, no funding/conflicts.
- [x] Data-availability statement: simulation data and code available on request; NC death-certificate data restricted (no repository link, author decision 2026-10-12).
- [x] Generative-assistance acknowledgement removed from the draft (author, 2026-10-12).
- [x] Current numerical/source verification and fresh rendered PDFs.
- [x] Separate shared-source two-column reading preview; submission layout retained.
- [x] Restore approved original title and revise six exhibits for two-column readability.
- [x] Author-approved order: Walther, Habib, Simpson, Lin; Lin corresponding.
- [x] Comprehensive review completed 2026-10-12 (see docs/plans/ctj_comprehensive_review_findings_2026-10-12.md).
- [ ] Author review of the revised full drafts and chosen exhibits.
- [ ] Ethics/data-use statement: CTJ requires an IRB statement in Methods for human data.
  Obtain the IRB determination for the death-certificate data (Habib/Simpson), including
  permission to publish 58-cluster aggregates, and add it to Methods.
- [x] Corresponding-author telephone (919) 966-9431 on the title page and cover letter
  (from Lin's UNC Gillings profile; the journal requires it).
- [ ] ORCIDs: Walther 0000-0002-1932-9653 is on the title page; collect the coauthors'
  ORCIDs in the submission system if required.
- [x] Author contributions statement added (author-supplied roles, 2026-10-12).
- [ ] Generative-AI disclosure: Sage policy requires one (methods or acknowledgements)
  for AI-assisted code, code checking and visualization. Author will write it before
  submission.
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
