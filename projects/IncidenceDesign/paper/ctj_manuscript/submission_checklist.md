# Clinical Trials submission checklist — current full draft

Checked 2026-10-02 against [journal instructions](https://journals.sagepub.com/author-instructions/ctj).

- [x] Full original-research draft derived from the canonical thesis chapter.
- [x] Body 3,362 prose/heading words (after the 2026-10-14 Monte Carlo definition); limit 3,500 applies to the body (author, 2026-10-14;
  end declarations, abstract, references and exhibits excluded). 3,628 with declarations
  (266 declaration words; recounted 2026-10-14, texcount text + headers, math excluded).
  Placeholder text (ethics statement, coauthor-approval sentence) is highlighted yellow.
  Structured abstract 394 words (limit 425); six keywords. Recounted 2026-10-14
  after the heading/abstract revision (informative subheadings, Discussion
  subsections, plainer abstract), using texcount on the body, declaration and
  abstract fragments.
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
- [x] Comprehensive review completed 2026-10-12 (see docs/reviews/2026-10-12-ctj-comprehensive-review-findings.md).
- [ ] Author review of the revised full drafts and chosen exhibits.
- [ ] Ethics/data-use statement: CTJ requires an IRB statement in Methods for human data.
  **Protocol filled in (2026-10-21):** UNC IRB 14-2036, "Sudden Unexpected Death in North
  Carolina" (PI Ross J. Simpson, Jr., Department of Medicine). The determination type and
  consent wording remain a highlighted placeholder pending Simpson (Sage also asks for the
  statement on the title page); confirm permission to publish 58-cluster aggregates.
  Body word count not recounted after the 2026-10-21 additions (Gan citation sentence,
  IRB protocol sentence; roughly +40 words).
- [x] Corresponding-author telephone (919) 966-9431 on the title page and cover letter
  (from Lin's UNC Gillings profile; the journal requires it).
- [ ] ORCIDs: Walther 0000-0002-1932-9653 is on the title page; collect the coauthors'
  ORCIDs in the submission system if required.
- [x] Author contributions statement added (author-supplied roles, 2026-10-12).
- [x] Generative-AI disclosure: Sage policy requires one (methods or acknowledgements)
  for AI-assisted code, code checking and visualization. Draft added as Acknowledgements
  (2026-10-14); author to confirm the tool list and wording before submission.
- [ ] Obtain every coauthor's approval of the final text, then delete
  "[Placeholder: pending coauthor approval.]" after the approval sentence in Author contributions.
- [ ] Confirm originality, exclusivity and permissions; disclose dissertation availability
  as applicable. No essential material from the preceding accepted study was copied.
- [ ] Complete the author-side cover letter, date/signature and any requested reviewer fields.
- [ ] Submit after author approval; no submission or public release was performed.
  (Repository commits are pushed to `origin/main`, verified 2026-10-11.)

This is a completed simulation-based design study using historical geographic
inputs. It reports no implemented education trial or clinical efficacy results;
trial registration/CONSORT results must not be fabricated. Interval undercoverage
and finite allocation-tail uncertainty remain scientific limitations in the text.

Files: `manuscript_ctj/CTJ_Manuscript.{tex,pdf}`, `supplementary_info/Supplementary_Information.{tex,pdf}`,
`figures/Figure_Legends.md`, the four referenced vector PDFs in `figures/compact/`, the
cover letter in `cover_letter/` and the shared bibliography. Rebuild with TinyTeX on PATH:
run `latexmk` on each TeX source inside its folder. The local `latexmkrc` and
template symlinks need no extra path setup.

`CTJ_Reading.tex`/PDF is an author preview, not a submission file or publisher
proof. It selects a layout branch in `CTJ_Manuscript.tex`; both layouts load
`CTJ_Exhibits.tex`. Rebuild both after changing the master or exhibit definitions.
