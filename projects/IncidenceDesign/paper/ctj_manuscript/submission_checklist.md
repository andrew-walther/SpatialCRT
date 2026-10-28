# Clinical Trials submission checklist — current full draft

Checked 2026-10-02 against [journal instructions](https://journals.sagepub.com/author-instructions/ctj). The 2026-10-27 proofread could not load the page (HTTP 403); re-check it before submission.

- [x] Full original-research draft derived from the canonical thesis chapter.
- [x] Body 3,498 prose/heading words (recounted 2026-10-27 after the V5 revision, Introduction through
  Discussion, texcount text + headers, math excluded); limit 3,500 applies to the body
  (author, 2026-10-14; end declarations, abstract, references and exhibits excluded).
  Declarations not recounted since "Ethical considerations" was added (2026-10-27).
  Placeholder text (ethics statement, coauthor-approval sentence) is highlighted yellow.
  Structured abstract 404 words (recounted 2026-10-27; was 411 on 2026-10-26; limit 425); six keywords. Earlier count 394 on 2026-10-14
  after the heading/abstract revision (informative subheadings, Discussion
  subsections, plainer abstract), using texcount on the body, declaration and
  abstract fragments.
  The body count excludes the abstract, declarations, references and exhibits.
- [x] Six main exhibits: two tables and four figures; extra detail in supplement. Main text
  shows the six detailed designs and SRS; the two related variants are in the supplement
  (author 2026-10-12). Table 1 is a hybrid table with one example map per design
  (author, 2026-10-25); it stays Table 1 and the count stays six.
- [x] Standalone introduction and citation to verified preceding BMC publication.
- [x] Author/title page, affiliations, corresponding email and short running head.
- [x] Double-spaced single-column review text; end-positioned exhibits and separate legends.
- [x] Confirmed current author list, corresponding author, no funding/conflicts.
- [x] Data-availability statement: simulation data and code available on request; NC death-certificate data restricted (no repository link, author decision 2026-10-12).
- [x] Generative-assistance acknowledgement removed from the draft (author, 2026-10-12); superseded by the Acknowledgements disclosure (2026-10-14; tools Claude Code and Codex, confirmed 2026-10-27).
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
  2026-10-27: a short IRB sentence stays in the application Methods; the full statement and
  placeholder moved to an "Ethical considerations" declaration.
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
`figures/Figure_Legends.md`, the four referenced vector PDFs in `figures/compact/` and the eight Table 1 thumbnails in `figures/compact/design_thumbs/`, the
cover letter in `cover_letter/` and the shared bibliography. Rebuild with TinyTeX on PATH:
run `latexmk` on each TeX source inside its folder. The local `latexmkrc` and
template symlinks need no extra path setup.

`CTJ_Reading.tex`/PDF is an author preview, not a submission file or publisher
proof. It selects a layout branch in `CTJ_Manuscript.tex`; both layouts load
`CTJ_Exhibits.tex`. Rebuild both after changing the master or exhibit definitions.

## Before submission (from the 2026-10-27 final proofread, `docs/reviews/2026-10-27-ctj-v5-final-proofread.md`)

- [ ] Remove draft labels from the submission file: the title-page line "Clinical Trials
  submission draft. Final submission approval remains pending." and the running header
  (`\volumenumber`/`\issuenumber`).
- [ ] Ethics statement on the title page once Simpson supplies the determination wording.
- [ ] Confirm with the editorial office whether Table 1 with embedded example maps counts as a
  table (six exhibits) or a figure (seven).
- [ ] Reference style: convert title-case article titles (refs 3, 11, 12, 13, 17) to sentence
  case in both the shared and the bios-dissertation master bib (title is a sync-checked
  field); edition out of the Hayes title; check the BMC paper's exact title punctuation; drop
  `url` where a DOI exists.
- [ ] Figure 1 panel A prints "2018-2021" with a hyphen (R's pdf device cannot print an en
  dash); acceptable, or switch that export to a Unicode-capable device.
- [ ] Body is at 3,498 of 3,500 words: any addition needs an equal cut.
- [ ] Coauthor ORCIDs; funding confirmation with Habib and Simpson (the Funding comment).

