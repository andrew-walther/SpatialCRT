# Clinical Trials submission checklist — current full draft

Checked 2026-10-02 against [journal instructions](https://journals.sagepub.com/author-instructions/ctj).

- [x] Full original-research draft derived from the canonical thesis chapter.
- [x] Body 2,217 words; structured abstract 283 words; six keywords.
- [x] Six main exhibits: two tables and four figures; extra detail in supplement.
- [x] Standalone introduction and citation to verified preceding BMC publication.
- [x] Author/title page, affiliations, corresponding email and short running head.
- [x] Double-spaced single-column review text; end-positioned exhibits and separate legends.
- [x] Confirmed current author list, corresponding author, no funding/conflicts.
- [x] Data-availability statement preserves restricted sources; aggregate outputs linked.
- [x] Generative-assistance disclosure draft included for author review.
- [x] Current numerical/source verification and fresh rendered PDFs.
- [ ] Replace author-approved ethics/data-use placeholder with institutional determination,
  appropriate protocol/approval/waiver details, consent determination and source/aggregate
  publication permissions. No approval or exemption is presently asserted.
- [ ] Supply Lin's telephone and author ORCIDs for submission.
- [ ] Confirm contribution wording and obtain every coauthor's approval of the final text.
- [ ] Review assistance disclosure for accuracy/completeness, including any other tools used.
- [ ] Confirm originality, exclusivity and permissions; disclose dissertation availability
  as applicable. No essential material from the preceding accepted study was copied.
- [ ] Complete the author-side cover letter, date/signature and any requested reviewer fields.
- [ ] Submit after author approval; no submission, public release or Git push was performed.

This is a completed simulation-based design study using historical geographic
inputs. It reports no implemented education trial or clinical efficacy results;
trial registration/CONSORT results must not be fabricated. Interval undercoverage
and finite allocation-tail uncertainty remain scientific limitations in the text.

Files: `CTJ_Manuscript.{tex,pdf}`, `Supplementary_Information.{tex,pdf}`,
`Figure_Legends.md`, the four referenced vector PDFs and the shared bibliography.
Rebuild from this directory with TinyTeX on PATH, `TEXINPUTS` and `BSTINPUTS`
including `../SAGE_Journal_Template`; run `latexmk -pdf` on each TeX source.
