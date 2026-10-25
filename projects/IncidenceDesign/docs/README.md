# docs/

Project documents, sorted by what kind of document they are. Each folder has a
README that lists its documents with a one-line purpose and a status.

| Folder | What goes there |
|---|---|
| [`plans/`](plans/README.md) | Plans and specs, written before the work and approved by the author |
| [`findings/`](findings/README.md) | Results, completion checkpoints, reports and verification logs from work that ran |
| [`reviews/`](reviews/README.md) | Reports from independent reviewers |
| [`prompts/`](prompts/README.md) | Paste-in prompts, review briefs and session handoffs |
| [`meetings/`](meetings/README.md) | Meeting records |
| [`TODO.md`](TODO.md) | Open items (living document) |

## Naming

`YYYY-MM-DD-topic-type.md`, all hyphens. The date is the day the document was
written. The date-first form keeps each folder in time order and sorts correctly
across years.

When you add a document, give it a row in its folder's README.

## Status words

- **active**: still in use, or still guides open work.
- **complete**: the work is done. The document is the record.
- **superseded**: a later document replaces it (named in the row).
- **not started**: approved or drafted, waiting on something.

## Renamed on 2026-10-25

Before 2026-10-25 every document sat in `docs/plans/`, with the date at the end of
the name. A few files still cite the old paths and were deliberately left alone:

- `code/04_estimation.R` and `code/05_run_simulation.R`. Their hashes are in the
  checkpoint manifest, so editing even a comment blocks loading the checkpoints.
- `application/code/methods_ext_core.R`, `application/code/run_methods_ext_queue_v2.sh`
  and the generated `application/report/real_sud_companion.html`.
- `paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex` (a comment on line 2).
  Committing it would trigger the CTJ sync into bios-dissertation.
- Frozen records in `paper/archive/`, `results/archive/` and the repository-level `archive/`.
- Files in the bios-dissertation repository.

Use the table to find the new name for an old path.

| Old path | New path |
|---|---|
| `docs/plans/TODO.md` | `docs/TODO.md` |
| `docs/superpowers/plans/2026-04-07-tau-se-power-plan.md` | `docs/plans/2026-04-07-tau-se-power-plan.md` |
| `docs/superpowers/specs/2026-04-07-tau-sensitivity-design.md` | `docs/plans/2026-04-07-tau-sensitivity-design-spec.md` |
| `docs/plans/manuscript-unification-plan.md` | `docs/plans/2026-09-18-manuscript-unification-plan.md` |
| `docs/plans/simulation-revision-plan.md` | `docs/plans/2026-09-24-simulation-revision-plan.md` |
| `docs/plans/simulation-revision-spec.md` | `docs/plans/2026-09-24-simulation-revision-spec.md` |
| `docs/plans/allocation_risk_plan_2026-09-27.md` | `docs/plans/2026-09-27-allocation-risk-plan.md` |
| `docs/plans/nc_sud_application_plan_2026-10-02.md` | `docs/plans/2026-10-02-nc-sud-application-plan.md` |
| `docs/plans/nc_sud_implementation_plan_2026-10-02.md` | `docs/plans/2026-10-02-nc-sud-implementation-plan.md` |
| `docs/plans/chapter3_length_review_plan_2026-10-12.md` | `docs/plans/2026-10-12-chapter3-length-review-plan.md` |
| `docs/plans/methods_extensions_plan_2026-10-20.md` | `docs/plans/2026-10-20-methods-extensions-plan.md` |
| `docs/plans/overnight_plan_2026-10-21.md` | `docs/plans/2026-10-21-overnight-plan.md` |
| `docs/plans/ctj_structure_and_prose_review_plan.md` | `docs/plans/2026-10-22-ctj-structure-and-prose-review-plan.md` |
| `docs/plans/pilot_2026-09-24/` | `docs/findings/2026-09-24-pilot-m1/` |
| `docs/plans/allocation_risk_findings_2026-09-27.md` | `docs/findings/2026-09-27-allocation-risk-findings.md` |
| `docs/plans/ch2_check_findings_2026-09-27.md` | `docs/findings/2026-09-27-ch2-check-findings.md` |
| `docs/plans/nc_sud_application_findings_2026-10-02.md` | `docs/findings/2026-10-02-nc-sud-application-findings.md` |
| `docs/plans/nc_sud_implementation_findings_2026-10-02.md` | `docs/findings/2026-10-02-nc-sud-implementation-findings.md` |
| `docs/plans/manuscript_completion_findings_2026-10-02.md` | `docs/findings/2026-10-02-manuscript-completion-findings.md` |
| `docs/plans/ctj_exhibit_revision_findings_2026-10-02.md` | `docs/findings/2026-10-02-ctj-exhibit-revision-findings.md` |
| `docs/plans/ctj_layout_review_2026-10-02.md` | `docs/findings/2026-10-02-ctj-layout-review.md` |
| `docs/plans/prelim_assembly_ctj_reference_review_2026-10-02.md` | `docs/findings/2026-10-02-prelim-assembly-ctj-reference-review.md` |
| `docs/plans/manuscript_verification_2026-10-02.txt` | `docs/findings/2026-10-02-manuscript-verification.txt` |
| `docs/plans/manuscript_exhibit_verification_2026-10-02.txt` | `docs/findings/2026-10-02-manuscript-exhibit-verification.txt` |
| `docs/plans/methods_extensions_findings_2026-10-20.md` | `docs/findings/2026-10-20-methods-extensions-findings.md` |
| `docs/plans/methods_extensions_report_2026-10-21.md` | `docs/findings/2026-10-21-methods-extensions-report.md` |
| `docs/plans/step0-review-findings-2026-09-23.md` | `docs/reviews/2026-09-23-step0-review-findings.md` |
| `docs/plans/nc_sud_independent_review_2026-10-02.md` | `docs/reviews/2026-10-02-nc-sud-independent-review.md` |
| `docs/plans/ctj_comprehensive_review_findings_2026-10-12.md` | `docs/reviews/2026-10-12-ctj-comprehensive-review-findings.md` |
| `docs/plans/step0-chapter-review-prompt.md` | `docs/prompts/2026-09-23-step0-chapter-review-prompt.md` |
| `docs/plans/real-data-aggregation-session-prompt.md` | `docs/prompts/2026-09-25-real-data-aggregation-session-prompt.md` |
| `docs/plans/handoff_2026-09-27_srs_ch2_audit.md` | `docs/prompts/2026-09-27-handoff-srs-ch2-audit.md` |
| `docs/plans/continuation_prompt_2026-10-02.md` | `docs/prompts/2026-10-02-continuation-prompt.md` |
| `docs/plans/claude_comprehensive_review_prompt_2026-10-02.md` | `docs/prompts/2026-10-02-claude-comprehensive-review-prompt.md` |
| `docs/plans/claude_comprehensive_review_prompt_2026-10-12.md` | `docs/prompts/2026-10-12-claude-comprehensive-review-prompt.md` |
| `docs/plans/talk_2026-10-29_kickoff_prompt.md` | `docs/prompts/2026-10-25-talk-kickoff-prompt.md` |
| `docs/plans/advisor_meeting_2026-10-14.md` | `docs/meetings/2026-10-14-advisor-meeting.md` |
