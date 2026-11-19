# Post-submission planning prompt: CTJ, committee prep, SI curation and the prelim talk

Author: Andrew Walther. Written 2026-11-19, four days after the written prelim was submitted
(Sunday 2026-11-15) and 18 days before the oral exam (Monday 2026-12-07).

**How to use.** Start a new session in `SpatialCRT/projects/IncidenceDesign/` and add the
folders `bios-dissertation/` and `SpatialCRT/projects/SpillSpatialDepSim/`. Paste the fenced
block below. The session plans only; it writes one plan file and stops for decisions. Each
workstream is then run in its own session from that plan.

```text
You are planning, not executing. The goal of this session is one dated, sequenced plan,
docs/plans/2026-11-<dd>-post-prelim-plan.md (with a row in docs/plans/README.md), that takes
the IncidenceDesign project and its prelim integration from the submitted written prelim to
the oral exam on Monday 2026-12-07 and then to CTJ submission. Use the brainstorming skill
first. Do not edit manuscripts, code, slides or results in this session, and do not tag,
branch or change hooks; Workstream 0 proposes those for my approval.

## The hard constraint

The oral prelim is Monday 2026-12-07. The room is booked for 2 hours, the exam will likely
take about 90 minutes, and my talk is about 60 minutes (Parts I, II and III, about 20 minutes
each). Committee: Feng-Chang Lin (advisor, Project 2
corresponding author), Naim Rashid, Yusha Liu, Michael Hudgens, Tonya Van Deinse. Thanksgiving
falls on 2026-11-26. Rehearsals are fixed: practice with Dr. Lin on 2026-12-02 to 12-04, and
I present Part III at the Rashid Lab meeting on 2026-12-03. So Parts I and II must be
rehearsable by 2026-12-01. Estimate in real hours and days; do not pad the plan to fill the
window.

The order is mine and settled: (1) finish this iteration of the CTJ manuscript (revision check,
the structure and prose review, the word cut), (2) review, clean up and reorganize the SI,
then (3) the slides and script. The Chapter 2 review and the sensitivity study (Workstream B)
fit where the plan finds room, but the sensitivity results must exist before the backup slides
are built if they are to appear there. Test whether this order fits before 2026-12-01; if it
does not, say what slips or shrinks and leave the choice to me.

## Read first

- Instructions: AGENTS.md here, ../../AGENTS.md, bios-dissertation/AGENTS.md, DECISIONS.md
  here and in bios-dissertation. Settled choices are not reopened.
- Status: docs/TODO.md, ROADMAP.md (progress log), bios-dissertation/ROADMAP.md.
- CTJ: paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/ (FCL feedback PDF and
  register FCL-001 to FCL-060, proposed-edits.md, word-limit-trim-plan.md),
  docs/meetings/2026-11-14-ctj-feedback-meeting-notes.md,
  docs/plans/2026-11-14-chapter3-review-session-plan.md (dispositions),
  paper/ctj_manuscript/submission_checklist.md,
  manuscript_ctj/revision_markup/ (CTJ_Reading_Markup_V6.pdf is a frozen copy taken after
  Dr. Lin's feedback was applied; it is not marked up, and the TODO line asking for my markup
  of it is wrong. The text is not yet cut to 3,500 words).
- Committee prep: presentation/prelim-qa-prep.md (Q-01 to Q-20),
  bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/committee-proposals.md.
  There is no equivalent question bank for Chapter 2 yet.
- SI: paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex (65 pages, 15
  sections), paper/ctj_manuscript/notes/si-reorg-map.md and si-ch3-alignment-log.md (SI edits
  not yet carried to Chapter 3 or CTJ), paper/thesis_chapter/notes/trimmed-material.md.
- Talks: SpillSpatialDepSim/presentation/ (Part I deck, OUTLINE.md, reviews rounds 1 to 4),
  presentation/ here (SUDDEN talk of 2026-10-29, built for clinicians and epidemiologists,
  OUTLINE.md), bios-dissertation/slides/prelim/ (README.md assembly rules, OUTLINE.md not yet
  approved, test/ theme scaffold, part-3/ mirrored from multiomicsGEP).
- The submitted record: bios-dissertation/prelim/full-draft/prelim.pdf and the commit that
  built it. Committee questions and errata start from the submitted Chapters 2 and 3, since
  that is what the committee read. The decks may use the current state of each project
  where it has moved on (most likely Chapter 4).
- Part III question prep, for the rehearsal plan only:
  multiomicsGEP/presentation/prelim-slides/prelim-qa-prep.md.
- Two things that trip readers of older notes. Notation changed on 2026-11-15 (β on Z, τ on
  X, ψ for spillover); the 2026-11-14 meeting notes and the SUDDEN deck still say γ or use the
  old letters. Word counts in the repo use different rules: the checklist's 3,498 (texcount,
  math excluded, before round 5) and the trim plan's 3,838 (math as one word, after round 5)
  cannot be compared directly.

## Workstream 0: freeze the submitted record (do first; small)

The plan must start here, because later edits can silently change what the committee read.
- A frozen snapshot of the whole bios-dissertation/prelim/ folder as submitted, dated
  2026-11-15 (the uploaded PDF and every source and assembly file), somewhere no hook or build
  script writes to. Propose its location.
- Work then continues in bios-dissertation/thesis/, which becomes an updated, expanded and
  improved version of the prelim, not a rebuild. Committee feedback on the written prelim
  goes there, and Chapter 4 will change most. thesis/ now holds only chapter stubs and a
  README, so propose how the prelim copy replaces or absorbs them, and how the sync hooks
  (Chapter 3 from SpatialCRT, Chapter 4 from multiomicsGEP) are pointed at thesis/ instead of
  prelim/.
- Propose tags for the submitted state in bios-dissertation and SpatialCRT (and
  multiomicsGEP), after confirming which commits produced the uploaded PDF. I approve the
  tags; pushing them is a separate ask.
- The SpatialCRT post-commit hook re-renders Chapter 3 into bios-dissertation and commits
  there whenever Thesis_Chapter.qmd, its figures, tools/ or the bib change. Decide with me
  whether post-submission Chapter 3 edits should flow into the prelim draft now, wait on a
  branch, or pause the hook until after the exam. The same question applies to the
  ctj-to-bios hook and the Part III slide mirror.
- Start an errata log (docs/findings/) for errors anyone finds in the submitted Chapters 2
  and 3. These get a spoken correction or a slide, never a quiet fix to the submitted PDF.

## Workstream A: CTJ revision check and word budget

- Verify the round 5 revision item by item against the raw sources: the FCL PDF annotations,
  my meeting notes (decisions D-1 to D-5; items M-, S-, A-, I-, W-, Q-) and the sign-off
  dispositions. Produce a table: item, disposition, where it landed in CTJ_Manuscript.tex /
  CTJ_Exhibits.tex, whether it reached Chapter 3, status (done, partial, missing, deviates).
  Check both layouts were rebuilt after the last edit.
- Re-establish the count. The trim plan counts 3,838 body words (math as one word) against a
  3,500 limit; the checklist's 3,498 was before round 5. Look up the current Clinical Trials
  author guidelines (body limit, how they count, abstract limit; checklist says 425, current
  abstract about 404) and pick one counting rule, scripted so it can be rerun after every edit.
- Turn word-limit-trim-plan.md (T1 to T14) into a per-section budget: words now, target,
  candidate cuts, where each cut goes (SI section with a pointer, or trimmed-material.md).
  T1 conflicts with FCL-032; flag it for me. Note any cut that changes an SI pointer, since
  that couples Workstream A to Workstream C.
- Edits are proposed as one sign-off list in the Chapter 1 round 5 format (ID, source,
  location, current to proposed text) before anything is applied, then a second list for
  Chapter 3 carry-over. After the trim, a reading-flow pass, so the cut paper reads as one
  piece.
- Blockers that do not gate the trim but gate submission: IRB determination and consent
  wording (Simpson), permission to publish 58-cluster aggregates and maps, who is assessed in
  control colleges and what the education outcome is, the pilot description, coauthor
  approvals, Lin's open prep questions (docs/meetings/2026-10-14-advisor-meeting.md). Put
  them on a separate list with the person to ask and the earliest useful date.

## Workstream B: committee questions and new analyses

- Build a question bank per committee member from the submitted Chapters 2 and 3 (and the
  CTJ draft), extending Q-01 to Q-20. Think about each member's work. Hudgens will likely
  ask about estimands under interference and how the saturation designs relate to
  two-stage randomized saturation designs (Hudgens and Halloran 2008; Baird et al. 2018):
  what β estimates, whether the spillover model is a partial-interference assumption, and
  what a design-based analysis would give. Van Deinse will likely ask about the applied
  settings (probation in Project 1, the community-college education trial in Project 2) and
  implementation. Rashid, Liu and Lin: the plan should propose angles, not guess
  answers.
- Chapter 2 gets the same two-reviewer treatment Chapters 3 and 4 had in round 5 (see
  round5/ch3/ and round5/ch4/: committee proposals, math check, consistency check), run on
  the submitted Chapter 2, to produce feedback and a question bank like prelim-qa-prep.md.
  Project 1 is final and read-only: findings become spoken answers, slides or errata, not
  edits to the submitted chapter, and any new Project 1 computation needs my approval and
  lives in SpillSpatialDepSim/review/ or presentation/, never in the published code.
  Propose where the outputs go (round5/ch2/ or a new folder).
- Triage every question: spoken answer, backup slide, main slide, or new analysis. A new
  analysis beyond the sensitivity study below is justified only if its answer could change a
  stated conclusion or a committee member could not be answered without it. Say which, and
  give me the candidates as a list to choose from.
- The Incidence-Guided Saturation sensitivity study is mine to run; the plan sizes and scopes
  it rather than deciding whether it happens. The question it answers: we tested one version
  (four regions, 80/60/40/20%), and we do not know whether other numbers of regions or other
  saturation sets perform similarly relative to SRS. It should be a simple analysis. Note that prelim-qa-prep.md Q-16 currently
  answers it as future work, so that answer changes once results exist.
- Pre-specify the Incidence-Guided Saturation sensitivity study (Q-16) before anything runs,
  as a short spec for me to approve:
  - Factors: number of regions and the saturation set. The 10×10 grid splits evenly into 2
    halves, 4 quadrants of 5×5 or 25 blocks of 2×2; other counts (9, 16) need unequal
    regions, so say how they would be drawn. On the NC map, regions follow the existing
    at-least-2-colleges rule. Candidate saturation sets: flatter (0.6/0.4 or 0.65/0.55/0.45/0.35),
    current 0.8/0.6/0.4/0.2, steeper (0.9 to 0.1). Saturations of 0 or 1 leave a region with no
    within-region contrast; say whether that is allowed.
  - Hold the share treated at exactly 50% so designs stay comparable, and state the rounding
    rule that does it (0.65 of a 25-cluster quadrant is 16.25 clusters; a 2×2 block allows only
    quarters). Prefer equal-size partitions; justify any unequal one. Say how a saturation set
    maps onto the regions when the counts differ (for example four values over 25 blocks).
  - For each guided configuration, run the matching unguided Saturation counterpart (same
    regions and saturations, random rank), so the incidence-guidance effect is separated from
    the region-count and saturation effects. Name the contrasts in advance. Primary: each
    version's MSE relative to SRS, computed as in the existing SRS benchmark (code/15).
    Secondary: guided minus matched unguided MSE, paired by surface. With 10 surfaces per configuration and several
    factor levels, say how many contrasts are reported and how multiplicity is handled.
  - Implement the new designs as a stand-alone function in the new script with new seed-key
    tags. Do not add a design ID to get_designs() in 03 (that would change the hashed code);
    prelim-qa-prep.md Q-16 suggests "one more design ID", which this overrides.
  - Grid, NC or both; spillover regimes; queen primary; oracle estimator; t + df intervals;
    which metrics (MSE relative to SRS, bias, coverage, burden share treated).
  - Reuse the stored incidence surfaces and key-seeded draws; new script (next free number,
    29) and a new output folder; never edit the hashed 01 to 05 or production outputs.
  - Compute estimate from the timing of the existing runs, and a smoke and pilot stage first.
  - Where results go: SI and a backup slide, not the CTJ main text.
- Consider and cost, without committing to, other analyses the question bank suggests, for
  example covariate-constrained randomization on incidence (Q-12), misspecified exposure
  beyond the existing non-oracle and spillover-form fits (Q-05), or more surfaces per
  configuration (Q-13). For each, first list existing results that already answer it.

## Workstream C: SI curation and prose review

- Build two documents from one source: a full-length SI that keeps every method and result
  (the archive), and an abridged SI for submission. Decide with me whether the full-length
  version is a new file or is Chapter 3's Appendix B, which already holds the long form.
- Method: make a claim inventory. For each SI section and exhibit, record which CTJ
  sentence points to it, what claim it supports, and whether a reviewer would miss it.
  Sections that no main-text pointer reaches are candidates for full-length only. Propose
  a target length for the abridged SI and an order (methods before results, one numbered
  list, as in the prelim appendices).
- Resolve si-ch3-alignment-log.md (SI edits not in Chapter 3 or CTJ) during the curation.
- Prose review after curation, not before: human-writer pass, then a statistician reviewer
  and a clinician reviewer reading the abridged SI cold. Run verify_manuscripts.py with
  --pending-si on candidate copies, and code/tests/test_manuscript_exhibits.R.
- Fold in the Workstream B results once they exist.

## Workstream D: Parts I and II, the script, backups, assembly and rehearsal

- The script is the engine. For each part, write the spoken script first (one paragraph per
  slide with its time and running total), then make each slide support what is said: a claim
  as the title, one figure or table, few words. Nothing on a main slide that the narration
  does not use. Existing decks keep read-aloud scripts in speaker notes; keep that, and also
  produce a readable script document per part.
- All three parts are first drafts and each needs a full review and update before the exam.
  Part I was drafted for the prelim and will later be adapted for Dr. Lin's conference talk.
  Part II began as the SUDDEN talk. Part III is presented for practice at the Rashid Lab
  meeting on 2026-12-03 and is edited in multiomicsGEP only (its review is planned there;
  this plan schedules it and the assembly).
- Part I (SpillSpatialDepSim deck, four review rounds, corrected results): what changes for a
  biostatistics committee, figure sharpening, fit to 20 minutes.
- Part II: keep the SUDDEN deck as presented on 2026-10-29, unchanged, and build the prelim
  Part II deck as a new file from a copy of it. Re-aim it from clinicians to the committee and bring it up to the
  submitted Chapter 3: β on Z, τ on X, ψ for spillover; "Case Study" instead of
  "Application"; the canonical design names; t + df intervals; SRS in every design test.
  Check every number and figure against the submitted chapter.
- Backups: one point per slide, indexed by question ID from Workstream B, more detail
  allowed but still one figure or table each. Include the new sensitivity results if ready.
  Flag any slide that depends on an open external fact. One known case: the
  case-identification flowchart (412,514 deaths, then 111,665 aged 18 to 64, then 23,523).
  My understanding is that the algorithm is Gan et al.'s, Habib (2026) applied it to get the
  counts, and the counts were later corrected by reclassifying cardiac cases, which raised
  the total. So the slide should cite Gan et al. for the method and Habib (2026) for the
  counts. Check that the two upstream counts belong to the same corrected count series as
  23,523, and record the answer in presentation/data/case_identification_counts.csv.
- Assembly in bios-dissertation/slides/prelim/ under its README rules: shared title and
  agenda, Parts I to III with part footers, one thank-you slide, three backup sections, one
  theme. OUTLINE.md there needs my approval; propose what goes in the shared opening and
  the closing synthesis across the three projects. Part III is edited in multiomicsGEP only.
- Rehearsal: timed solo runs per part (recorded) before 2026-12-02; practice with Dr. Lin on
  2026-12-02 to 12-04; Part III at the Rashid Lab meeting on 2026-12-03; a full timed run
  and Q&A drills from the question banks before 2026-12-07. Leave time after the Lin
  sessions for their changes.
- After Part I is final: the standalone Project 1 talk is Part I with its own title and
  agenda slides. Plan it for after the exam.

## Workstream E: housekeeping and things I may have forgotten

Check TODO.md, both ROADMAPs and the memory index for open items and place each in the plan
or say why it waits: the read-only comprehensive review prompt (2026-10-26), the
reproducibility plan (on hold), the stale application report input
(downside_by_design.csv), the optional proportional-versus-on/off spillover sentence, the
known minor code issues, the Simpson questions (who is assessed in control colleges, who the
pilot trains and what it measures; no new information, so they stay open and must not block
the slides), Part III rehearsal and Q&A alongside Parts I and II, and the final defense.
List anything else you find that I did not mention.

The post-exam section of the plan back-schedules from these dates (from me, 2026-11-19):
- Final exam: Friday 2027-04-09, 10am to noon (primary); Wednesday 2027-04-07 (backup).
  The AIM requires 4 months between the prelim and the final exam, so 2027-04-07 is the
  earliest date.
- Abstract to internal.bios.unc.edu 5 weeks before the exam: 11:59pm 2027-03-04.
- Thesis to internal.bios.unc.edu 3 weeks before the exam: 11:59pm 2027-03-18.
- Graduate School dissertation (ETD) deadline: Wednesday 2027-04-14, five days after the
  primary exam date, so committee revisions must be small or done before the exam.
  The ETD abstract is limited to 350 words (bios-dissertation ROADMAP).
Include a step that updates the bios-dissertation ROADMAP's admin items and ETD risk flag with
these dates.

Two TODO items are probably superseded and should be confirmed and closed, not redone: my
review of the step 7b draft (2026-10-25) and my review of the concise chapter with both CTJ
layouts and the SI (2026-10-13). The round 5 review and the submitted prelim covered that
material. The CTJ structure and prose review plan
(docs/plans/2026-10-22-ctj-structure-and-prose-review-plan.md) is not superseded: it belongs
in Workstream A, before the word cut, and its SI parts in Workstream C.

## What the plan file must contain

1. A calendar from today to 2026-12-07 and a post-exam section, with hours per task.
2. For each workstream: steps, the raw sources each step starts from, the tools and checks,
   the output files, and a verification that does not rely on you or me (verifier scripts,
   rebuilt PDFs, a scripted word count, a timed rehearsal, a cold reviewer, my sign-off).
3. Dependencies between workstreams (CTJ cuts change SI pointers; analyses feed SI and
   backups; the question bank feeds backups and the analyses).
4. Which steps can run in parallel sessions and which must run one at a time because they
   touch the same files.
5. Every decision that is mine, gathered in one list at the end, each with your
   recommendation and what it blocks.

## Stops

Stop and wait for me after (1) the freeze decisions in Workstream 0, (2) the draft calendar
(including whether my order fits before 2026-12-01), and (3) the full plan. Before writing the plan file, get an independent review of
the draft plan from a fresh reviewer who sees this prompt and the plan but not your reasoning.
```
