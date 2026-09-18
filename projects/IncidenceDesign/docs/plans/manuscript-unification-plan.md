# Plan: Write Project 2 long-form first; derive CTJ + supplementary
# material from it

**Written in a bios-dissertation session on 2026-09-18** (that's where
this planning conversation happened) and copied here. If this file and
`bios-dissertation/prelim/project-proposals/notes/project2-unification-plan-DRAFT-for-SpatialCRT.md`
ever diverge, this copy — the one actually next to the code and the
manuscripts it plans for — is authoritative; update or delete the other.

**Correction made during drafting, same session:** an earlier version of this plan
had the direction backward — CTJ + its Supplementary Information as the
primary source, with the dissertation chapter derived from it. That
matched Project 1's pattern, but Project 1 is the repo's own *named
exception* (its chapter is copy/paste from an already-accepted
manuscript, precisely because it couldn't be rewritten). `ROADMAP.md`'s
"Manuscript & Chapter Writing Strategy" already settles this for Project
2, and had before this session touched it: **write one long-form version,
unrestricted in length; the 3,500-word CTJ submission is a condensed
derivative, with material cut for length preserved in a
supplementary-material file rather than discarded.** This version of the
plan follows that.

## The four deliverables and their hierarchy

Four end products, not two documents each doing double duty:

1. **Dissertation chapter** (Chapter 3 body) — long-form, unrestricted.
2. **Dissertation chapter appendix** — exhaustive detail the chapter body
   doesn't need to carry to make its argument.
3. **CTJ manuscript** — ≤3,500 words / 6 exhibits, condensed from (1).
4. **CTJ supplementary material** — what gets cut going from (1)+(2) down
   to (3), formatted for journal submission, not just a scratch file.

The hierarchy is one direction only: **(1) and (2) are built and settled
first; (3) and (4) are derived from them afterward.** Both (1)/(2) and
(3)/(4) involve the same kind of decision — what does a reader need to
follow the argument, versus what would they only look up — applied twice,
independently, at two different length budgets. Getting (1) vs. (2) right
is not optional groundwork for getting (3) vs. (4) right later; a chapter
that dumps undigested detail into its body because there's no word cap to
force discipline will produce a worse chapter *and* a harder compression
step, since step 5 would then be compressing content that should have
been appendix material from the start rather than genuinely cutting for
length.

## Build section by section, not as one pass — components have
different fates in different documents

Not every piece of a manuscript belongs in all four deliverables. The
clearest existing precedent is Project 1's own port (see
`prelim/project-proposals/notes/chapter2-integration-plan.md`): the
journal abstract and the full `Declarations` block (ethics, consent,
competing interests, author contributions, funding, data availability)
came **out of the dissertation chapter entirely**, preserved verbatim in
a side file (`manuscript-declarations.md`) rather than deleted. Amber
Young's own dissertation, which that decision followed, keeps exactly two
of those as unnumbered sections after the chapter's Discussion — Funding
and Data Availability — and drops the rest.

Project 2 needs the same per-component decision, made explicit before
drafting rather than discovered by accident later. Working table, to be
filled in and corrected as each component is actually drafted:

| Component | Chapter body | Chapter appendix | CTJ manuscript | CTJ supplementary |
|---|---|---|---|---|
| Standalone abstract | **No** — the dissertation's own front matter carries its abstract; a per-chapter journal-style abstract is redundant with it (Project 1 precedent) | — | **Yes** — required by the journal | — |
| Introduction / Methods / Results / Discussion | Yes, in full | Overflow only (per the body-vs-appendix test in step 3) | Yes, condensed | Overflow from condensing |
| Full parameter-grid tables, robustness breakdowns | No | Yes | No | Yes |
| Application section | Yes, reframed as a proposal (step 2) | — | Placeholder, until real data | — |
| Ethics / consent / competing interests / author contributions | No | No | Yes (journal requires them) | No |
| Funding / Data Availability Statement | Unnumbered sections after Discussion, per Amber's precedent | — | Yes | — |
| Acknowledgments | Case-by-case | — | Yes | — |

Build in the order the steps below already imply — chapter body, then
chapter appendix, then the CTJ condensation — but **decide each
component's fate at the point it's drafted**, not retroactively once
everything exists. A component whose fate is "chapter only" doesn't need
journal-declaration boilerplate written into it in the first place; a
component whose fate is "CTJ only" (the abstract, the declarations)
shouldn't be drafted as part of the chapter and then deleted later.

## Context

Confirmed by direct investigation (2026-09-18): `paper/` in IncidenceDesign
has had zero commits since 2026-07-03. CTJ has not been submitted, blocked
on the real SUDDEN dataset (still pending IRB/access) and on a
"consolidated review pass across all three documents" that was flagged
open and never done. `Dissertation_Chapter.qmd` (7,305 words, renders to
28 pp. in its own format) is ~90% a restatement of content already in
`CTJ_Manuscript.tex` (2,458 words body) + `Supplementary_Information.tex`
(14 pp., complete) — the two genuinely new items anywhere in the long
version are a coverage-by-ρ table for Checkerboard and one sentence
describing the pilot intervention itself.

**Decision:** `Dissertation_Chapter.qmd` becomes the primary, long-form
source — expanded and completed, not just reformatted. `CTJ_Manuscript.tex`
and its Supplementary Information become **derived** from it, condensed by
hand to fit the journal's 3,500-word / 6-exhibit limit. The existing CTJ
prose isn't wasted — it's a strong first draft of already-compressed
phrasing for several claims and is worth reusing where the compression it
already found is good — but it stops being the source of truth once the
long-form chapter is complete.

**Sizing, calibrated against real renders, not guessed:**
Project 1's chapter rendered at 42 pp. body + 4 pp. appendix (13,902 words
of source prose, ≈302 words/page in `bios-prelim.cls`). The current
`Dissertation_Chapter.qmd` is already close in scale (7,305 words, ≈24 pp.
at that ratio) and is the better starting point than CTJ's 2,458 words —
expanding a document that's already ~90% of the way to a full argument is
far less work than expanding CTJ's compressed 8-page version twice (once
to a full chapter, then trimming a *different* subset back down for the
journal). Project 2's simulation study is also narrower in scope than
Project 1's (one design comparison vs. two spillover mechanisms across a
larger design space), so landing shorter than Project 1's 42 pp. is
expected, not a shortfall to correct.

## Steps

**0. Give `Dissertation_Chapter.qmd` the same review pass the literature
review sections got — it hasn't had one yet.**
This document was written independently, not reviewed against a citation
audit trail or checked for voice against Project 1's accepted manuscript
the way the lit review was. Before expanding it, verify what's already
there: every citation resolves and supports its claim, no mannered prose,
voice consistent with the rest of this dissertation.
→ *Verify:* a fresh reviewer pass (subagent or otherwise) returns no
unresolved findings on the existing content before new content is added
on top of it.

**1. Fold in the two genuinely new items already identified**, plus
whatever else the review in step 0 or a closer read surfaces: the
ρ-stratified Checkerboard coverage table, and the pilot-intervention
description sentence. Both are small, confirmed additions — this step is
mechanical, not open-ended drafting.
→ *Verify:* both appear in the chapter; neither duplicates content
already present elsewhere in it.

**2. Reframe the application section as a forward-looking proposal.**
The scaffolding already exists almost verbatim in the current text:
dataset description (SUDDEN algorithm, 2018–2021, NC death certificates,
IRB status), planned covariates, geography (58 CC service areas, Lenoir
County pilot site + Kinston/Jones satellites), and the hoped-for insight
(process measures, not incidence rate, on a pilot timescale). This is a
**deletion-and-reframing edit**, not new drafting: cut the synthetic
MSE/coverage numeric run and its figure/table; keep everything else; shift
the surrounding prose from "placeholder result shown, to be superseded"
to "analysis plan stated, not yet run." Mention that the R ingestion
pipeline already exists and is waiting on data access, so a reader knows
this isn't stalled for lack of a plan. This applies to the long-form
chapter now; CTJ's own application section stays as-is (still a
placeholder) until the real dataset is available and CTJ is actually
being prepared for submission — the two documents don't need to be
resolved on this point at the same time.
→ *Verify:* the section reads as an honest in-progress item with no
fabricated or placeholder numbers.

**3. Decide body vs. appendix for the chapter's own supplementary
content, using a real test, not a page-ceiling-driven default.**
A generous word budget is a reason the chapter body *can* carry more
detail than CTJ's SI does; it is not a reason it *should* carry
everything the SI contains. The test: **does a reader need this content
to follow the argument, or would they only consult it?**
- **Body:** the headline results and the comparisons that support the
  chapter's central claims — the MSE/coverage ranking, the omnibus
  Friedman result, the bias-variance decomposition explaining *why*
  Checkerboard fails, the reframed application section from step 2.
- **Appendix:** exhaustive breakdowns a reader looks up rather than reads
  through — the full per-incidence-mode / per-parameter (ρ, γ, spillover
  regime) / τ-level grids, the 8-design comparison (the chapter's main
  argument is scoped to 6 by design), reproducibility detail (seeding
  protocol, exact metric formulas). Much of this can come from
  `Supplementary_Information.tex` directly rather than being rewritten —
  it's already built for this role.
- **Judgment call, decide once and record it:** the two extra application
  maps. Body if a reader needs them to picture the geography the
  application section describes; appendix if they're one of several
  equivalent views.
For each candidate, write down which side of the test it landed on and
why — one line each — so the reasoning survives, not just the outcome.
→ *Verify:* re-render after each batch; running body and appendix page
counts recorded; each placement decision has its one-line reasoning
somewhere (a source comment is sufficient).

**4. Render the completed chapter + appendix in `bios-prelim.cls` and
get the real page count.**
Don't estimate past this point — Project 1 and the literature review both
hit real, non-obvious page-count effects from format changes alone, and
this step replaces every estimate above with an actual number.
→ *Verify:* PDF renders; body and appendix page counts recorded.

**5. Derive CTJ from the finished chapter, by hand, condensing to
3,500 words / 6 exhibits.**
Work from the chapter's finished prose, not from the old CTJ file's
prose — where the two would say the same thing, prefer the chapter's
version and compress it, rather than reverting to CTJ's already-compressed
phrasing (which may have compressed content the chapter's review in step
0 changed). Where the old CTJ prose already found a better-compressed
phrasing for something the chapter says at more length, reusing it is
fine — the point is that the chapter is the arbiter when they disagree,
not that the old file is worthless.

Everything cut in this step goes into a supplementary-material file
(`draft/trimmed-material.md` or equivalent — same convention already used
for the literature review) rather than being discarded, per `ROADMAP.md`'s
own instruction. Confirm whether that cut material becomes CTJ's own
online Supplementary Information submission, replaces the current
`Supplementary_Information.tex`, or sits alongside it — decide this once
step 5 is underway and it's clear how much of the current SI survives
unchanged versus needs rebuilding from the chapter's condensed leftovers.
→ *Verify:* CTJ renders within 3,500 words / 6 exhibits with a small
margin, not flush against the cap; every claim in CTJ traces to the
chapter; nothing cut is silently lost.

**6. Retire the old `CTJ_Manuscript.tex` / `Supplementary_Information.tex`
content once step 5's derived versions are stable** — mark superseded
rather than delete, consistent with this dissertation's "preserve cut
material" convention, since the old files may still hold reference value
for phrasing even after they stop being edited independently.
→ *Verify:* the derived CTJ + SI are what's live going forward; the old
files are clearly marked as superseded, not silently left ambiguous.

**7. Copy the finished chapter + appendix into bios-dissertation.**
Becomes Chapter 3 (Project 2) in `prelim/project-proposals/project2-incidence/`,
following the same directory/symlink/master-bib pattern already used for
Chapter 2 (Project 1) — `bios-prelim.cls` symlinked in, `bibliography:`
pointed at the shared `prelim/references.bib`, not a local copy. Per the
component-fate table above, write a `manuscript-declarations.md` alongside
it (matching `project1-spillover/manuscript-declarations.md`'s form)
holding CTJ's abstract and full Declarations block verbatim, and confirm
whether Funding + Data Availability get the same unnumbered-section
treatment Amber's precedent uses.
→ *Verify:* renders cleanly against the master bib; no duplicate or
conflicting citekeys introduced; nothing dropped from the chapter is
undocumented (either in the appendix or in `manuscript-declarations.md`).

## Open items carried forward, not resolved by this plan

- Whether CTJ's own application section gets the proposal-style reframing
  before submission, or waits for the real dataset — depends on timing
  not yet known, and on whether the real data arrives before or after
  CTJ is otherwise ready to submit.
- The real SUDDEN dataset / IRB access itself — external dependency,
  tracked in that repo's own CLAUDE.md, not something this plan can move
  forward.
- Whether the derived Supplementary Information (step 5) replaces or
  supplements the current `Supplementary_Information.tex` — deferred to
  when step 5 is actually underway.
