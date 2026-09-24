# Plan: Write Project 2 long-form first; derive CTJ + supplementary
# material from it

**Written in a bios-dissertation session on 2026-09-18** (that's where
this planning conversation happened) and copied here; this is the copy of
record.

**Correction made during drafting, same session:** an earlier version of
this plan had the direction backward — CTJ + its Supplementary
Information as the primary source, with the dissertation chapter derived
from it. That matched Project 1's pattern, but Project 1 is the repo's
own *named exception* (its chapter is copy/paste from an already-accepted
manuscript, precisely because it couldn't be rewritten). `ROADMAP.md`'s
"Manuscript & Chapter Writing Strategy" already settles this for Project
2, and had before this session touched it: **write one long-form
version, unrestricted in length; the 3,500-word CTJ submission is a
condensed derivative, with material cut for length preserved in a
supplementary-material file rather than discarded.** This version of the
plan follows that.

**Reviewed adversarially by a fresh Opus pass, 2026-09-18, same session
(after this plan was first written by a Sonnet session).** 16 findings,
all applied below — mostly operational gaps (a no-op step, a missing bib
step, no deadline checkpoint, figures never addressed) and two places
where the plan mischaracterized its own cited precedent. Nothing found
undermined the chapter-first decision itself; the review was explicitly
scoped not to re-litigate that.

## Deadline

The prelim package is due to the committee **2026-11-16**. Only steps
0–4 and 7 below feed that deadline — **steps 5 and 6 (deriving and
retiring the old CTJ files) are explicitly post-deadline work**; CTJ
submission is a soft target (`ROADMAP.md`'s own Per-Project Status table),
not something to spend November on at the expense of the chapter. Target
steps 0–4 and 7 complete by **~2026-11-01** to leave room for full prelim
assembly (lit review + Chapter 2 + Chapter 3 + front matter as one
document, which has never been attempted) before the 11/16 deadline.

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

**Known and accepted:** the BIOS ~40 pp. combined guideline for the
prelim (lit review + proposal, excluding bibliography) is already
exceeded by the lit review (22 pp.) + Chapter 2 (42 pp. body + 4 pp.
appendix) alone — ~68 pp. before Chapter 3 exists. Chapter 3's length is
not constrained by that guideline, per the same user decision already
recorded in `chapter2-integration-plan.md`. This is worth restating here
because a plan that only ever discusses *not being too short* could
otherwise read as implying no ceiling exists at all; the ceiling exists
and has already been knowingly exceeded.

**Advisor feedback, 2026-09-23 — revises step 3's balance, not the
hierarchy itself.** An advisor reviewing this project emphasized that the
dissertation chapter should carry the bulk of the relevant content, so a
reader isn't sent back and forth between chapter and appendix for this
project specifically. This doesn't change the four-deliverable hierarchy
above or the chapter-first decision — it changes where step 3's
body-vs-appendix line falls. The original test ("does a reader need this
to follow the argument, or would they only consult it?") was already
correct in principle but had been applied too readily toward appendix for
genuinely relevant results; step 3 below is revised to bias toward body.
The appendix is now reserved for what's mechanical rather than merely
detailed — see step 3.

## Build section by section, not as one pass — components have
different fates in different documents

Not every piece of a manuscript belongs in all four deliverables. The
clearest existing precedent is Project 1's own port (see
`prelim/project-proposals/notes/chapter2-integration-plan.md`): the
journal abstract and the full `Declarations` block (ethics, consent,
competing interests, author contributions, funding, data availability)
came **out of the dissertation chapter entirely for the prelim**,
preserved verbatim in a side file (`manuscript-declarations.md`) rather
than deleted. **All seven were dropped, with no exceptions** — Amber
Young's own dissertation keeps two of the seven (Funding, Data
Availability) as unnumbered sections after the chapter's Discussion, but
that was recorded as a *thesis-stage option*, not adopted for the prelim
port. Read `manuscript-declarations.md` directly before treating this
table as settled; it also warns that reinstating a Data Availability
Statement later requires a real edit (the code's actual public location),
not a mechanical copy.

Project 2 needs the same per-component decision, made explicit before
drafting rather than discovered by accident later. Working table, to be
corrected as each component is actually drafted:

| Component | Chapter body | Chapter appendix | CTJ manuscript | CTJ supplementary |
|---|---|---|---|---|
| Standalone abstract | **No** — remove the chapter's current ~430-word abstract explicitly (step 0); the dissertation's own front matter carries its abstract | — | **Yes** — but SAGE's guidelines call for an **unstructured abstract ≤250 words**, and the current draft is a ~300-word structured one; this needs rewriting, not trimming (double-check against the live journal page, not this plan, since guidance pages go stale) | — |
| Introduction / Methods / Results / Discussion | Yes, in full | Overflow only (per the body-vs-appendix test in step 3) | Yes, condensed | Overflow from condensing |
| Full parameter-grid tables, robustness breakdowns | No | Yes | No | Yes |
| Application section | Yes, reframed as a proposal (step 2) | — | Placeholder, until real data | — |
| Ethics / consent / competing interests / author contributions | No | No | **To be written**, if and when SAGE requires it — CTJ's own Declarations currently has only "Declaration of conflicting interests" and "Funding," no ethics/consent/authorship content exists to copy from anywhere; verify against SAGE's actual checklist, don't assume Project 1's BMC set applies | No |
| Funding / Data Availability Statement | **No for the prelim** — matching Chapter 2, which dropped all seven; Amber's two-unnumbered-section treatment is a thesis-stage decision to make for Chapters 2 and 3 together, not a default to apply here alone | — | Yes (CTJ already has a Funding section) | — |
| Acknowledgments | **No** — Project 1's chapter dropped this too (it thanked the editor/reviewers, meaningless outside a journal submission) | — | Yes | — |
| Title page / author list + affiliations / keywords / running head | N/A — dissertation front matter handles this | — | Yes — mostly already present in `CTJ_Manuscript.tex`, confirm nothing's missing | — |
| Cover letter | N/A | N/A | **Yes, but not a manuscript component** — a separate submission requirement this table can't track; add to a submission checklist when step 5 is actually underway | N/A |

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
open and never done. `Dissertation_Chapter.qmd` (~7,305 words, `wc -w` on
the raw source including YAML/markup, ~89% of that is actual prose) is
~90% a restatement of content already in `CTJ_Manuscript.tex` (body word
count doesn't reproduce cleanly across counting methods — estimates from
~2,100 to ~2,540 depending on whether captions/markup are included;
**recount with one stated method before relying on the 3,500-word margin
in step 5**) + `Supplementary_Information.tex` (**15 pp.**, complete, 12
tables + 10 figures) — the two genuinely new items anywhere in the long
version, both already present, not something to add, are a coverage-by-ρ
table for Checkerboard (its Table 6) and one sentence describing the
pilot intervention itself.

**Decision:** `Dissertation_Chapter.qmd` becomes the primary, long-form
source — expanded and completed, not just reformatted. `CTJ_Manuscript.tex`
and its Supplementary Information become **derived** from it, condensed by
hand to fit the journal's 3,500-word / 6-exhibit limit. The existing CTJ
prose isn't wasted — it's a strong first draft of already-compressed
phrasing for several claims and is worth reusing where the compression it
already found is good — but it stops being the source of truth once the
long-form chapter is complete.

**Sizing, calibrated against real renders, not guessed — body only, the
appendix has its own target below:**
Project 1's chapter rendered at 42 pp. body + 4 pp. appendix (13,902
words of source prose by the same counting method, ≈302 words/page in
`bios-prelim.cls`). The current `Dissertation_Chapter.qmd` is already
close in scale (~7,305 words, ≈24 pp. of *body* at that ratio) and is the
better starting point than CTJ's ~8-page version — expanding a document
that's already ~90% of the way to a full argument is far less work than
expanding CTJ's compressed version twice (once to a full chapter, then
trimming a *different* subset back down for the journal). Project 2's
simulation study is also narrower in scope than Project 1's (one design
comparison vs. two spillover mechanisms across a larger design space), so
landing shorter than Project 1's 42 pp. body is expected, not a shortfall
to correct.

**Both halves of this estimate move with the 2026-09-23 advisor
feedback, in opposite directions from what an earlier version of this
plan expected.** Step 3's body-biased test keeps most of the per-
incidence-mode/per-parameter detail in the body rather than the
appendix, so the ≈24 pp. body figure above is a floor, not the expected
final count — expect it to grow once that content moves in. The appendix
correspondingly shrinks: rather than the "well above 4 pp." this plan
estimated before the advisor feedback (reasoning from the SI's own 15
pp./22 exhibits), expect something closer to Project 1's 4-pp./no-floats
scale, since most of the SI's content is now staying in body. Give both
their own explicit target once step 3 is underway rather than trusting
either number in advance.

## Steps

**0. First, snapshot today's versions before touching anything — this
plan edits `Dissertation_Chapter.qmd` in place from here on.**
Copy today's `paper/dissertation_chapter/Dissertation_Chapter.qmd`,
`paper/ctj_manuscript/CTJ_Manuscript.tex`, and
`paper/ctj_manuscript/Supplementary_Information.tex` into this repo's
existing `archive/` folder (it already holds `OutcomeIncidenceDesign_Legacy`,
`PreliminarySpatialSim`, `SpatialSim_Unified` — same convention, one more
entry), dated. This is a safeguard in addition to git history, not a
replacement for it — the point is a version someone can find without
archaeology. Nothing described in this plan deletes content; by the end
there are four live deliverables (chapter body, chapter appendix, CTJ
manuscript, CTJ supplementary material) and the current CTJ/SI additionally
persist here as a dated snapshot once step 6 marks the live copies
superseded.

**Then give `Dissertation_Chapter.qmd` the review pass it's never had, and
fix its bibliography and abstract before adding anything else.**

Three things, all before any new content goes in:

- **Citations and voice.** Check every citation resolves and supports
  its claim, no mannered prose, voice consistent with the rest of this
  dissertation. Use `prelim/references.bib` as the citation authority —
  it is already verified; do not re-derive citations from CrossRef or
  similar, which is how a stale entry entered `IncidenceDesign_shared.bib`
  in the first place (see the bibliography point below). For numeric
  claims, treat `results/six_design_manuscript/{six_design_summary,
  dissertation_results_extract,application_table_6design}.txt` as ground
  truth. Also apply `CLAUDE.md`'s standing rule that manuscript prose
  never mentions DIM. Consider writing a companion execution prompt for
  this step the way Chapter 2 had one, with an explicit "report findings,
  don't fix them yourself yet" brief.
- **Bibliography.** Repoint `Dissertation_Chapter.qmd`'s bibliography at
  `prelim/references.bib` (the shared master) rather than
  `IncidenceDesign_shared.bib`. Three keys the chapter currently cites
  diverge between the two files, and the master's version is the correct
  one in each case — `baird_optimal_2018` (the shared bib has the
  superseded working-paper metadata; the master has the corrected
  published-article version, deliberately fixed during Chapter 2's own
  bibliography work) and `leung_rateoptimal_2022` /
  `leung_clusterrandomized_2025` (differing `pages` fields). Do not
  concatenate the two `.bib` files — Chapter 2's own plan warned that
  colliding keys make pandoc silently keep one arbitrarily. All 23 keys
  the chapter currently cites already exist in the master, so this is a
  pointer change today, not a merge; that stops being true once step 3
  imports SI content, so do this now while it's still free.
- **Abstract.** Remove the chapter's current ~430-word `\begin{abstract}`
  block from the body — move its content into a
  `manuscript-declarations.md` file (see step 7) rather than deleting it.
  The dissertation's own front matter carries its abstract; a per-chapter
  one is redundant with it, matching Chapter 2.

→ *Verify:* a fresh reviewer pass returns no unresolved findings;
`Dissertation_Chapter.qmd` renders against `prelim/references.bib` with
no missing keys; the in-body abstract is gone and its content is saved,
not discarded.

**Step 0 progress (2026-09-24).** Done: the snapshot (`archive/IncidenceDesign_Manuscripts_PreUnification_2026-09-23/`),
the review brief (`step0-chapter-review-prompt.md`), and the first review pass
(`step0-review-findings-2026-09-23.md`: 17 errors, 31 warnings, 17 style). Not
done: the bibliography, abstract, and prose fixes. **Two changes to this step,
by author decision:** (1) the bibliography is one self-contained file in this
repo, `paper/SpatialCRT_IncidenceDesign.bib` (renamed from
`IncidenceDesign_shared.bib`, with the three diverging entries corrected from
`prelim/references.bib`), not a pointer into bios-dissertation; step 7
checks entry-by-entry agreement with the master before porting. (2) The
review showed that several results come from simulation oversights, not
wording, so Methods and Results are rewritten only after step 0.5.

**0.5. [Added 2026-09-24] Revise the simulation and re-run it before
rewriting Methods/Results.** Full plan, decisions, and verification:
`simulation-revision-plan.md`. In short: designs and outcomes use the
same incidence surface; Poisson clusters get 100,000 people; ties are
broken randomly; seeding is explicit; Monte Carlo SEs are honest;
aliasing is flagged; a non-oracle sensitivity estimator is added; a
validated fast estimator is used. Then downstream results are
regenerated and step 0's remaining prose work continues against the new
ground truth. The chapter-first order of steps 1–7 is unchanged.

**Step 0.5 progress (2026-09-24).**
- **Done:**
  - Phase A: bib, abstract, and the results-independent prose fixes.
  - Phase B: spec, code, tests, pilot, and the full re-run. It's verified, and a 1%
    `lagsarlm` cross-check passed.
- **Decisions after the pilot:**
  - Keep the 6 designs.
  - Balanced Quartiles treats exactly 50.
  - Framing: which designs work under heterogeneous incidence, not whether knowing
    incidence is necessary.
- **In progress:** Phase C (downstream regeneration).
- **Next:** Phase D, which needs the author's choice of the one recommended design.

**1. Confirm the two items the chapter needs are already present — this
is a verification step, not drafting.**
Both the ρ-stratified Checkerboard coverage table and the
pilot-intervention description sentence identified during the
CTJ-vs-chapter comparison are **already in `Dissertation_Chapter.qmd`**
(Table 6; the "workforce development and disaster-response-style
training" sentence near the application section). Confirm both are still
present and accurate after step 0's edits — don't re-add them, and don't
skip this check just because nothing needs adding.
→ *Verify:* both located in the current text; neither was accidentally
removed by step 0.

**2. Reframe the application section as a forward-looking proposal.**
The scaffolding already exists almost verbatim in the current text:
dataset description (SUDDEN algorithm, 2018–2021, NC death certificates,
IRB status), planned covariates, geography (58 CC service areas, Lenoir
County pilot site + Kinston/Jones satellites), and the hoped-for insight
(process measures, not incidence rate, on a pilot timescale). This is a
**deletion-and-reframing edit**, not new drafting: cut the synthetic
MSE/coverage numeric run and the placeholder-caveat paragraph around it
(the chapter's application section currently has no dedicated
figure/table of its own to remove — its three figures are all results
figures; the application *maps* are handled separately in step 3, and not
all of them carry the placeholder problem this step is fixing — see
step 3); keep everything else; shift the surrounding prose
from "placeholder result shown, to be superseded" to "analysis plan
stated, not yet run." Mention that the R ingestion pipeline already
exists and is waiting on data access, so a reader knows this isn't
stalled for lack of a plan. This applies to the long-form chapter now;
CTJ's own application section stays as-is (still a placeholder) until the
real dataset is available and CTJ is actually being prepared for
submission — the two documents don't need to be resolved on this point
at the same time.
→ *Verify:* the section reads as an honest in-progress item with no
fabricated or placeholder numbers.

**If the real SUDDEN dataset arrives before this plan is finished
(genuinely possible — access was described as imminent as of
2026-09-23):** don't wait for it, and don't restart this step once it
does. The proposal-framing content above (dataset description,
covariates, geography, pilot site, hoped-for insight) is not throwaway
scaffolding — it is the setup paragraph the real-results version will
also need almost verbatim. If data lands mid-plan, replace the "not yet
run" framing with the actual results in the same section, in place, as a
small targeted update once the real analysis is done — this is cheaper
than delaying the whole chapter to write the application section once,
under deadline pressure, and it's exactly why this step exists rather
than waiting.

**That backfill is a batch of new content, not a numbers swap — expect
it to need its own pass through step 3.** As of 2026-09-23, only one
application figure exists that isn't tied to placeholder data —
`community_college_service_area_clusters.png`, showing which counties
cluster into each NC community college's service area. Once the real
study runs, expect at minimum:

- a **county-level** real incidence map — genuinely new, not a
  replacement of anything that exists even as a placeholder today. The
  current synthetic pipeline (`synthetic_incidence_map.png`) only ever
  produces incidence at the **cluster** (58 CC-service-area) level, so
  this is an additional exhibit, not a swap.
- a **cluster-level** real incidence map, the direct real-data analogue
  of today's `synthetic_incidence_map.png` — this one *is* a like-for-
  like replacement.
- **design-comparison results applied to that real cluster-level
  surface** — the real-data analogue of the existing six/eight-design
  MSE/coverage comparison, run against actual incidence instead of the
  synthetic placeholder.

Confirm this list against whatever the real study actually produces
rather than treating it as exhaustive — it's a floor, not a ceiling, on
what to expect. When that happens, run the application section's own
figures and tables through step 3's body-vs-appendix test and the
exhibit-fate table the same way the rest of the chapter's content already
went through it, rather than dropping them in wherever is convenient
because the rest of the chapter is by then already finished.

**3. Decide body vs. appendix for the chapter's own supplementary
content and figures alike, using a real test, not a page-ceiling-driven
default — biased toward body, per the 2026-09-23 advisor feedback above.**
A generous word budget is a reason the chapter body *can* carry more
detail than CTJ's SI does; for this project specifically, an advisor has
also said it *should* — the chapter is meant to be readable without
flipping to the appendix, so "detailed" is no longer by itself a reason
to push something out. The test, revised: **does the appendix version
exist only to spare the reader mechanical detail, or does it contain a
result the reader needs to trust or understand the chapter's central
claims?** Only the former goes to appendix now; "exhaustive" and
"irrelevant to the argument" are no longer treated as the same thing.

*Prose and tables:*
- **Body — including detail the earlier version of this test would have
  pushed out:** the headline results and the comparisons that support
  the chapter's central claims (MSE/coverage ranking, the omnibus
  Friedman result, the reframed application section from step 2), *and,
  per the advisor feedback,* the per-incidence-mode / per-parameter (ρ,
  γ, spillover regime) / τ-level breakdowns — these are genuine results
  a reader would otherwise have to leave the chapter to see, not
  reproducibility mechanics. The 8-design comparison is a judgment call
  under the revised test: it extends past the chapter's main 6-design
  scope, but if it's presented as a robustness check *of* the 6-design
  argument rather than a separate result, that argues for keeping a
  compact version in body too, with only the full 8-design table in
  appendix. Decide this one explicitly and record which way it went.
- **Appendix — narrowed to what's genuinely mechanical, not just
  detailed:** seeding protocol, exact metric formulas, and other
  reproducibility material a reader consults to verify the pipeline
  rather than to understand a result. Much of this can come from
  `Supplementary_Information.tex` directly rather than being rewritten —
  it's already built for this role. If in doubt whether something is a
  "result" (body) or "mechanics" (appendix), default to body per the
  advisor's direction — the cost of an appendix that's slightly too thin
  is much lower here than the cost of a chapter a reader can't follow
  without leaving it.

*Figures — not yet addressed anywhere before this plan; do this
explicitly, since figures don't sort into body/appendix the same way
prose does. Three application-related figures, not two, exist across
`application/report/figures/` and `ctj_manuscript/figures/
application_maps/` — inventory both locations, not just
`results/six_design_manuscript/`:*
- The chapter's current figures are **stale relative to CTJ's** — verify
  by checksum against `results/six_design_manuscript/` before reusing
  anything the chapter currently has; at least one figure predates a
  reordering the chapter's own prose already assumes.
- The **bias-variance decomposition** figure exists only in
  `ctj_manuscript/figures/` and is not discussed anywhere in the current
  chapter text. If it belongs in the body (it explains *why* Checkerboard
  fails, which is argument-critical, not supplementary), that requires
  writing new prose around it, not just placing an existing figure — be
  honest that this is new writing, not a placement decision, when
  scoping the work.
- **Two of the three application figures are placeholder-tainted; one is
  not — don't treat them as a single group.** `synthetic_incidence_map.png`
  and `kmeans_regions_map.png` are built from the synthetic placeholder
  data step 2 is reframing away from; give either the same
  placeholder-caveat treatment as the application section's prose if
  imported, or leave them out until real data exists.
  `community_college_service_area_clusters.png`
  (`application/report/figures/`) is different in kind — it shows the
  real geography (which counties feed which CC service area), is not
  data-dependent, and doesn't need a placeholder caveat at all. It is
  currently used only in CTJ's Supplementary Information at reduced
  width (`0.75\linewidth`) and doesn't appear in the chapter at all —
  add it to the chapter's application section (or its appendix) at full
  size regardless of when the real dataset arrives.
- **A fourth consideration, not covered by the body-vs-appendix content
  test above: physical size and exhibit budget, independent of
  relevance.** `community_college_service_area_clusters.png` is
  3600×2400px with a legend dense enough that CTJ's SI already renders
  it at reduced width, and it would likely need simplifying further (or
  omitting, with a cross-reference to the chapter/SI) to fit legibly
  inside CTJ's 6-exhibit, narrow-column budget — while the same figure
  at full size and full legend fits comfortably in the chapter's
  single-column dissertation layout with no such constraint. This is a
  format decision, not a content-relevance one: a figure can be
  argument-critical (per the test above) and still need an abbreviated
  version for CTJ specifically, purely because CTJ's page format can't
  hold what the chapter's can. Apply this same check to every other
  figure/table as it's inventoried, not just this one.
- Fold in the repo's own already-flagged, still-open to-do: a considered
  review of the full main-text + SI figure list to decide what to
  keep/drop/combine (`CLAUDE.md` describes the current set as "a quick
  fit," not a deliberate final selection) — and use this same pass to
  surface other results/figures/tables from `results/` or `application/`
  that never made it into any manuscript but are worth the chapter
  appendix, now that the appendix is no longer assumed to be minimal.
- Build a parallel exhibit-fate table (mirroring the component-fate table
  above) once the figure review above is done, so the appendix's exhibit
  count and page estimate (below) rest on a real list, not a guess. Give
  it a fifth column, or a note per row, for the format-constraint check
  above — a figure's fate can differ between "which document" and "what
  size/version" independently.

*Judgment call, decide once and record it:* whether the chapter appendix
holds all of the SI's 22 exhibits or a curated subset. With the revised,
body-biased test above, expect the appendix to land closer to Project
1's 4-pp./no-floats scale than the "well above 4 pp." estimate this plan
carried before the 2026-09-23 advisor feedback — most of the SI's content
is being kept in body now, not moved to appendix wholesale. Give it an
explicit target once this step's figure review is done rather than
assuming either scale in advance.

For each candidate — prose or figure — write down which side of the test
it landed on and why, one line each, so the reasoning survives, not just
the outcome.
→ *Verify:* re-render after each batch; running body and appendix page
counts recorded against this step's own stated target, not the body-only
estimate in Context; every figure reused is checksum-verified against its
canonical source; each placement decision has its one-line reasoning
somewhere (a source comment is sufficient).

**4. Render the completed chapter + appendix in `bios-prelim.cls` and
get the real page count.**
Don't estimate past this point — Project 1 and the literature review both
hit real, non-obvious page-count effects from format changes alone, and
this step replaces every estimate above with an actual number.
→ *Verify:* PDF renders; body and appendix page counts recorded.

**5. [Post-deadline — see "Deadline" above.] Derive CTJ from the finished
chapter, by hand, condensing to a body in the ~3,000–3,300 word range
(not maximizing toward 3,500), plus a real abstract rewrite, not a trim.**

**Target word count, set explicitly (2026-09-23): ~3,000–3,300 words of
body, not the cap itself.** The current CTJ draft sits at roughly
2,400–2,540 words — meaningfully under-using the 3,500-word budget SAGE
allows. Once this step has the finished chapter's full argument to draw
from, use enough of that headroom to produce a substantive manuscript,
while stopping short of the cap on purpose: the goal is a robust
submission with real room to add content or clarification during peer
review, not one so close to 3,500 that a single reviewer request forces
a cut elsewhere to make space. Landing flush against the cap defeats
that purpose as much as landing at 2,400 does — both leave no room to
respond to review.

**The exhibit limit is 6 tables and figures combined**, already verified
against SAGE's *Clinical Trials* author guidelines (see the Context
section above) — same number, whether counted as "figures" or exhibits
generally; there's no separate, larger figure-specific allowance on top
of it.

Work from the chapter's finished prose, not from the old CTJ file's
prose — where the two would say the same thing, prefer the chapter's
version and compress it, rather than reverting to CTJ's already-compressed
phrasing (which may have compressed content the chapter's review in step
0 changed). Where the old CTJ prose already found a better-compressed
phrasing for something the chapter says at more length, reusing it is
fine — the point is that the chapter is the arbiter when they disagree,
not that the old file is worthless. The abstract needs separate handling
per the component-fate table: SAGE's current guidelines call for an
unstructured abstract ≤250 words, and the existing draft is a longer,
structured one — confirm the live requirement before finalizing, and
rewrite rather than trim it into shape.

Everything cut in this step goes into a supplementary-material file
(`draft/trimmed-material.md` or equivalent — same convention already used
for the literature review) rather than being discarded, per `ROADMAP.md`'s
own instruction. Confirm whether that cut material becomes CTJ's own
online Supplementary Information submission, replaces the current
`Supplementary_Information.tex`, or sits alongside it — decide this once
step 5 is underway and it's clear how much of the current SI survives
unchanged versus needs rebuilding from the chapter's condensed leftovers.
Also assemble the non-component submission checklist flagged in the
component-fate table (cover letter, any journal-specific forms) — this
plan can name that it's needed but can't fill it in this far ahead of
submission.
→ *Verify:* CTJ body lands in the ~3,000–3,300 word range (using a
word-count method stated explicitly this time, not the ~2,400–2,540
estimate range this plan couldn't pin down earlier) — not flush against
3,500 and not left under-using the budget the way the current draft is;
exactly 6 exhibits or fewer; abstract ≤250 words, unstructured, confirmed
against the live journal page; every claim in CTJ traces to the chapter;
nothing cut is silently lost.

**6. [Post-deadline.] Retire the old `CTJ_Manuscript.tex` /
`Supplementary_Information.tex` content once step 5's derived versions
are stable** — mark superseded rather than delete, consistent with this
dissertation's "preserve cut material" convention, since the old files
may still hold reference value for phrasing even after they stop being
edited independently.
→ *Verify:* the derived CTJ + SI are what's live going forward; the old
files are clearly marked as superseded, not silently left ambiguous.

**7. Copy the finished chapter + appendix into bios-dissertation.**
`prelim/project-proposals/project2-incidence/` currently holds only a
pointer `.qmd` — no `draft/` directory, no `bios-prelim.cls` symlink, no
`figures/` — so this step is scaffold-then-port-then-render, not a
straight copy:

- Scaffold the directory the way Chapter 2's was scaffolded, reusing its
  verify checklist directly (adapted to `CHAPTER 3:` / `APPENDIX B:`):
  the chapter heading reads `CHAPTER 3:` at 2 in from the top; any
  `\setcounter` for chapter numbering is placed *after* `\mainmatter` or
  it silently reads as Chapter 1; the class's `CSLReferences`
  patch is documented FRAGILE with only prior use to test against, so
  render early and often rather than assuming it will just work.
- Convert the citation machinery: the current chapter uses
  `cite-method: natbib`, but the prelim class expects citeproc-style
  `@key` citations — this needs converting, not just a bibliography
  pointer change.
- Convert or resolve the executable `knitr::kable` R chunk (around lines
  113–125 of the current file) — either inline it as a static raw-LaTeX
  table (consistent with everything else being static prose by this
  point) or confirm R is actually available at render time in this
  environment.
- Symlink `bios-prelim.cls` in and point `bibliography:` at the shared
  `prelim/references.bib`, not a local copy (already done as of step 0,
  just confirm it survived the port).
- Per the component-fate table, write a `manuscript-declarations.md`
  alongside it (matching `project1-spillover/manuscript-declarations.md`'s
  form) holding the removed abstract (from step 0) and CTJ's existing
  Declarations content verbatim — currently just "Declaration of
  conflicting interests" and "Funding," not the full seven-part set
  Project 1's had, since CTJ's Declarations section is that short as of
  now.

→ *Verify:* renders cleanly against the master bib; no duplicate or
conflicting citekeys introduced; the chapter heading and numbering match
Chapter 2's checklist exactly; nothing dropped from the chapter is
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
- The exact figure-fate list (step 3) and the CTJ submission checklist
  (step 5) — both explicitly deferred to when those steps are underway,
  not filled in speculatively here.
