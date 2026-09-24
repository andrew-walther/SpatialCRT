# Step 0 review brief — `Dissertation_Chapter.qmd` citation / numeric / voice pass

Companion to `manuscript-unification-plan.md`, step 0 ("Citations and
voice"). Written 2026-09-23; updated 2026-09-24 (Phase D block 2) to check against
the revised simulation (full re-run 2026-09-24). Changes: numbers ground truth is now
the revised-run summaries plus read-only R on the 20260924_025509 .rds files
(`application_table_6design.txt` dropped: it is labeled STALE); method authority is
`docs/plans/simulation-revision-spec.md`, with the code as secondary; citations are
checked against the chapter bib, with the bios-dissertation master still authoritative
for shared keys; "Project 1" is now "Chapter 2"; the abstract line is removed (the
abstract is gone); a "Decisions already made" section lists settled author choices.
Paste the block below into a fresh agent/session.
It is a **review-only** brief: the reviewer reports findings, the executing
session decides and applies fixes. Re-run it unchanged after fixes are
applied; step 0 is done when a fresh run returns no unresolved findings.

Paths below are absolute because the brief spans two repos.

---

```
You are reviewing one file and reporting findings. Do NOT edit any file.

Target: /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/paper/dissertation_chapter/Dissertation_Chapter.qmd
(Chapter 3 of a BIOS PhD dissertation: 6 treatment-assignment designs for
spatial cluster randomized trials.)

## Authorities — use these, not your own memory or the web
- Citations: the chapter's bibliography is
  /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/paper/SpatialCRT_IncidenceDesign.bib
  (read by the chapter through the symlink paper/dissertation_chapter/shared-refs.bib).
  Every @key must resolve there. For keys that are also in
  /Users/ajwalther/GithubProjects/bios-dissertation/prelim/references.bib, that
  file remains the verified master: the two entries must agree, and the master
  settles any disagreement. For keys only in the chapter bib (e.g., the five
  test citations friedman_use_1937, nemenyi_distributionfree_1963,
  demsar_statistical_2006, wilcoxon_individual_1945, holm_simple_1979), check
  that the entry's metadata is plausible and that the citation supports the
  claim; flag problems, don't reject the key. Do NOT look anything up on
  CrossRef, Google Scholar, or the web. For "does the source support the
  claim", use the entry's title/journal/year/abstract/note fields, plus how the
  same key is used in the already-reviewed literature review:
  /Users/ajwalther/GithubProjects/bios-dissertation/prelim/literature-review/draft/literature-review-draft.qmd
  If the bib + lit review can't settle it, report it as UNVERIFIABLE rather
  than guessing.
- Numbers (ground truth), under
  /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/results/:
  six_design_manuscript/six_design_summary.txt,
  six_design_manuscript/dissertation_results_extract.txt,
  eight_design_supplementary/eight_design_summary.txt,
  sim_data/full_run_verification.txt. You may also run read-only R
  computations on sim_data/*_20260924_025509.rds and
  sim_data/surface_results_*_20260924_025509.rds (do not write any file).
  Do not use application_table_6design.txt: it is labeled STALE. Check every
  number in the prose and every table cell, including derived claims ("less
  than 0.02", "13 of 15", counts like 9,600 / 1,920 / 320 / 160, "every
  design", "every configuration").
- Method facts: /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/docs/plans/simulation-revision-spec.md
  is the method authority. The code
  (/Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/code/
  01_spatial_setup.R ... 05_run_simulation.R) and the project CLAUDE.md at
  /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/CLAUDE.md
  are secondary. Check design descriptions (Table 1), neighbor definitions,
  replication counts, deterministic-design handling, incidence generation,
  seeding, and estimator specification.
- Chapter 2 claims: the chapter describes Chapter 2 (small-grid proof of
  concept). Check those claims against
  /Users/ajwalther/GithubProjects/bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.qmd
  (model name, grids, designs, what its discussion anticipated).

## What to check
1. CITATIONS: every @key resolves in the chapter bib (and agrees with the master where shared); each citation supports
   the specific claim it's attached to (not just the topic); attribution is
   correct (e.g., who "formalized" or "introduced" something).
2. NUMBERS: as above. Flag rounding only if it changes a stated comparison.
3. METHOD / FACT ACCURACY vs the revision spec (code secondary) and the Chapter 2 draft.
4. INTERNAL CONSISTENCY: claims that contradict each other within the file;
   table/figure numbering and references (Table N / Figure N as rendered in
   order of appearance); captions that don't describe their figure (open the
   PDFs in paper/dissertation_chapter/figures/ if you can, otherwise report
   the caption as unverified).
5. OVERCLAIMS: conclusions the evidence doesn't support (e.g., a Friedman
   rejection at each tau shows designs differ at each tau, not that the
   *ranking* is stable).
6. STANDING RULES:
   - Manuscript prose never mentions DIM / difference-in-means as a comparator.
   - Designs are named descriptively, never "Design 1"/"D8" shorthand.
   - No mannered prose, as defined in /Users/ajwalther/GithubProjects/bios-dissertation/CLAUDE.md
     ("Remove all mannered prose" section). Also flag padded/inflated
     phrasing ("remarkably", "it is worth noting", stacked hedges,
     sentences that restate the previous sentence).
7. VOICE: compare against the Chapter 2 draft (path above) and the literature
   review draft — person ("we" for actions, "this chapter" for the
   document), references to other chapters ("Chapter 2", never "Project 1"),
   tense (past for results). Report inconsistencies; don't rewrite.
8. SOURCE-ONLY FACTS: claims that no file above can verify (e.g., pilot-site
   geography, dataset sizes, IRB status). List them separately as
   NEEDS-AUTHOR-CONFIRMATION with the exact sentence; do not judge them.

## Decisions already made — do not flag as errors
These are settled author choices. Report a finding only if the text states
one of them inaccurately or inconsistently, not because you would decide
otherwise.
- The chapter recommends Incidence-Guided Saturation Quadrants. Balanced
  Quartiles is presented as another design that works well (marginally better
  only under the both-arms, i.e. bidirectional, spillover regime).
- Queen contiguity is primary and rook contiguity is the sensitivity analysis;
  numbers are never pooled over neighbor types without a "pooled" label.
- Under rook contiguity, Checkerboard's tau is not identified (WZ = 1 - Z).
  This is linked to Chapter 2 through its collinearity result
  beta-hat ~ beta - psi. "Checkerboard paradox" is used only as Chapter 2
  used it: amplified bias under bidirectional spillover with spatial
  dependence.
- Queen Checkerboard's high MSE is variance (WZ = 1/2 at interior clusters).
- Framing: the question is which designs work under heterogeneous incidence,
  NOT that "incidence knowledge is needed".
- Saturation Quadrants being indistinguishable from Incidence-Guided
  Saturation Quadrants is a finding, not a problem.
- "Block Stratified Sampling" and "Checkerboard" are interchangeable names;
  the chapter states this once.
- Coverage is about 0.94 for every design except Checkerboard under rook, and
  it does not separate the designs.
- Links between Chapters 2 and 3 are kept brief by author choice; don't flag
  missing cross-chapter connections.
- Aliased rook Checkerboard estimates are kept and flagged, by decision.
- The Application section will be reframed in a later step (step 2), and its
  placeholder numeric run is expected to be removed then. Report Application
  findings, but mark each one "step 2".

## Output
A numbered list. For each finding: line number(s), category (1-8), severity
(ERROR = factually wrong / contradicts ground truth; WARN = overclaim,
unsupported, inconsistent; STYLE = voice/mannered), the exact quoted text,
the evidence (file + value), and a one-line suggested direction (not a
rewrite). End with a count per severity. If a category has no findings, say
so explicitly so silence isn't ambiguous.
```
