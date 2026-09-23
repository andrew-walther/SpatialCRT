# Step 0 review brief — `Dissertation_Chapter.qmd` citation / numeric / voice pass

Companion to `manuscript-unification-plan.md`, step 0 ("Citations and
voice"). Written 2026-09-23. Paste the block below into a fresh agent/session.
It is a **review-only** brief: the reviewer reports findings, the executing
session decides and applies fixes. Re-run it unchanged after fixes are
applied; step 0 is done when a fresh run returns no unresolved findings.

Paths below are absolute because the brief spans two repos.

---

```
You are reviewing one file and reporting findings. Do NOT edit any file.

Target: /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/paper/dissertation_chapter/Dissertation_Chapter.qmd
(Project 2 of a BIOS PhD dissertation: 6 treatment-assignment designs for
spatial cluster randomized trials; it will become Chapter 3 of the prelim.)

Ignore the \begin{abstract}...\end{abstract} block; it is being removed.

## Authorities — use these, not your own memory or the web
- Citations: /Users/ajwalther/GithubProjects/bios-dissertation/prelim/references.bib
  is the verified master bibliography. Do NOT look anything up on CrossRef,
  Google Scholar, or the web. For "does the source support the claim", use the
  master entry's title/journal/year/abstract/note fields, plus how the same key
  is used in the already-reviewed literature review:
  /Users/ajwalther/GithubProjects/bios-dissertation/prelim/literature-review/draft/literature-review-draft.qmd
  If the bib + lit review can't settle it, report it as UNVERIFIABLE rather
  than guessing.
- Numbers: /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/results/six_design_manuscript/
  six_design_summary.txt, dissertation_results_extract.txt,
  application_table_6design.txt. These are ground truth. Check every number
  in the prose and every table cell, including derived claims ("less than
  0.03", "more than four times", "13 of 15", counts like 9,600 / 384 / 320,
  "every design", "every configuration").
- Method facts (what the code actually does): /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/code/
  01_spatial_setup.R ... 05_run_simulation.R, and the project CLAUDE.md at
  /Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/CLAUDE.md.
  Check design descriptions (Table 1), neighbor definitions, replication
  counts, deterministic-design handling, incidence generation, estimator
  specification against the code.
- Project 1 claims: the chapter describes "Project 1" (small-grid proof of
  concept). Check those claims against
  /Users/ajwalther/GithubProjects/bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.qmd
  (model name, grids, designs, what its discussion anticipated).

## What to check
1. CITATIONS: every @key resolves in the master bib; each citation supports
   the specific claim it's attached to (not just the topic); attribution is
   correct (e.g., who "formalized" or "introduced" something).
2. NUMBERS: as above. Flag rounding only if it changes a stated comparison.
3. METHOD / FACT ACCURACY vs code and Project 1 draft.
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
   review draft — person ("we" vs. "this chapter"), how the dissertation
   refers to other chapters/projects (e.g., "Project 1 of this dissertation"
   vs. "Chapter 2"), tense. Report inconsistencies; don't rewrite.
8. SOURCE-ONLY FACTS: claims that no file above can verify (e.g., pilot-site
   geography, dataset sizes, IRB status). List them separately as
   NEEDS-AUTHOR-CONFIRMATION with the exact sentence; do not judge them.

## Output
A numbered list. For each finding: line number(s), category (1-8), severity
(ERROR = factually wrong / contradicts ground truth; WARN = overclaim,
unsupported, inconsistent; STYLE = voice/mannered), the exact quoted text,
the evidence (file + value), and a one-line suggested direction (not a
rewrite). End with a count per severity. If a category has no findings, say
so explicitly so silence isn't ambiguous.
```
