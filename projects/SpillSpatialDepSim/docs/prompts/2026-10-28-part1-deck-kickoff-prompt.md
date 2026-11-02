# Kickoff prompt: Project 1 talk deck (prelim Part I)

Paste everything below the line into a new Claude Code session opened in
`/Users/ajwalther/GithubProjects/SpatialCRT/projects/SpillSpatialDepSim`.

---

## Task

Two connected pieces of work for my dissertation Project 1, "Spatial cluster randomized
trials: spillover and spatial dependence" (SpillSpatialDepSim):

1. **Put the corrected Chapter 2 into the prelim documents** in bios-dissertation, so the
   chapter text and the slides agree.
2. **Plan and build a roughly 20-minute talk deck with a full speaker script.** It is
   **Part I of my prelim exam presentation** (Monday 2026-12-07), the first of three
   ~20-minute parts (Parts I/II/III, about an hour of talk in a 90-minute exam). The deck
   must also work **standalone**, since a trimmed copy becomes Dr. Lin's conference talk
   in December.

Work in FAST mode:
- read first;
- ask me the open questions in one batch;
- present a numbered plan covering both pieces;
- wait for my approval before changing files.

**Audience:** the prelim committee, four statisticians and one non-statistician. This deck
can be more technical than my Part II deck. Show the model and the estimation properly,
while keeping the motivation and take-home clear for the non-statistician.

## Piece 1: corrected Chapter 2 into the prelim

The corrected chapter is ready in `review/thesis-chapter/`. Its `README.md` has the copy
steps and the follow-up checks. I approved the integration on 2026-10-28. The end state
is:

- **`bios-dissertation/prelim/project-proposals/project1-spillover/draft/`** holds the
  corrected chapter in the thesis format.
- **`project1-spillover/draft-uncorrected/`** holds the current, uncorrected chapter. Move
  it there first, so nothing is lost.
- **`bios-dissertation/prelim/full-draft/`** is rebuilt so its Chapter 2 is the corrected
  version. It keeps no uncorrected copy; that copy lives only in `draft-uncorrected/`.

Steps:
- Move the current `draft/` contents to `project1-spillover/draft-uncorrected/`.
- Copy in the corrected `.qmd`/`.tex`/`.pdf`, the two replaced 3×3/3×4 figures and
  `manuscript-declarations.md`, as the README describes. Carry over every other file
  `draft/` needs to render, such as the figures and the class file.
- Re-render the chapter in `draft/` and confirm it is still 52 pages with no unresolved
  references.
- Rebuild the full prelim with `python3 prelim/full-draft/tools/build_prelim.py`. Check
  the build log, the page count, and that Chapter 2 in `prelim.pdf` shows the corrected
  3×3/3×4 tables.
- Note the correction and `draft-uncorrected/` in `project1-spillover/README.md`.

**Checks**

- `prelim/project-proposals/notes/check-prose-fidelity.py` will now report the intended
  differences. That exception is recorded in bios-dissertation `DECISIONS.md` (2026-10-25).
- Update bios-dissertation `ROADMAP.md` (Chapter 2 status).

**Funding name fixed (2026-10-28)**

"National Institutes of Mental Health" was corrected to "National Institute of Mental
Health" in `review/manuscript/Revisions_V2c_corrected.tex` and
`review/thesis-chapter/manuscript-declarations.md`, so the copied declarations carry the
fix.
- Check whether the chapter `.qmd` or the full draft has its own copy of the funding text.
- Re-render the corrected manuscript PDF and its latexdiff markup in `review/manuscript/`
  so they match.

## Piece 2: the Part I deck

### Where things stand

- **Outline.** `presentation/OUTLINE.md` (Draft 2, 2026-10-25) has a visual version
  `outline.html`, a `figures/` folder, `tools/make_figures.R` and
  `tools/build_outline_html.py`. Treat it as a strong starting point to review, refine
  and improve, not a final plan. Its "Open questions for Andrew" section lists six items.
  Bring the ones my answers below don't settle to me in your first batch, with your own.
- **Source of truth.** Use the *corrected* manuscript
  `review/manuscript/Revisions_V2c_corrected.tex` and the corrected Chapter 2 for content
  and numbers.
- **Read before planning.**
  - `AGENTS.md`, `README.md` and `../../AGENTS.md`;
  - my memory notes on Project 1 (BSS framing, the correction, and the rule that
    Chapter 2 is final);
  - `bios-dissertation/slides/prelim/README.md`, for the prelim deck structure and theme.

### Model to follow: the Part II deck (built 2026-10-28)

Study it before planning: `../IncidenceDesign/presentation/` (`README.md`, `OUTLINE.md`,
`sudden_talk_2026-10-29.qmd`, `theme/`, `tools/`), plus its reviews in
`../IncidenceDesign/docs/reviews/2026-10-28-sudden-talk-*`. Reuse what worked:

1. **Format.**
   - Quarto revealjs at 1920×1080, rendered to one self-contained HTML file (fonts,
     figures and MathJax embedded) so it runs offline.
   - Render with RStudio's Quarto
     (`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`), because
     `/usr/local/bin/quarto` is broken on this machine.
   - **S** opens the speaker view.
2. **Theme.** Copy `unc-gillings.css`, `prelim-overrides.css` and `footer-part.html`
   unchanged from `bios-dissertation/slides/prelim/test/`. Put anything specific to this
   deck in `theme/deck.css`, so Parts I and II look identical when the prelim deck is
   assembled.
3. **Sections.** Use the same five, written into each slide's footer through
   `data-footer-part`, with no divider slides: **Background**, **Methods**,
   **Simulation Results**, **Application**, **Discussion**, then **Backup**.
4. **Slides.**
   - Each slide title is a claim written as a sentence.
   - Figures carry no in-figure titles or notes, and they print true hyphens.
   - Where a manuscript figure is dense, make a slide-sized version with larger text from
     the frozen results, using scripts in `presentation/tools/`.
5. **Model notation slide.** Show the full spatial lag model and a glossary that defines
   every symbol: y_i, α (baseline), β (intervention effect), x_i, ψ (spillover), z_i, ρ
   (spatial dependence), W with w_ij (row-normalized rook neighbors) and ε_i.
   - Show side by side how the binary spillover indicator z_i is defined under each
     spillover type (Control Clusters Only; Control & Intervention Clusters).
   - State the level the model works at, as the manuscript describes it: subjects are
     sampled within clusters and the model uses cluster-level, subject-averaged
     quantities. Check the exact wording against the manuscript.
   - Define MSE, bias and variance on a slide or in the glossary.
   - Backup can carry the reduced form and the likelihood.
6. **Speaker notes.**
   - A verbatim read-aloud script, one sentence per bullet, with per-slide timing that
     totals about 20 minutes.
   - Backup slides get short "use when asked" notes.
   - Write it with the `human-writer` skill: flowing sentences, few colons, no sentence
     opening with "Because".
7. **Title slide.** List me only.
8. **Backup section.** Generous, since Q&A is long. Candidates:
   - the 3×3 and 3×4 results;
   - where BSS falls among random allocations on every grid;
   - full result tables;
   - reduced form and likelihood;
   - the weights matrix;
   - spillover versus spatial dependence;
   - why inference is model-based;
   - buffer-zone designs;
   - structural versus causal estimands;
   - the simulation algorithm.
9. **Outline kept in sync.** Regenerate the `OUTLINE.md` slide entries from the deck
   source (adapt `outline_from_deck.py` from Part II) and rebuild `outline.html`. List
   every file and the build commands in a `presentation/README.md`. Keep the outline's
   both/prelim/conference tags.
10. **Review loop.**
    - After the first full render, run two fresh, independent reviewers in parallel:
      - a **statistician**: technical accuracy; every symbol defined; numbers traced to
        the corrected manuscript, Chapter 2 or `results/`; whether the claims hold up
        under committee questions;
      - a **non-statistician**: motivation, clarity, and what the take-home means for
        probation trials.
    - Fix their findings, then repeat with fresh reviewers until neither reports a
      blocking issue. Part II needed four rounds.
    - Save the reports in `review/` with dated names.
    - If a reviewer objects to wording I chose, keep my wording and raise it with me.
11. **Prelim adaptation notes.** End `presentation/README.md` with an "Adapting this deck
    for the prelim (Part I)" section, modeled on the one in
    `../IncidenceDesign/presentation/README.md`. Check it against the prelim skeleton,
    `bios-dissertation/slides/prelim/test/slide-deck-test/prelim-test.qmd`. The section
    covers:
    - **Standalone frame:** which slides to drop (title, outline, thank-you). The
      acknowledgements go on the single prelim closing slide after Part III, grouped with
      Parts II and III. The paper reference goes on the last Part I Discussion slide.
    - **Divider and footers:** the Part I divider; footer labels prefixed
      `Part I · <Section>`; menu titles prefixed `I:` (author approved, 2026-10-28).
    - **Carry-over tags:** keep, rewrite or drop each slide by its
      both/prelim/conference tag.
    - **Background:** identify the content most relevant to Project 1, since Part I
      comes first and sets up spillover, spatial dependence and the spatial lag model.
      Part II then builds on it without repeating it.
    - **Backup:** slides moved under the `Backup · Part I` divider. The prelim has three
      backup sections, one per part, after Part III ends.
    - **Theme and figures:** keep the theme consistent with Parts II and III, merge
      `deck.css` without class clashes against Part II's, and copy the figures.
    - **Timing:** about 20 minutes per part, about 60 in total. Parts I and II may run a
      little longer while Part III is still brief.
12. **Ask early.** Batch content decisions before building. My answers forced a third
    Part II draft, and asking first avoids that.

### My decisions already made (2026-10-28)

- **The correction.** Don't point it out anywhere in the slides. The corrected 3×3 and
  3×4 results are simply the results. Leave out a correction backup slide and any
  correction status line.
- **Publication.** Reference the BMC Medical Research Methodology paper at the end of the
  Part I Discussion (with the implications, limitations and future work). The prelim has a
  single thank-you/Q&A slide at the very end, not one per part.
- **Standalone closing slide.**
  - For standalone use, the deck still ends with a "Thank you / Questions" slide, tagged
    so it can be dropped from the prelim.
  - It acknowledges Dr. Feng-Chang Lin (advisor) and Dr. Tonya Van Deinse (collaborator)
    and states the funding as the corrected declarations do: National Institute of Mental
    Health K01MH129619 (Dr. Van Deinse), and the Bureau of Justice Assistance
    15PBJA-24-GG-01896-PRJH award to the NC Department of Adult Correction.
  - Note in the outline that these acknowledgements move to the shared prelim closing
    slide.
- **Notation.** It differs between Parts I and II by design, because the models differ.
  Part II adds incidence and models at the cluster level; Part I uses subject-averaged
  cluster quantities. Each part defines its own model and notation in full. No
  harmonizing, and no bridge slide.
- **Framing.**
  - BSS is a good design, neither the best nor the worst, compared with "rolling the dice"
    on one SRS allocation. Lead with where BSS falls among the random allocations.
  - Be plain that under spillover into both arms most random allocations beat BSS,
    because in a checkerboard spillover exposure is the exact opposite of treatment.
  - Per-scenario "optimal allocation" counts are not the message.
  - The phrase "checkerboard paradox" is retired.
- **Prelim structure.** The prelim deck has a brief shared motivation, then each project
  with its own prior work and gap, and no standalone literature review. Ask me whether
  Part I, coming first, should carry some of that shared motivation.

## Constraints

- **Frozen files.** Project 1's `code/`, `results/` and `paper/` are frozen. In `review/`,
  only the re-renders named above are allowed. Everything else new goes in
  `presentation/`.
- **Authorship.** Script headers credit "Author: Andrew Walther" only. No Claude or Codex
  mentions in the deck or in documents.
- **Git.** Commit or push nothing in either repo without my OK. Copy the deck into
  `bios-dissertation/slides/prelim/part-1/` only after I approve it.

## Done means

- `project1-spillover/draft/` and the full-draft prelim both carry the corrected Chapter 2,
  the uncorrected copy is in `project1-spillover/draft-uncorrected/`, both re-render
  cleanly, and the README and ROADMAP are updated.
- The deck:
  - `OUTLINE.md` is revised;
  - it renders to self-contained HTML with footer section names, about 20 minutes of main
    slides, and a backup section;
  - the speaker view carries the full script with timing;
  - the review reports are saved, and the last round has no blocking issues;
  - `presentation/README.md` lists the files and build commands and ends with the prelim
    adaptation section.
- You report the paths and any questions left for me.
