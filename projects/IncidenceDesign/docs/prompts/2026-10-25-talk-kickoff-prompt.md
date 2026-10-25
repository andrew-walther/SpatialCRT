# Kickoff prompt: 20-minute project talk (Thursday 2026-10-29)

Paste the block below into a new Claude Code session started in
`SpatialCRT/projects/IncidenceDesign/`.

---

I'm giving a 20-minute talk on Thursday 2026-10-29 about the IncidenceDesign project
(Spatial CRT, Project 2 / Chapter 3). The audience is clinicians and epidemiologists,
so it should be relatively non-technical. The deck should also be a good first version
of "Part II" of my prelim deck, so it can later be adapted and extended for a technical
biostatistics audience (prelim on 2026-12-07). Help me plan it, then build it.

Scope: a full overview of the project, based on the thesis chapter
(`paper/thesis_chapter/Thesis_Chapter.qmd`) and the CTJ manuscript draft
(`paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex`), covering in order:
1. Motivation and background literature.
2. The models used (outcome model, spatial dependence, spillover), kept to what this
   audience needs.
3. The simulation study and the allocation designs considered.
4. The NC SUD application study.
5. Key conclusions and implications, limitations, and future or open work.

Weight the time and the discussion toward the application study: this audience is best
placed to give feedback on it (the SUD-prevention education intervention, scaling from
a single-county pilot to a statewide community-college trial, colleges needed, partial
participation, incidence-guided allocation). Model detail beyond the essentials goes to
backup slides that the prelim version can promote.

Read first:
- `AGENTS.md` (standing rules: design names, "Saturation Regions", no single
  recommended design, β = 0 primary with the education endpoint) and `DECISIONS.md`.
- `docs/findings/2026-10-21-methods-extensions-report.md` (Sections 2 and 5) and the
  application report, `application/report/IncidenceDesign_Application_Study_Report.qmd`.
- `bios-dissertation/slides/prelim/OUTLINE.md`, Part II, including the "Methods
  extensions" entries P2-X1 to P2-X13. This talk should reuse and feed that outline.
- `bios-dissertation/slides/prelim/notes/` (planning kickoff, reference-deck synthesis)
  and `bios-dissertation/templates/` for the slide style and format.

Constraints:
- Data: the NC data and findings are fine to share with this group (the approval for
  them comes from this group). The restricted-data handling rules in
  `application/AGENTS.md` still apply to files and code.
- Location: build this deck in the project repo, in the `presentation/` folder (the
  placeholder `presentation/OUTLINE.md` is there for the outline), using the same style
  and format as the bios-dissertation prelim template. It is later copied or adapted
  into `bios-dissertation/slides/` for prelim Part II; this deck is the exception to the
  "slides live in bios-dissertation" rule.
- New figures from methods-extension step 7a (power vs number of colleges,
  colleges-needed bar chart, spread-vs-gain scatter) may or may not exist yet; check
  `docs/TODO.md` and tell me if the talk is waiting on them.
- Follow my global workflow: present a numbered plan and wait for approval before
  building anything.

A first draft outline already exists (2026-10-25): `presentation/OUTLINE.md` (19 main
slides, 8 backups, each with title, figure, slide text, notes, minutes and a
prelim/technical/talk-only tag) and `presentation/outline.html`, a visual version with
the figures inline, rebuilt from the markdown by `python3 presentation/tools/build_outline_html.py`.
Its PNG figures are in `presentation/figures/`. Start by going through its "Open
questions for Andrew" with me, revise the outline from my answers, and only then build
the deck.
