# Codex prompt: Phase 2g, prose review of the curated SI

Author: Andrew Walther. Written 2026-11-25. Paste everything below the line into a fresh Codex
session opened in `SpatialCRT/projects/IncidenceDesign`.

---

You are running Phase 2g of `docs/plans/2026-11-25-project2-phased-plan.md`: the prose review of the
curated Supplementary Information (SI) for the Clinical Trials (CTJ) paper on sampling designs for
spatial cluster randomized trials. Do only 2g. Items 2h (final rebuild and page counts) and 2i
(Chapter 3 alignment) belong to a later session; collect what 2i will need, but do not do it.

## Read first

- `AGENTS.md` (standing rules; follow every one) and the 2026-11-25 entries in `DECISIONS.md`.
- The hand-over from Phases 2a–2f:
  `paper/ctj_manuscript/revisions/2026-11-25-si-curation/` (every sign-off file, the decisions at
  the end of each, `apply_2c.py`, `apply_2d.py`).
- `docs/reviews/2026-11-24-si-review.md`, sections (c) and (d), and
  `docs/reviews/2026-11-24-si-claim-inventory.md` ("Cross-section duplications"). They predate
  2a–2f, so their line and S-numbers are stale; find each passage by its text. Items already
  fixed in 2c–2e (recommended-design wording, SRS in the tests, "identifiable", stale statements,
  figure labels, Table S18's minus sign, the S7/S12 merge) need no second pass.
- The human-writer skill: `~/.claude/skills/human-writer/SKILL.md`. Apply it to every sentence you
  propose.
- The CTJ itself (`paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`),
  read only, so the SI complements it rather than repeating it.

## How the SI is built

- One source, two builds. `paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex`
  builds the **curated** SI that is submitted (43 pages now). `Supplementary_Information_Full.tex`
  sets `\SIFullVersion` and builds the **unabridged** SI (64 pages). Material that only the
  unabridged build prints sits between `\iffullsi` and `\fi`; curated replacements use
  `\iffullsi <unabridged text>\else <curated text>\fi`.
- Review the **curated** text (what prints in `Supplementary_Information.pdf`). Text shared by both
  builds changes in both, which is fine. Leave unabridged-only blocks alone unless a shared edit
  forces a change there.
- The author's rule: the curated SI and the CTJ are companions, and the curated SI is
  self-contained. It must never mention or point to the unabridged version (treat that version as
  if it did not exist; cut material instead). It may cite the article's own exhibits through the
  `M-` labels (for example `Figure~\ref{M-fig:nc-mse} of the article`).
- The CTJ prints SI numbers through `\ref{SI-<label>}`, and the SI cites the article through
  `\ref{M-<label>}` (both `xr-hyper`). Never hard-code an SI or article number. Build order:
  in `supplementary_info/` run `latexmk Supplementary_Information.tex` and
  `latexmk Supplementary_Information_Full.tex`; in `manuscript_ctj/` run `latexmk CTJ_Manuscript.tex`
  and `latexmk CTJ_Reading.tex`; then build the curated SI once more.

## What to do

Work in three stages and stop for the author's sign-off after each. Never edit SI, CTJ or Chapter
3 text without sign-off on the exact wording.

1. **Human-writer pass on the curated SI.** The SI grew by accretion: caveats repeated in up to
   15 places, 400-word paragraphs, pipeline jargon ("stored fits", "manifest", "production
   allocations", "aggregate output tables", "Spill"), prose that describes a figure instead of
   letting the caption do it, and hedges stacked three deep. State each caveat once where it
   belongs and refer back after that. Split and cut the long paragraphs. Keep every number, claim
   and citation. If a sentence's number has to change, stop and ask.
   Known items from the review to cover: "fitted" used as a verb (use "fit"; "fitted" only before
   a noun); abbreviations used before they are defined (SRS, ML, MDE, "Jaccard overlap");
   "strategies" versus "designs" (use "designs", with SRS named separately); prose colons (avoid
   them unless they are clearly the best framing); colon-led captions and the "Case study: …"
   section titles; Table S5's whole-number counts printed as "27.00"; Table S20's row label
   "Baseline / adjusted"; the SI 5.2 rounding sentence that now overlaps the sentence after it;
   one sentence that the case study's "2×2 Blocking" is the spatial-blocking adaptation the CTJ
   names. The curated overview keeps "the full methods" (author decision).
2. **Two independent reads of the curated SI after stage 1 is applied**: one as a biostatistician
   (correctness of methods, inference, Monte Carlo precision, whether each claim follows from its
   table or figure) and one as a clinician or epidemiologist planning a trial (can they follow the
   case study and use the planning guidance). Do the two reads separately, so that neither sees
   the other's notes, and do not anchor on stage 1's reasoning. Each read lists findings; you then
   turn the ones worth acting on into one sign-off list.
3. **Apply** approved items by ID with a script (pattern: `apply_2c.py`; each edit must match
   exactly the expected number of times, or the script stops), rebuild in the order above, and run
   the checks.

## Sign-off format (author rules)

- Batches of about 8 to 12 IDs, with your recommendation first and a one-line reason.
- Every item shows the full current text and the full proposed text, not a description. Put the
  lists in `paper/ctj_manuscript/revisions/2026-11-25-si-curation/2g-signoff.md` (stage 1) and
  `2g-review-signoff.md` (stage 2), and record the author's decisions at the end of each file.
- Style: American English; "However" and "Therefore" open their sentences; careful articles; never
  a bare "NC" as a descriptor ("the NC case study" or "North Carolina"); design names exactly as
  in AGENTS.md; MSE and MSE-derived values to three decimals, ratios to SRS to two; notation β on
  **Z**, τ on **X**, ψ for spillover (the code still says `tau`, `beta`, `gamma`).

## Checks after every applied batch

- `python3 paper/tools/verify_manuscripts.py` must PASS. It also fails on any undefined reference
  in either SI build. Its notice that `tab:adjacent-pairs` is pending is expected until 2i.
- `Rscript code/tests/test_manuscript_exhibits.R` must pass.
- `python3 paper/tools/ctj_wordcount.py` stays at 3,435 unless the CTJ changed.
- `pdftotext` the curated PDF: no "??", and no "unabridged" in any case.
- Report failures exactly as printed.

## Rules

- Numbers come from saved results or executed code, never retyped. A figure change goes through
  the script that produces it, never a hand edit.
- Stage and commit only your own paths, by name (never `git add -A`, `git add .` or
  `git commit -a`). Other sessions work in this checkout. Commit after each applied batch with an
  imperative summary and an itemized body. Author: Andrew Walther; no AI attribution anywhere,
  including commit messages. Do not push, and do not tick boxes in the phased plan.
- Never edit `paper/ctj_manuscript/manuscript_ctj/revision_markup/`, the review package in
  `revisions/2026-11-24-round6-word-limit/for-lin/`, `code/01`–`05`, or anything in
  `projects/SpillSpatialDepSim/`. Do not edit `paper/thesis_chapter/`.
- Keep a running list of SI changes that Chapter 3's appendix will need, appended to the 2i list
  in the shared hand-over: add a section "2g changes for Chapter 3" to
  `paper/ctj_manuscript/revisions/2026-11-25-si-curation/2g-signoff.md`. The earlier 2i items
  (C3-1 to C3-6) are in `.../2c-signoff.md`, `2d-signoff.md` and the 2e files.
- Leave the two new-content slots empty: 2.2 and a results section before the case study (Phase
  3, Incidence-Guided Saturation sensitivity) and 3.2 (Phase 4e, estimating the incidence
  coefficient τ). They are drafted elsewhere.

Done means: stage 1 and stage 2 sign-off lists decided and applied, both SI builds clean, all
checks passing, the Chapter 3 list written, and a short report of commits, check results and the
curated SI's page count.
