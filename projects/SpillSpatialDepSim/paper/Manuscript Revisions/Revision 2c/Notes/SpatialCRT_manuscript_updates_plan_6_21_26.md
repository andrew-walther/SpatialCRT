# SpatialCRT Manuscript Updates — Implementation Plan (Revision 2c)

**Date:** 2026-06-21
**Driver:** 6/19/2026 Zoom with Feng-Chang Lin + handwritten markup on `Revisions_V2b_markup.pdf`
**Source notes:** `Revision 2c/Notes/manuscript_revision_notes_06_19_2026.md`
**Working folder:** `Revision 2c/`
**Baseline being edited:** `Revisions_V2b.tex` (the submitted V2b manuscript) → becomes `Revisions_V2c.tex`
**Annotated-diff baseline (unchanged):** `Revisions_V1.tex`

---

## 1. Goal

Apply a targeted notation/clarity pass to the model section and appendix, add an explicit ICC (intra-cluster correlation) statement + limitation, and fix punctuation/indentation nits. Then refresh the annotated diff and the response letter, and have a reviewer subagent confirm V2c vs. V2b before handoff.

This round is **notation and framing only** — no new scientific results, no changes to figures, tables, or simulation numbers.

## 2. Decisions locked (from 6/21 clarifying questions)

1. **Versioning:** `V2c` suffix. Folder `Walther_SpatialCRT_LaTeX_Revisions_V2c/`, files `Revisions_V2c.tex/.pdf`, `Revisions_Annotated_V2c.tex/.pdf`, response letter `response_letter_V6`. Annotated diff still baselined against `Revisions_V1.tex`.
2. **Short track:** Ignore/drop this round. The `_short` folder and `response_letter_V5_short.*` are removed/shelved during cleanup; only the full version gets edits.
3. **Simulation verification gate — ✅ resolved (see Phase 1):** the R code is **individual/point-level** (per-individual iid error), *not* cluster-level. Adopted compromise: present cluster-level `y_i` for reader clarity with the honest "shared cluster mean, independent iid error, no within-cluster correlation" wording — **never claim individuals share an identical outcome**. ICC is therefore not parameterized (limitation + future work).

## 2a. Execution quick reference (operating manual)

**Run in FAST mode.** This plan is the authoritative spec — work the phases 0→9 in order.

**Read first:**
- This plan.
- Meeting notes: `Revision 2c/Notes/manuscript_revision_notes_06_19_2026.md`
- Memory: `project_revision2c_clusterlevel_decision.md`, `project_revision2b_status.md` (locked C1–C6 — do not regress), `project_latexdiff_workflow.md`

**Locked — do NOT re-ask (resolved 6/21):**
- **Versioning V2c.** Rename folder → `Walther_SpatialCRT_LaTeX_Revisions_V2c/`; files → `Revisions_V2c.tex/.pdf`; annotated → `Revisions_Annotated_V2c`; letter → `response_letter_V6`.
- **Short track dropped.** Archive 2b artifacts (session handoffs, both `v5*_letter_build.js`, `_short` folder/letter, old C1–C6 diffs) to `Revision 2c/_archive_2b/` — archive, don't delete. `Revision 2b/` stays untouched.
- **Sim is individual/point-level** (verified `02_simulation_core.R`: `n_pts×n_pts` W, `rnorm(n_pts)` iid error). Present cluster-level `y_i` for clarity but state honestly: modeled at individual level; within a cluster the response & covariates are identical so the ONLY within-cluster variation is iid error; no within-cluster error correlation. **Never write "identical outcomes."**
- **Model equation = Option (a):** displayed equation purely cluster-level with a single `ε_i`; individual detail in prose only.
- **ICC:** iid errors, no cluster random effect → ICC not parameterized. Limitation (§5.2) + future work (§5.3, add a cluster random effect).
- **D3 = Option A:** delete the redundant displayed estimation equation; keep §2.7 as lean prose; repoint `\ref{Equation:SARestimation}` → `Equation:SARmodel`; drop NO `\cite` keys.
- **D4 indentation:** `\usepackage{indentfirst}` for first-paragraph-after-heading; `\noindent` on mid-sentence post-equation continuations.
- **Response letter V6:** start from the AUTHOR-EDITED `response_letter_V5.docx` in `Revision 2b/revision response/` (canonical; `v5_letter_build.js` is STALE — do not re-run). Edit in place, preserve salutation/addressing verbatim (no name changes), add items additively.

**Phase map:** 0 setup/cleanup → 1 sim-gate (✅ resolved) → 2 notation/equation restructuring → 3 ICC → 4 punctuation/indent → 5 appendix numbering → 6 compile + grep/citation gates → 7 annotated diff (vs. V1) → 8 response letter V6 → 9 reviewer subagent (V2c vs V2b).

**Deliverables:** (1) `response_letter_V6.{docx,pdf}`; (2) `Revisions_Annotated_V2c.{tex,pdf}` (V1→V2c); (3) `Walther_SpatialCRT_LaTeX_Revisions_V2c/` with `Revisions_V2c.{tex,pdf}` + `references.bib` + `figures/`.

---

## 3. Key reframing — and the tension it resolves

The 6/19 meeting **reverses the prior Revision 2b "Option B" decision** (which retained individual-subject notation `y_{ik}`, `y_{jk'}` with a clarifying sentence). Feng-Chang now wants the main-text model written at the **cluster level** (`y_i`, drop `k` and `k'`), because:

- The reduced/equilibrium form `y = (I−ρW)⁻¹(α + Xβ + Zψ + ε)` and the weights matrix `W` are already **n×n cluster-level** objects (n = number of clusters). The `y_{ik}` individual notation was inconsistent with them.
- Modeling the cluster-level outcome sidesteps ICC: all individuals in a cluster share the same structural outcome, so within-cluster correlation is not part of the estimand.

**Tension — RESOLVED by Phase 1 (see below).** The code is in fact individual/point-level (iid per-individual error), so the meeting note's literal "all individuals share the same outcome" is *not* what the code does. Author-directed compromise: present the model at the cluster level (`y_i`) for reader clarity, while honestly stating that sampled individuals share the cluster mean and differ only by an independent error. This keeps the locked C2/C4 "within-cluster variation" prose correct and makes the ICC limitation precise (errors iid → no ICC parameterized).

---

## Phase 0 — Folder setup & cleanup (Revision 2c)

The `Revision 2c/` folder is already a copy of `Revision 2b/`. Clean it and rename to the V2c convention.

**Rename (manuscript working folder + files):**
- `Walther_SpatialCRT_LaTeX_Revisions_V2b/` → `Walther_SpatialCRT_LaTeX_Revisions_V2c/`
- Inside it: `Revisions_V2b.tex` → `Revisions_V2c.tex` (this is the file we edit). Delete the stale compiled `Revisions_V2b.{pdf,aux,log,out}` and `Revisions_Annotated_V2b.*` (regenerated fresh as V2c).
- Keep: `references.bib`, `figures/`, and `Revisions_V2b_markup.pdf` (rename → `Revisions_V2b_markup_REFERENCE.pdf`; it is the edit source — keep until implementation is verified, then it can be archived).

**Drop (2b-specific, no longer needed in 2c):**
- `Notes/session_1..6_handoff.md` (2b session log)
- `Notes/v5_letter_build.js` and `Notes/v5_short_letter_build.js` (stale build scripts — V6 is now produced by editing the author-edited V5 docx in place, not by re-running these)
- `Walther_SpatialCRT_LaTeX_Revisions_V2b_short/` (entire short-version folder)
- `revision response/response_letter_V5_short.{docx,pdf}` (short track)
- `Diffs/pdfs/` + `Diffs/snapshots/` + `Diffs/short/` 2b per-item C1–C6 artifacts (clear; V2c per-item diffs start fresh)

**Keep:**
- `Notes/manuscript_revision_notes_06_19_2026.md` and this plan
- `revision response/response_letter_V5.{docx,pdf}` — the author-edited V5 letter (V6 is built by editing this in place; see Phase 8)
- `Original & Annotated Manuscripts/.../Revisions_V1.tex` (annotated-diff baseline — **do not touch**)
- `journal comments/` (reviewer/editor docx — reference)
- `Diffs/create_clean_diff.py` only if still wanted as reference (latexdiff is the actual tool)

**Verify:** `ls` the renamed folder; confirm `Revisions_V2c.tex` opens and that `figures/` + `references.bib` are present so it can compile.

> **D1 RESOLVED (author, 6/21):** archive — move the 2b-specific artifacts (session handoffs, short-track folder/letter, old per-item C1–C6 diffs) to `Revision 2c/_archive_2b/` rather than deleting. The original `Revision 2b/` folder also stays in place as the canonical 2b reference.

---

## Phase 1 — Simulation-code verification gate — ✅ RESOLVED (2026-06-21)

Code read: `~/GithubProjects/SpatialCRT/projects/SpillSpatialDepSim/code/UnifiedSpatialSim/02_simulation_core.R` (lines 121–142).

**Finding (the code is INDIVIDUAL/POINT-level, not cluster-level):**
- `W` is built `n_pts × n_pts` (160/180/240 points), with block structure: every point connects to every point in neighboring districts.
- `epsilon <- rnorm(n_pts, mean=0, sd=sd)` — one **independent** error per individual. No cluster-specific random effect.
- `linear_resp = α + β·intervention + ψ·spillover + ε`; `response = solve(I − ρW, linear_resp)` — reduced form at the point level.
- Fit: `lagsarlm(response ~ intervention + spillover, listw = mat2listw(W, "W"))` — fit on individual points.
- Individuals in a cluster share the same cluster-level covariates (`intervention`, `spillover`) → the **same structural mean** — but differ through their independent error draw. They are **not identical**.

**Adopted resolution (author-directed compromise, 6/21):** present the model at the **cluster level (`y_i`) for reader clarity**, while being **honest** that the simulation samples individuals who share the cluster mean and differ only by an independent error term. Do **not** claim individuals have identical outcomes.

**Implications for prose (locks the wording in Phase 2/3):**
- "where" paragraph: `y_i` = cluster-level outcome (modeled quantity); add one honest sentence — individuals in a cluster share the same structural mean and differ only via an independent error.
- The locked **C2 sentence about within-cluster variation contributing to the error-variance estimate is CORRECT per the code → keep/align it**, do not delete (reverses my earlier concern).
- ICC statement is precise: errors are iid with no cluster random effect, so ICC is not parameterized; adding a cluster random effect is the future extension.
- `W` presented as cluster-level `n×n` is a clarity abstraction (code builds it point-level). Acceptable; do not over-claim identical outcomes alongside it.

---

## Phase 2 — Notation & equation restructuring (core change)

All edits in `Revisions_V2c.tex`, §2.5–§2.7 (current lines ~275–327) and Appendix A/B (~726–793).

### Target structure of §2.5 "Spatial Autoregressive (SAR) Model"

**N2.1 — Drop the "general form" block (old Eq 1 + old Eq 2).**
- Delete current Eq (1) general SAR (`Equation:SARmodel`, lines 279–281) **as a separate general equation**, and delete the matrix-form Eq (2) (lines 285–287) and its lead-in "This can also be expressed in matrix notation as:".
- Per markup p13 (strikeouts on "This can also be expressed in matrix notation as" and on `y = ρWy + Xβ + ε`) and notes ("drop the general section with equations 1 & 2").

**N2.2 — Promote the study model to the primary equation (new Eq 1).**
- Old Eq (3) (lines 289–291) becomes the first equation, carrying the `\label{Equation:SARmodel}` (so the C1 estimands `\ref{Equation:SARmodel}` on line 298 stays valid).
- Rewrite at **cluster level** and in the **preferred term order** (α + linear predictor + spatial dependence last):

  `y_i = α + β x_i + ψ z_i + ρ \sum_{j \in \mathcal{N}(i)} w_{ij} y_j + ε_i`

  - Drop `k` and `k'` everywhere (`y_{ik}`→`y_i`, `y_{jk'}`→`y_j`, `ε_{ik}`→`ε_i`). **Option (a) confirmed:** the displayed equation is purely cluster-level with a single cluster-level error `ε_i`; the individual-level detail lives in prose, not the equation.
  - `\sum_{j \neq i}` → `\sum_{j \in \mathcal{N}(i)}` (markup p12 "j in N(i)").
- **Retain the introductory framing prose** (do not delete it with old Eq 1). Keep a lead-in such as: *"Outcomes in the presence of spatial dependence and spillover effects can be considered via modeling with a spatial autoregressive model\cite{...}. For the design considered in this study, the spatial autoregressive model for the cluster-level outcome has the form:"* — so the section still motivates the model before presenting it. Preserve the three existing citations from the deleted intro sentence.

**N2.3 — Add the structural matrix/vector form (Feng-Chang's "vector form of Eq 2 with the Zψ term").**
- Immediately after new Eq (1), add the matrix analog (unnumbered `equation*` — not cross-referenced):

  `y = α + Xβ + Zψ + ρWy + ε`

- Lead-in: "In matrix form, …". This replaces the deleted bare `y = ρWy + Xβ + ε` with the *full* model including spillover.

**N2.4 — Rewrite the "where" definition paragraph (cluster-level `y_i`, no k/k', with honest sampling note).**
- **Migrate ALL still-relevant definitions** that lived in the deleted Eq 1/Eq 2 "where" paragraph so the new primary equation is **fully self-defined** — nothing orphaned. Define every symbol in the new Eq (1): `y_i` = **cluster-level outcome for cluster i** (the modeled quantity); `α` baseline; `β` intervention effect; `x_i ∈ {0,1}` cluster treatment indicator; `ψ` scalar spillover parameter; `z_i ∈ {0,1}` cluster spillover indicator; `\mathcal{N}(i)` = set of clusters sharing a Rook edge with cluster i; `w_{ij}` = elements of the row-normalized weights matrix `W` (= 1/|𝒩(i)| for Rook neighbors, 0 otherwise; full construction Appendix A); `ρ` spatial autocorrelation coefficient; `ε_i` error term. Cross-check against the old paragraphs that nothing defined-once-and-still-used is lost.
- **Honest individual-level sentence (author wording, 6/21 — do NOT say "identical outcomes"):** outcomes are *modeled at the individual level*, but within any given cluster the specified response and covariates (`x_i`, `z_i`) are the same, so the **only** within-cluster variation is through an **independent (iid) random error**; because there is no error correlation within clusters, the outcome is presented at the cluster level (`y_i`) for clarity.
- Remove `k`, `k'` from the equation/notation; the honesty note is prose-level (no `k` subscript needed).
- Keep N(i) **defined and now actually used** (resolves the markup note "N(i) isn't referenced anywhere").
- Begin with `\noindent` ("no indent" markup, p12) — it continues the sentence after the equation, so no new-paragraph indent (see N4.3).

**N2.5 — Reconcile downstream prose (R1.2, C1, C2-in-§3.2) with cluster-level `y_i`.**
- The R1.2 ρ-vs-ψ paragraph (line 295) and the C1 estimands paragraph (lines 298, 302) reference `y_{ik}`, "individuals", "within that cluster", "individual responses differ only through ε_{ik}". Switch notation to cluster-level `y_i`/`ε_i`, but **preserve the substance** — the C1 point that individual responses differ only through the error term is true and should be retained, reworded cleanly.
- The C1 logic ("response on both sides → reduced form") **still holds** because the matrix form (N2.3) has `ρWy` on the RHS. Keep Eq (4) `Equation:ReducedForm` and its surrounding argument intact except for notation consistency.
- The C2 sentences in §3.2 about within-cluster variation contributing to the error-variance estimate are **correct per the code → keep and align notation**; do not delete. They support the honest framing.

### §2.6 "Spatial Dependence"

**N2.6 — Eq (SpatialDepTerm) cluster-level + N(i).**
- `ρ\mathbf{Wy} = ρ\sum_{j \in \mathcal{N}(i)} w_{ij} y_j` (drop `y_{jk'}`→`y_j`; `j≠i`→`j∈𝒩(i)`).
- "where j indexes clusters … " paragraph: drop `k'` / "subject in cluster j"; begin with `\noindent` (markup p15 "no indent").

### §2.7 "Estimation for SAR Model"

**N2.7 — Consolidate the redundant estimation equation → ✅ Option A LOCKED (author, 6/21).**
- Eq (6) `Equation:SARestimation` (lines 323–325) is the **same equation as Eq 3/new Eq 1**, only term-ordered differently (spatial lag first vs. last). Once standardized it is character-identical → redundant.
- **Action:** **Delete the displayed Eq (6)**, but **keep §2.7 "Estimation for SAR Model" as a lean prose subsection** (do NOT delete the whole section — it carries the estimation method, not a repeat). Surviving prose states: MLE is used to estimate the cluster-level spatial lag model in Equation~(\ref{Equation:SARmodel}); full likelihood/log-likelihood + MLE derivation in Appendix~B; estimates obtained via `\texttt{lagsarlm()}` from `\textbf{spdep}`.
- **Repoint** the `\ref{Equation:SARestimation}` at line 354 → `Equation:SARmodel`. Reads correctly: the structural Eq (1) is exactly what `lagsarlm()` fits, while the reduced form (`Equation:ReducedForm`) generates the data — the §3.2 distinction stays clean.
- **PRESERVE ALL METHOD CITATIONS (author emphasis, 6/21)** — even though we focus on the specific implementation:
  - §2.5 retained model intro keeps the **SAR-model** cites: `\cite{waller_applied_2004}`, `\cite{fotheringham_sage_2009}`, `\cite{ullah_handbook_1998}`.
  - §2.7 retained estimation prose keeps the **MLE** cites: `\cite{ord_estimation_1975}`, `\cite{bivand_comparing_2015}`, `\cite{bivand_computing_2013}`, and the **tool** cite `\cite{bivand_r_2022}` (lagsarlm/spdep).
  - Net: no `\cite` key may be dropped by this consolidation. Verify bibliography is unchanged (Phase 6 / N6.5).
- The "where all terms are as defined…" line is removed along with the equation; remaining §2.7 prose begins as a normal indented paragraph (first-after-heading, so `indentfirst` applies).

---

## Phase 3 — ICC clarification + limitation/future work

**N3.1 — ICC statement at the model (§2.5).** In the new "where" paragraph (N2.4), add one sentence: because the within-cluster errors are independent (no cluster-specific random effect), the specification does not parameterize intra-cluster correlation (ICC); the cluster-level outcome `y_i` captures the structural/design signal while individual variation is treated as independent noise.

**N3.2 — Limitation (§5.2).** Add a short sentence: the current specification models cluster-level outcomes and therefore does not estimate or account for ICC; incorporating individual-level outcomes with within-cluster correlation is left to future work.

**N3.3 — Future Considerations (§5.3).** One sentence flagging an individual-level / ICC-aware extension as a future direction (can pair with the existing buffer-zone / individual-location future-work text).

**N3.4 — Response letter:** mirror this (Phase 8).

---

## Phase 4 — Punctuation, indentation & formatting nits

**N4.1 — Commas/periods after equations.** Where an equation is mid-sentence (followed by "where …"), end the equation line with a comma; where it ends a sentence, a period. Apply to new Eq (1), the matrix form, SpatialDepTerm, and any appendix equations that are mid-sentence. (`Equation:ReducedForm` already has its trailing comma.)

**N4.2 — `\noindent` on "where" continuations.** Ensure every post-equation "where …" line that continues the sentence starts with `\noindent` (markup flagged: §2.5 line 292, §2.6 line 316, §2.7 line 326; line 283's `\noindent` is removed with old Eq 1).

**N4.3 — Indent the first paragraph of every section/subsection (author-confirmed).** Current behavior: first paragraphs after a heading are *not* indented while subsequent paragraphs are — the standard LaTeX default the author wants reversed. Fix document-wide with `\usepackage{indentfirst}` in the preamble (one consistent fix, not scattered `\indent`). Markup p12 "indent" at §2.5 confirms.
- **Critical distinction (author, 6/21):** an indent marks a *new paragraph*. When the text breaks only for a **displayed equation but the sentence continues** (the "where …" continuation), that is **not** a new paragraph → it must **not** indent. Enforce via `\noindent` on every such continuation line (N4.2). So: `indentfirst` handles genuine new paragraphs (incl. first-after-heading); `\noindent` suppresses the false indent after mid-sentence equations. These two rules together must be consistent everywhere.

**N4.4 — Indentation consistency sweep.** After N4.2/N4.3, scan §2.5–§2.7 and the rest of the document for: (a) every first-paragraph-after-heading now indents; (b) every post-equation "where/with" continuation uses `\noindent`; (c) no genuine new paragraph is missing its indent.

**N4.5 — Em-dashes.** Already clean — the 19 remaining `---` are all in `%` comments / `% ===` dividers, none in prose. **Verify only** (grep prose for `---` and unicode `—`); no edit expected.

---

## Phase 5 — Appendix equation numbering

**N5.1 — Unnumber appendix equations not referenced in the main text.** Appendix A (Eqs for `w̃_{ij}`, row-normalization) and Appendix B (cluster-level model, residuals, likelihood, log-likelihood) are **not** `\ref`'d anywhere (confirmed: only `Equation:SARmodel`, `Equation:ReducedForm`, `Equation:SARestimation` are referenced). Convert these appendix `equation`/`align` environments to starred `equation*`/`align*` (notes: "Add `*` to LaTeX equation calls"). Markup checkmarks confirm Appendix A and B content is otherwise correct.

**N5.2 — `j∈N(i)` in appendix sums.** Row-normalization denominators in Appendix A (markup p43/p44 underline on `\sum_{j≠i}`) → `\sum_{j \in \mathcal{N}(i)}` for consistency. Apply to Appendix B sums too (`\sum_{j≠i} w_{ij} y_j`).

**N5.3 — Appendix B term order.** Reorder the cluster-level scalar model in Appendix B (`y_i = ρΣ… + α + βx_i + ψz_i + ε_i`) to α-first, spatial-last to match new Eq (1). The note says vector forms stay as-is; only individual/scalar forms are reordered.

---

## Phase 6 — Compile & per-item verification

**N6.1 — Compile** `Revisions_V2c.tex` with the full sequence **including bibtex** (author emphasis): `pdflatex → bibtex → pdflatex → pdflatex`. `references.bib` must be present in the V2c folder so all in-text `\cite{}` references and the end bibliography render (a single pdflatex pass leaves `[?]` citations and drops the ~5-page bibliography). Confirm 0 errors, **0 undefined citations/references**, and that the bibliography is present in the PDF (the two benign `\float@end` warnings are expected). Page count noted but not a constraint this round.

**N6.2 — Grep gates:**
- No `y_{ik}`, `y_{jk'}`, `k'`, `\varepsilon_{ik}`, `x_{ik}`, `z_{ik}` remaining in the **main text** model sections (cluster-level only). (Appendix C subject-sampling step legitimately keeps subject language — confirm scope.)
- `\sum_{j \neq i}` eliminated in favor of `\sum_{j \in \mathcal{N}(i)}` everywhere.
- "causal" still appears only in the one permitted §5.2 scope note (preserve the 2b invariant).
- No `(Rx.x)` reviewer tags in body text.

**N6.3 — Cross-reference integrity:** every `\ref{Equation:*}` resolves (especially after the SARestimation consolidation). Verify the equation numbers cited in prose (reduced form, estimation) match the recompiled numbering.

**N6.4 — Citation integrity (author emphasis):** confirm the consolidation dropped **no** `\cite` keys. Diff the set of cited keys before/after; the SAR-model cites (waller/fotheringham/ullah) and MLE+tool cites (ord/bivand_comparing_2015/bivand_computing_2013/bivand_r_2022) must all still appear and resolve. 0 bibtex warnings; bibliography entry count unchanged.

**N6.5 — Per-item diff (optional, 2b-style):** snapshot before/after for the notation pass into `Revision 2c/Diffs/` if the author wants the incremental record.

---

## Phase 7 — Annotated diff (Revisions_V1 → Revisions_V2c)

Regenerate the full-markup annotated manuscript using the documented latexdiff workflow:
- `export PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$PATH"`
- Baseline: `Revision 2c/Original & Annotated Manuscripts/.../Revisions_V1.tex`
- `latexdiff --flatten --config "PICTUREENV=(?:picture|DIFnomarkup|tikzpicture|tabular)" "$V1" Revisions_V2c.tex > Revisions_Annotated_V2c.tex`
- Apply the known fixes: float-option leak regex; hardcode the two deleted-text labels that render `??` (`Equation:DataGen`→6, `appendix:MoransI`→D).
- Compile ×3 (bibtex between). Confirm 0 undefined refs and both add/del markup present.
- Clean build artifacts (keep `.tex` + `.pdf` only).

---

## Phase 8 — Response letter V6

**Source of truth = the author-edited `response_letter_V5.docx` in `Revision 2b/revision response/`** (NOT the `v5_letter_build.js` script — the author hand-edited the docx after that script ran, so the script is stale). Start from the edited docx, **edit it in place** (use the docx skill / python-docx, which preserves the existing 3-column tables and formatting), and save as `response_letter_V6`. Render docx→pdf with LibreOffice headless.

**Preserve verbatim from the edited V5 (do NOT rebuild or alter):**
- The salutation / addressing exactly as the author wrote it — **do not re-introduce or change editor/reviewer names**; carry forward whatever the edited letter currently says.
- The B&W 3-column layout and all existing C1–C6 / editor response blocks.

**Content updates (additive — alignment with V2c edits):**
- **Notation:** add/extend a response item noting the model is now presented at the **cluster level** (`y_i`), redundant general/matrix equations consolidated, summations expressed over the Rook-neighbor set `𝒩(i)`, and term order standardized.
- **ICC:** new statement that the model targets the cluster-level outcome (individuals share the cluster mean, only iid error varies; no within-cluster error correlation) so ICC is not parameterized — now explicitly flagged as a limitation and future direction (point to §2.5, §5.2, §5.3).
- **Formatting:** brief note that equation punctuation/indentation were corrected for consistency.
- Update any section/equation/page references that shifted due to the equation restructuring.

**Verify:** letter's cited section/equation numbers match the recompiled V2c PDF; addressing/names unchanged from the edited V5.

---

## Phase 9 — Reviewer subagent review (V2c vs V2b) — REQUIRED

Dispatch a `superpowers:code-reviewer` (or general review) subagent to compare **Revisions_V2c against Revisions_V2b** and confirm:
1. Every meeting-note item is addressed (cluster-level `y_i`; k/k' dropped; `j∈𝒩(i)`; general Eqs 1–2 dropped; matrix form with Zψ added; estimation eq consolidated/reordered; ICC statement + limitation + future work; appendix equations unnumbered; appendix sums use 𝒩(i); punctuation/indent fixed).
2. No regressions to locked 2b content (C1–C6 scientific framing intact; "causal" only in §5.2; no reviewer tags in body; reduced-form argument still coherent).
3. Prose matches the Phase 1 sim-code finding (no contradictory within-cluster-variation claims).
4. V2c compiles clean; annotated diff has 0 `??` and both add/del markup; response letter references are consistent.

Subagent returns a READY / NOT-READY verdict with a checklist. Surface it to the author.

---

## Deliverables (end state of Revision 2c)

1. **Response letter** — `revision response/response_letter_V6.{docx,pdf}`
2. **Annotated manuscript** (V1→V2c, adds/deletes) — `Walther_SpatialCRT_LaTeX_Revisions_V2c/Revisions_Annotated_V2c.{tex,pdf}`
3. **Final manuscript folder** — `Walther_SpatialCRT_LaTeX_Revisions_V2c/` with `Revisions_V2c.{tex,pdf}`, `references.bib`, `figures/`

## Decision points — all resolved (6/21)

- **D1 (Phase 0):** ✅ Archive 2b artifacts to `Revision 2c/_archive_2b/`; keep the `Revision 2b/` folder as canonical reference.
- **D2 (Phase 1):** ✅ Sim code is individual/point-level (iid per-individual error). Compromise: present cluster-level `y_i` for clarity + honest "modeled at individual level, only iid error varies within a cluster, no within-cluster error correlation" note; ICC framed as not-parameterized. C2 within-cluster-variation prose retained.
- **D3 (Phase 2, N2.7 — estimation equation):** ✅ Option A LOCKED. Delete the redundant displayed Eq (6); keep §2.7 as a lean prose subsection (MLE + Appendix B pointer + lagsarlm/spdep); repoint `\ref{Equation:SARestimation}` → `Equation:SARmodel`. Preserve ALL method citations (SAR-model + MLE + tool); drop no `\cite` keys.
- **D4 (Phase 4, N4.3 — indentation):** ✅ `indentfirst` for first-paragraph-after-heading indent, AND `\noindent` on mid-sentence post-equation continuations (no indent when it is not a new paragraph).
- **Model equation form:** ✅ Option (a) — displayed equation purely cluster-level (`ε_i`); individual-level detail carried in prose.

## Verification checklist (Karpathy pre-completion)

- [ ] Sim-code gate finding obtained before prose rewording (Phase 1)
- [ ] All main-text model notation cluster-level; grep gates pass (N6.2)
- [ ] Cross-references resolve; equation numbers in prose/letter match PDF (N6.3, Phase 8)
- [ ] No `\cite` keys dropped; bibliography unchanged (N6.4)
- [ ] No regression to locked C1–C6 framing (Phase 9)
- [ ] Three deliverables produced and named V2c/V6
- [ ] Reviewer subagent returns READY verdict

---

## Session Kickoff Prompt (paste into a fresh session)

The full operating detail now lives in this plan (§2a quick reference + phases 0→9), so the kickoff can be short:

```
Implement the Revision 2c manuscript updates for the Spatial CRT paper. Work in FAST mode.

Read and execute "Revision 2c/Notes/SpatialCRT_manuscript_updates_plan_6_21_26.md" end to end — it is the authoritative spec. Start at §2a (Execution quick reference) for the read-first list and the locked do-not-re-ask decisions, then work phases 0→9 in order.

All decisions are locked in the plan (versioning V2c; short track dropped/archived; cluster-level y_i with honest individual-level/iid-error wording; Option (a) equation; D3 estimation-equation consolidation; indentfirst + \noindent; response letter V6 built by editing the author-edited V5 docx in place). Do not re-ask them.

Hard requirements to honor exactly as written in the plan:
- Compile with bibtex (pdflatex → bibtex → pdflatex → pdflatex) so all \cite{} + the bibliography render; drop NO \cite keys vs. V2b.
- Run the Phase 6 grep + citation gates.
- FINAL: dispatch a superpowers:code-reviewer subagent comparing Revisions_V2c vs Revisions_V2b and return its READY/NOT-READY verdict against the Phase 9 checklist.

End by updating the status memory and giving a PIR summary.
```

> The long-form version of this prompt (with every equation string and `\cite` key inline) is reconstructable from §2a + phases 2/6/8 if a future session prefers a fully self-contained prompt over a pointer.
