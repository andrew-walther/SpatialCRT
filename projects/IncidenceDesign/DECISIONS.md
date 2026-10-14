# DECISIONS.md — IncidenceDesign

Author decisions with their reasons, newest first. `AGENTS.md` keeps only the
one-line rule each decision implies.

---

## 2026-10-14 — Best designs by scenario, not one recommended design

**Decision.** The chapter and manuscripts report the best designs by scenario,
mainly by spillover regime, instead of naming one design for all settings.

**Why.** The rankings change with the spillover regime. Checkerboard and 2×2
Blocking are near the best designs under both-arms spillover and the worst two
under control-only spillover, while Incidence-Guided Saturation stays first or
second under both. A trial team may not know the regime in advance, so an
answer given for each scenario is more robust than a single recommendation.

**Replaces.** The "one recommended design" in
`docs/plans/simulation-revision-plan.md` (left as a dated record) and in the
Phase D ROADMAP item.

**Files.** `paper/thesis_chapter/Thesis_Chapter.qmd` (Discussion, guidance by
spillover scenario, added 2026-10-11).

---

## 2026-10-02 — Application outcome: no SUD baseline term in the primary model

**Decision.** Primary continuous
education outcomes have no SUD-incidence baseline term (β = 0); observed rates
inform allocation. An explicitly hypothetical β = 1 baseline sensitivity is
accepted with rank-scaled X and an X-adjusted fit; primary ML fits treatment and
spillover without X. τ is constant within
a scenario in both cases. This supersedes older matched-primary β = 1 proposals;
the completed grid study remains unchanged.

**Source.** Author clarification, 2026-10-02 (moved from AGENTS.md 2026-10-14).

---

## 2026-10-02 — Application implementation approved

**Decision.** Implementation approved October 2. Continue the
approved plan without re-interviewing settled choices. The companion is
`application/report/real_sud_companion.html`; computation lives in
`application/code/real_sud_{setup,simulation}.R` and `run_real_sud.R`.
Follow `application/README.md` and the implementation findings for run status;
smoke/pilot outputs are verification stages, not full-study conclusions.
The full chapter/appendix, CTJ/SI and comprehensive review (2026-10-12) are done.
Old deadline dates no longer govern work.

Application design authority: start with
`docs/plans/nc_sud_application_plan_2026-10-02.md`; its status records settled
versus proposed choices. The integrated prompt is
`docs/plans/continuation_prompt_2026-10-02.md`. Do not treat proposals as approved
model/design changes.

**Source.** Moved from AGENTS.md "Standing Author Rules" 2026-10-14.

---

## 2026-10-02 — Checkerboard stays in the main comparison

**Decision.** Retain Checkerboard in
Project 2's main design comparison, with detailed identification diagnostics in
the appendix. A brief Chapter 2 link motivates the extension; no direct numerical
cross-project comparison is required. CTJ must cite the verified Project 1
publication to establish the extension while remaining self-contained.

**Source.** Author-confirmed 2026-10-02 (moved from AGENTS.md 2026-10-14).

---

## 2026-10-02 — CTJ author declarations

**Decision.** Retain Walther, Habib, Simpson and Lin in that order,
with Lin corresponding. Author confirmed no current funding or conflicts.
The in-text ethics placeholder was removed (author, 2026-10-12); the IRB statement
is a TODO before submission. Never invent a protocol number or claim approval/exemption.

**Source.** Moved from AGENTS.md 2026-10-14.

---

## 2026-09-27 — How CTJ and the chapter refer to Project 1

**Decision.** Dissertation prose may refer directly
to "Chapter 2". CTJ must use self-contained wording (e.g., "In previous work...")
with a citation to the Project 1 BMC Medical Research Methodology paper, never
"Chapter 2" or "Project 1" as reader-facing references. Describe findings within
their studied conditions; verify final bibliographic details before adding the citation.

**Note (2026-10-14).** The Project 1 paper is now published (doi:10.1186/s12874-026-02968-0).

**Source.** Author-confirmed 2026-09-27 (moved from AGENTS.md 2026-10-14).
