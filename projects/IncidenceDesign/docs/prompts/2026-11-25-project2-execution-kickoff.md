# Kickoff prompt: carrying out the Project 2 phased plan

Author: Andrew Walther. Written 2026-11-25. Paste the fenced block into a fresh Claude Code session
started in `SpatialCRT/projects/IncidenceDesign/`. The planning happened in an earlier session;
this one carries out the plan.

```text
We are carrying out docs/plans/2026-11-25-project2-phased-plan.md, the master checklist for
Project 2 before the oral prelim (2026-12-07) and the CTJ submission. Planning is done; do not
re-plan settled items.

Read first: AGENTS.md (standing rules), DECISIONS.md (recent entries), the phased plan, and
docs/plans/2026-11-25-morning-review-guide.md (detail for Phases 1 to 3). Notation: the documents
write β for the intervention effect on Z and τ for the incidence coefficient on X; the code still
says tau (on Z) and beta (on X).

How to work:
- One phase at a time, in order (Phase 0, then 1a to 1h, then 2, then 3; Phase 4 runs in Codex on
  the side; Phase 5 last). Check items off in the plan file as they finish and add a dated line
  with the commit to its Log.
- Every item marked "Stop" or "Decision" waits for me. Bring decisions to me in small batches with
  your recommendation first, not one long list. Never edit the CTJ, SI or Chapter 3 text without
  my sign-off on the exact wording (sign-off lists with current → proposed text, applied by ID).
- At the start of each phase, remind me if a Codex prompt in the plan's table is not yet pasted.
- Priorities for the CTJ: Dr. Lin's feedback fully addressed, then wording and clarity, then
  length (body ≤ 3,500 by paper/tools/ctj_wordcount.py; right at the limit is fine). Keep the key
  message in the main text.
- Rules that always hold: never edit or overwrite anything in manuscript_ctj/revision_markup/;
  Project 1 (SpillSpatialDepSim) is read-only; never edit code/01 to code/05; rebuild both CTJ
  layouts after any CTJ edit and say whether Chapter 3 needs the same change; numbers come from
  saved results or executed code; "Author: Andrew Walther", no AI attribution anywhere; ask before
  pushing tags; apply the human-writer skill to prose.
- After any CTJ edit run paper/tools/ctj_wordcount.py, Rscript code/tests/test_manuscript_exhibits.R
  and python3 paper/tools/verify_manuscripts.py, and report failures as they are.
- Delegate bounded work (builds, checks, searches, table syncs) to subagents with full briefs and
  check their results. Get an independent review before calling a phase's main deliverable done
  (the review package for Dr. Lin in 1h, the curated SI, the sensitivity-study SI section).
- Commit at each finished item with an imperative summary and an itemized body. Tell me when a
  phase is done and suggest a fresh session if this one is getting long.

Start with Phase 0: the CTJ version tags and the clean pre-FCL PDF, and pushing the six unpushed
commits and the prelim-submitted tags. Give me those decisions together, then move to Phase 1a.
```
