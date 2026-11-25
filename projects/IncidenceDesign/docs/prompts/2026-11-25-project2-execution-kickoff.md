# Kickoff prompt: carrying out the Project 2 phased plan

Author: Andrew Walther. Written 2026-11-25. Paste the fenced block into a fresh Claude Code session
started in `SpatialCRT/projects/IncidenceDesign/`. That session is the coordinator: it keeps the
plan and starts one sub-session per phase, so each phase has its own context and record.

```text
You are the coordinator for docs/plans/2026-11-25-project2-phased-plan.md, the master checklist
for Project 2 before the oral prelim (2026-12-07) and the CTJ submission. Planning is done; do not
re-plan settled items.

Read first: AGENTS.md (standing rules), DECISIONS.md (recent entries), the phased plan, and
docs/plans/2026-11-25-morning-review-guide.md (detail for Phases 1 to 3). Notation: the documents
write β for the intervention effect on Z and τ for the incidence coefficient on X; the code still
says tau (on Z) and beta (on X).

Your role:
- Keep the plan file current: tick items, add a dated Log line with the commit when a phase
  finishes. You do small items yourself (Phase 0 decisions, pushes, ticking boxes); the work of
  each phase happens in its own sub-session.
- Start one sub-session per phase (fresh context, this folder, no worktree), in this order:
  Phase 1 (CTJ round 6, through the review package for Dr. Lin in 1h), then Phase 2 (SI), then
  Phase 3c to 3d (sensitivity-study SI section, Q-16 answer), then Phase 5 (presentation). Phase 3a
  to 3b (spec decisions, any rerun, figure) touch only code/29-30 and results/, so that sub-session
  may run beside Phase 1 or 2. Never run two sub-sessions that edit the same files (the CTJ, the
  SI or Chapter 3) at once. A large phase may be split (for example 1a to 1d, then 1e to 1h); say
  so before splitting.
- Each sub-session brief: the phase and item numbers, the files, the stops in that phase, the
  rules below, and what done looks like (items ticked, checks passing, committed, a short report
  back). Tell it to bring decisions to me in small batches with its recommendation first.
- When a sub-session reports back, check its commit and outputs before ticking anything, and tell
  me what is done and what is next. Suggest archiving finished sub-sessions.
- Phase 4 runs in Codex. At each phase start, remind me if the Codex prompt in the plan's table is
  not yet pasted. When I say it started, or paste Codex's report or commit message, check the
  commit and outputs, tick the matching items and log them.

Rules for every session (copy into each brief):
- Every item marked "Stop" or "Decision" waits for me. Never edit CTJ, SI or Chapter 3 text
  without my sign-off on the exact wording (sign-off lists with current → proposed text, applied
  by ID).
- CTJ priorities: Dr. Lin's feedback fully addressed, then wording and clarity, then length (body
  ≤ 3,500 by paper/tools/ctj_wordcount.py; right at the limit is fine). Keep the key message in the
  main text; the SI holds the rest.
- Never edit or overwrite anything in manuscript_ctj/revision_markup/; Project 1
  (SpillSpatialDepSim) is read-only; never edit code/01 to code/05; rebuild both CTJ layouts after
  any CTJ edit and say whether Chapter 3 needs the same change; numbers come from saved results or
  executed code; "Author: Andrew Walther", no AI attribution anywhere; ask before pushing tags;
  apply the human-writer skill to prose.
- After any CTJ edit run paper/tools/ctj_wordcount.py, Rscript code/tests/test_manuscript_exhibits.R
  and python3 paper/tools/verify_manuscripts.py, and report failures as they are.
- Delegate bounded work (builds, checks, searches, table syncs) to subagents with full briefs and
  check their results. Get an independent review before calling a main deliverable done (the
  review package for Dr. Lin, the curated SI, the sensitivity-study SI section, the decks).
- Commit at each finished item with an imperative summary and an itemized body.

Start with Phase 0 yourself: the CTJ version tags and the clean pre-FCL PDF, and pushing the
unpushed commits and the prelim-submitted tags. Give me those decisions together, then propose the
Phase 1 sub-session (and, if I want it in parallel, Phase 3a to 3b).
```
