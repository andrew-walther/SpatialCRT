# Overnight plan, 2026-10-21

**Status: approved by the author 2026-10-21 (about 03:45).** It follows the decisions
recorded in `DECISIONS.md` (2026-10-21 entry), made while reviewing
[../findings/2026-10-20-methods-extensions-findings.md](../findings/2026-10-20-methods-extensions-findings.md).
No manuscript, chapter, CTJ, SI or slide file is edited by this plan.

## What runs

| Step | Work | Output | Time |
|---|---|---|---|
| 0 | Grid burden share on the stored grid surfaces (`code/22_grid_burden_share.R`) | `results/grid_burden_share/` | done before launch (1 min) |
| 1 | v2 checks (`application/tests/test_methods_ext_v2.R`) | `verification/test_methods_ext_v2.txt` | 2 min |
| 2 | A3 rerun, profile `a3v2`: saturation regions of at least 2 colleges, and the t + df-corrected test's rejections and coverage stored per fit cell | `participation/a3v2/` | about 4 h |
| 3 | Reproduction check against the 2026-10-20 A3 run (`application/tests/check_a3v2_reproduction.R`) | `a3v2_reproduction.log` | minutes |
| 4 | Bayesian CAR at ρ = 0.5 with chains four times longer, on the same data (`a6bayeslong`) | `a6bayeslong/` | about 4 h |
| 5 | Summaries: A3 tables with t + df power beside size-adjusted and Wald (`summaries_v2/`); long-chain Bayes comparison (`sensitivity_summaries_v2/`) | as named | minutes |
| 6 | V9: production computation files unchanged | `v9_v2.log` | seconds |

Paths are under `application/results/methods_ext_2026-10/` unless they start with
`results/`. The queue is `application/code/run_methods_ext_queue_v2.sh`, and its
status log is `queue_v2.log`.

## Checks and stop rules

- Any step that exits non-zero, or that does not print its all-tasks-succeeded line,
  stops the queue with `QUEUEFAILED <step>`. Nothing after it runs.
- The reproduction check must pass. The rerun uses the same seeds and code paths as
  the 2026-10-20 run, so the Wald and size-adjusted counts must be bit-identical for
  every design. The exceptions are saturation subsets whose regions changed under the
  new rule: those subsets, and the size-adjusted counts of their tasks, are excluded.
- The t + df counts must reject no more often, and cover no less often, than Wald.
- The long-chain Bayes summary stops unless the SAR fits on the shared data are
  identical in both runs.
- V9 must pass.

Smoke checks before launch (2026-10-21): the v2 checks passed; `pilotv2` reproduced
`pilot` exactly; the summaries and headline tables ran with the t + df columns; and
the long-chain Bayes pipeline ran on the smoke profile.

## After the queue

If it ends with `QUEUEDONE`, the session that launched it:

1. checks the logs and the key numbers;
2. adds the 2026-10-21 results to the findings doc;
3. writes `../findings/2026-10-21-methods-extensions-report.md`, covering the findings, proposed
   text for the thesis appendix, the CTJ SI and the few main-text sentences (every
   number traced to its CSV), and proposed Part II slides;
4. adds those slides to `bios-dissertation/slides/prelim/OUTLINE.md`.

If it ends with `QUEUEFAILED`, the session records what failed and why in the
findings doc and drafts nothing.

Commits stay local, with no push. The author decides on pushing in the morning.
