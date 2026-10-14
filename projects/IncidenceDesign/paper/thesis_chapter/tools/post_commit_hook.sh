#!/bin/sh
# ============================================================
# Script: post_commit_hook.sh
# Purpose: SpatialCRT post-commit logic. If HEAD touched the Chapter 3
#          sources, run sync_to_prelim.sh; if it touched the CTJ manuscript,
#          run sync_ctj_to_bios.sh. Never blocks, fails, or undoes the
#          SpatialCRT commit (the commit already exists when this runs).
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: POSIX sh, git
# ============================================================
# Installed as a thin shim by tools/install_hooks.sh; the shim calls this file
# from the current worktree, so the logic stays versioned here.
#
# Trigger paths (relative to the SpatialCRT root):
#   projects/IncidenceDesign/paper/thesis_chapter/Thesis_Chapter.qmd
#   projects/IncidenceDesign/paper/thesis_chapter/figures/
#   projects/IncidenceDesign/paper/thesis_chapter/tools/
#   projects/IncidenceDesign/paper/SpatialCRT_IncidenceDesign.bib
# CTJ trigger path (added 2026-10-13):
#   projects/IncidenceDesign/paper/ctj_manuscript/manuscript_ctj/
# Output is appended to tools/sync.log (gitignored).

root=$(git rev-parse --show-toplevel 2>/dev/null) || exit 0
tools="$root/projects/IncidenceDesign/paper/thesis_chapter/tools"
log="$tools/sync.log"
sync="$tools/sync_to_prelim.sh"
sync_ctj="$tools/sync_ctj_to_bios.sh"

# Files changed by HEAD (first parent for merges; --root for a first commit).
changed=$(git diff-tree --no-commit-id --name-only -r --root -m --first-parent HEAD 2>/dev/null)
ch3=0; ctj=0
echo "$changed" | grep -q -E '^projects/IncidenceDesign/paper/(thesis_chapter/(Thesis_Chapter\.qmd$|figures/|tools/)|SpatialCRT_IncidenceDesign\.bib$)' && ch3=1
echo "$changed" | grep -q -E '^projects/IncidenceDesign/paper/ctj_manuscript/manuscript_ctj/' && ctj=1
[ "$ch3" -eq 1 ] || [ "$ctj" -eq 1 ] || exit 0

sha=$(git rev-parse --short HEAD 2>/dev/null)

warn() {
  {
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!! Chapter 3 prelim sync FAILED for SpatialCRT commit $sha."
    echo "!! $1"
    echo "!! Your SpatialCRT commit is fine; bios-dissertation was NOT updated."
    echo "!! Fix the problem, then run by hand:"
    echo "!!   $sync"
    echo "!! Log: $log"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""
  } >&2
}

# Mid-rebase every replayed commit would trigger a sync; skip and say so.
if [ -d "$(git rev-parse --git-path rebase-merge)" ] || [ -d "$(git rev-parse --git-path rebase-apply)" ]; then
  echo "post-commit: rebase in progress; bios-dissertation syncs skipped. Run $sync and/or $sync_ctj afterwards." >&2
  exit 0
fi

# CTJ manuscript copy (fast: copies committed files, no render). Its result
# line, including any stale-PDF WARNING, is shown in the terminal.
if [ "$ctj" -eq 1 ]; then
  echo "post-commit: $sha touched the CTJ manuscript; copying to bios-dissertation (log: $log) ..." >&2
  {
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') SpatialCRT $sha (CTJ)"
    "$sync_ctj"
    echo "=== exit $?"
  } >>"$log" 2>&1
  if [ "$(tail -1 "$log")" = "=== exit 0" ]; then
    # Only this run's lines: from its last "=== ... (CTJ)" header onward.
    awk '/^=== .* \(CTJ\)$/{buf=""} {buf=buf $0 "\n"} END{printf "%s", buf}' "$log" |
      grep -E '^sync_ctj_to_bios: (OK|WARNING)' | sed 's/^/post-commit: /' >&2
  else
    {
      echo ""
      echo "!! CTJ manuscript copy FAILED for SpatialCRT commit $sha."
      tail -3 "$log" | sed 's/^/!!   /'
      echo "!! Your SpatialCRT commit is fine; bios-dissertation was NOT updated."
      echo "!! Fix the problem, then run by hand: $sync_ctj"
      echo ""
    } >&2
  fi
fi
[ "$ch3" -eq 1 ] || exit 0

if [ ! -x "$sync" ]; then
  warn "sync script not found or not executable: $sync"
  exit 0
fi

echo "post-commit: $sha touched Chapter 3 sources; syncing to bios-dissertation (log: $log) ..." >&2
{
  echo "=== $(date '+%Y-%m-%d %H:%M:%S') SpatialCRT $sha"
  "$sync"
  rc=$?
  echo "=== exit $rc"
} >>"$log" 2>&1

last=$(tail -2 "$log" | head -1)
case "$(tail -1 "$log")" in
  "=== exit 0") echo "post-commit: $last" >&2 ;;
  *)
    # Quote this run's error lines (from its "===" header onward), since the
    # sync's own stderr went to the log, not the terminal.
    reason=$(awk '/^=== .* SpatialCRT /{buf=""} {buf=buf $0 "\n"} END{printf "%s", buf}' "$log" |
             grep -E 'ERROR|FAILED|^  - ' | sed 's/^/!!   /')
    warn "Reason (from the log):
$reason" ;;
esac
exit 0
