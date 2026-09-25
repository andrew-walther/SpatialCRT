#!/bin/sh
# ============================================================
# Script: post_commit_hook.sh
# Purpose: SpatialCRT post-commit logic. If HEAD touched the Chapter 3
#          sources, run sync_to_prelim.sh. Never blocks, fails, or undoes the
#          SpatialCRT commit (the commit already exists when this runs).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-25
# Dependencies: POSIX sh, git
# ============================================================
# Installed as a thin shim by tools/install_hooks.sh; the shim calls this file
# from the current worktree, so the logic stays versioned here.
#
# Trigger paths (relative to the SpatialCRT root):
#   projects/IncidenceDesign/paper/dissertation_chapter/Dissertation_Chapter.qmd
#   projects/IncidenceDesign/paper/dissertation_chapter/figures/
#   projects/IncidenceDesign/paper/dissertation_chapter/tools/
#   projects/IncidenceDesign/paper/SpatialCRT_IncidenceDesign.bib
# Output is appended to tools/sync.log (gitignored).

root=$(git rev-parse --show-toplevel 2>/dev/null) || exit 0
tools="$root/projects/IncidenceDesign/paper/dissertation_chapter/tools"
log="$tools/sync.log"
sync="$tools/sync_to_prelim.sh"

# Files changed by HEAD (first parent for merges; --root for a first commit).
changed=$(git diff-tree --no-commit-id --name-only -r --root -m --first-parent HEAD 2>/dev/null)
echo "$changed" | grep -q -E '^projects/IncidenceDesign/paper/(dissertation_chapter/(Dissertation_Chapter\.qmd$|figures/|tools/)|SpatialCRT_IncidenceDesign\.bib$)' || exit 0

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
  echo "post-commit: rebase in progress; Chapter 3 sync skipped. Run $sync afterwards." >&2
  exit 0
fi

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
