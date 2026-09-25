#!/usr/bin/env bash
# ============================================================
# Script: install_hooks.sh
# Purpose: Install SpatialCRT's post-commit hook, which keeps
#          bios-dissertation's Chapter 3 in sync (tools/post_commit_hook.sh).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-25
# Dependencies: bash, git
# ============================================================
# The hook goes in the repo's COMMON hooks dir (git rev-parse --git-common-dir),
# so every worktree shares it. The installed file is a thin shim that runs
# tools/post_commit_hook.sh from whichever worktree made the commit, so the
# hook's logic stays versioned. Re-running is safe: it replaces only a hook it
# installed itself (recognised by the marker line) and refuses to overwrite
# anyone else's post-commit hook.

set -euo pipefail

MARKER="# spatialcrt-chapter3-prelim-sync"
cd "$(dirname "${BASH_SOURCE[0]}")"

if hp="$(git config --get core.hooksPath)"; then
  echo "install_hooks: core.hooksPath is set ($hp); install the hook there by hand." >&2
  exit 1
fi

common="$(git rev-parse --path-format=absolute --git-common-dir)"
hooks="$common/hooks"
hook="$hooks/post-commit"
mkdir -p "$hooks"

if [ -e "$hook" ] && ! grep -q "^$MARKER\$" "$hook"; then
  echo "install_hooks: $hook already exists and was not installed by this script." >&2
  echo "  Not overwriting it. Add this line to it by hand instead:" >&2
  echo "    sh \"\$(git rev-parse --show-toplevel)/projects/IncidenceDesign/paper/dissertation_chapter/tools/post_commit_hook.sh\"" >&2
  exit 1
fi

cat >"$hook" <<EOF
#!/bin/sh
$MARKER
# Installed by projects/IncidenceDesign/paper/dissertation_chapter/tools/install_hooks.sh.
# Runs the versioned hook logic from the worktree that made the commit; it never
# blocks or undoes the commit.
h="\$(git rev-parse --show-toplevel)/projects/IncidenceDesign/paper/dissertation_chapter/tools/post_commit_hook.sh"
[ -f "\$h" ] && sh "\$h"
exit 0
EOF
chmod +x "$hook"
echo "install_hooks: installed $hook"
