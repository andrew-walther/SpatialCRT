#!/usr/bin/env bash
# ============================================================
# Script: sync_ctj_to_bios.sh
# Purpose: Copy the CTJ manuscript source and its two-column reading PDF into
#          bios-dissertation's manuscripts/project2-incidence/ as a reference
#          copy, and commit it there.
# Author: Andrew Walther
# Created: 2026-10-13
# Dependencies: bash, git
# ============================================================
#
# Usage:  tools/sync_ctj_to_bios.sh
#
# Normally run by the SpatialCRT post-commit hook (tools/post_commit_hook.sh)
# when a commit touches paper/ctj_manuscript/manuscript_ctj/; safe to run by
# hand at any time.
#
# Copies the COMMITTED versions (git show HEAD:<path>), not the working tree,
# so the copy always matches a SpatialCRT commit. Nothing is rendered: the PDF
# is whatever was last committed. If the manuscript source was committed more
# recently than CTJ_Reading.pdf, the PDF may be stale; the copy still goes
# ahead and the script prints a WARNING. (A comment-only .tex edit also
# triggers it; the check compares commit times, not content.)
#
#   CTJ_Manuscript.tex -> project_2_paper.tex
#   CTJ_Exhibits.tex   -> project_2_paper_exhibits.tex
#   CTJ_Reading.pdf    -> project_2_paper.pdf
#
# The copy does not build in bios-dissertation (sagej.cls, SageV.bst and the
# bib are symlinks in SpatialCRT, and \input names keep the CTJ_ file names).
# The submission-layout PDF (CTJ_Manuscript.pdf) is not copied; submissions go
# from SpatialCRT.
#
# Commits only those three paths, as "Sync CTJ manuscript from SpatialCRT
# <short-sha>", and only if something changed. Never pushes.
#
# Environment:
#   BIOS_DISSERTATION  bios-dissertation checkout (default
#                      ~/GithubProjects/bios-dissertation); override for tests.
# ============================================================

set -euo pipefail

fail() { echo "sync_ctj_to_bios: FAILED: $*" >&2; exit 1; }

# Paths -------------------------------------------------------------------
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPATIAL_ROOT="$(git -C "$TOOLS_DIR" rev-parse --show-toplevel)"
SRC_REL="projects/IncidenceDesign/paper/ctj_manuscript/manuscript_ctj"
SHA="$(git -C "$SPATIAL_ROOT" rev-parse --short HEAD)"

# Inside a git hook, GIT_DIR / GIT_INDEX_FILE etc. point at SpatialCRT; clear
# them so every git call below acts only on the repo named with -C.
unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_PREFIX GIT_OBJECT_DIRECTORY

BIOS="${BIOS_DISSERTATION:-$HOME/GithubProjects/bios-dissertation}"
TARGET_REL="manuscripts/project2-incidence"

[ -d "$BIOS/.git" ] || [ -f "$BIOS/.git" ] || fail "not a git checkout: $BIOS"
[ -d "$BIOS/$TARGET_REL" ] || fail "target directory missing: $BIOS/$TARGET_REL"

# source name : destination name (in TARGET_REL)
MAP=(
  "CTJ_Manuscript.tex:project_2_paper.tex"
  "CTJ_Exhibits.tex:project_2_paper_exhibits.tex"
  "CTJ_Reading.pdf:project_2_paper.pdf"
)

# 1. Stale-PDF check ----------------------------------------------------------
# Compare the last commit times of the source files and of the reading PDF.
src_time="$(git -C "$SPATIAL_ROOT" log -1 --format=%ct -- \
  "$SRC_REL/CTJ_Manuscript.tex" "$SRC_REL/CTJ_Exhibits.tex" "$SRC_REL/CTJ_Reading.tex")"
pdf_time="$(git -C "$SPATIAL_ROOT" log -1 --format=%ct -- "$SRC_REL/CTJ_Reading.pdf")"
STALE=0
if [ -n "$src_time" ] && [ -n "$pdf_time" ] && [ "$src_time" -gt "$pdf_time" ]; then
  STALE=1
  echo "sync_ctj_to_bios: WARNING: CTJ_Reading.pdf was last committed before the" >&2
  echo "  manuscript source, so the copied PDF may not show the latest text. Rebuild" >&2
  echo "  with 'latexmk CTJ_Reading.tex' in $SRC_REL and commit the PDF." >&2
fi

# 2. Copy the committed versions ----------------------------------------------
PATHS=()
for pair in "${MAP[@]}"; do
  src="${pair%%:*}"; dst="$TARGET_REL/${pair#*:}"
  git -C "$SPATIAL_ROOT" cat-file -e "HEAD:$SRC_REL/$src" 2>/dev/null ||
    fail "$SRC_REL/$src is not committed in SpatialCRT"
  git -C "$SPATIAL_ROOT" show "HEAD:$SRC_REL/$src" >"$BIOS/$dst"
  PATHS+=("$dst")
done

# 3. Commit only those paths ----------------------------------------------------
git -C "$BIOS" add -- "${PATHS[@]}"
BIOS_COMMIT="no commit (nothing changed)"
if ! git -C "$BIOS" diff --cached --quiet -- "${PATHS[@]}"; then
  git -C "$BIOS" commit -q -m "Sync CTJ manuscript from SpatialCRT $SHA" -- "${PATHS[@]}"
  BIOS_COMMIT="bios-dissertation $(git -C "$BIOS" rev-parse --short HEAD)"
fi

# 4. Summary ------------------------------------------------------------------
NOTE=""
[ "$STALE" -eq 0 ] || NOTE="; WARNING: reading PDF older than source"
echo "sync_ctj_to_bios: OK: CTJ manuscript from SpatialCRT $SHA -> $BIOS_COMMIT$NOTE"
