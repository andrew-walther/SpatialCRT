#!/usr/bin/env bash
# ============================================================
# Script: sync_to_prelim.sh
# Purpose: Regenerate bios-dissertation's Chapter 3 draft from
#          Dissertation_Chapter.qmd (the ONLY file anyone edits), render it,
#          and commit the generated copy in bios-dissertation.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-25
# Dependencies: bash, git, python3 (stdlib only), RStudio's bundled Quarto,
#               TinyTeX
# ============================================================
#
# Usage:  tools/sync_to_prelim.sh [--no-render]
#
# Normally run by the SpatialCRT post-commit hook (tools/install_hooks.sh);
# safe to run by hand at any time.
#
# Steps:
#   1. Bib check + transform (prelim_transform.py). On any citekey or bib-field
#      mismatch it exits non-zero BEFORE writing anything.
#   2. Write draft/project2-incidence-draft.qmd, sync draft/figures/, copy
#      manuscript-declarations.md, make sure the bios-prelim.cls symlink exists.
#   3. Render the draft with RStudio's Quarto (skipped with --no-render).
#   4. Stage ONLY the Chapter 3 paths under
#      prelim/project-proposals/project2-incidence/, explicitly by path, and
#      commit them with "Sync Chapter 3 from SpatialCRT <short-sha>" (only if
#      something changed). `git commit -- <paths>` commits those paths alone,
#      even if something else happens to be staged. Never pushes.
#   5. Print a one-line summary.
#
# Environment:
#   BIOS_DISSERTATION  bios-dissertation checkout (default
#                      ~/GithubProjects/bios-dissertation); override for tests.
# ============================================================

set -euo pipefail

NO_RENDER=0
for arg in "$@"; do
  case "$arg" in
    --no-render) NO_RENDER=1 ;;
    -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
    *) echo "sync_to_prelim: unknown argument: $arg" >&2; exit 64 ;;
  esac
done

fail() { echo "sync_to_prelim: FAILED: $*" >&2; exit 1; }

# Paths -------------------------------------------------------------------
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHAPTER_DIR="$(dirname "$TOOLS_DIR")"
PAPER_DIR="$(dirname "$CHAPTER_DIR")"
CHAPTER_QMD="$CHAPTER_DIR/Dissertation_Chapter.qmd"
CHAPTER_BIB="$PAPER_DIR/SpatialCRT_IncidenceDesign.bib"
DECLARATIONS="$CHAPTER_DIR/manuscript-declarations.md"
SPATIAL_ROOT="$(git -C "$CHAPTER_DIR" rev-parse --show-toplevel)"
CHAPTER_REL="${CHAPTER_QMD#"$SPATIAL_ROOT"/}"
SHA="$(git -C "$SPATIAL_ROOT" rev-parse --short HEAD)"

# Inside a git hook, GIT_DIR / GIT_INDEX_FILE etc. point at SpatialCRT; clear
# them so every git call below acts only on the repo named with -C.
unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_PREFIX GIT_OBJECT_DIRECTORY

BIOS="${BIOS_DISSERTATION:-$HOME/GithubProjects/bios-dissertation}"
TARGET_REL="prelim/project-proposals/project2-incidence"
TARGET="$BIOS/$TARGET_REL"
DRAFT_DIR="$TARGET/draft"
MASTER_BIB="$BIOS/prelim/references.bib"
QUARTO="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"

[ -f "$CHAPTER_QMD" ] || fail "chapter not found: $CHAPTER_QMD"
[ -d "$BIOS/.git" ] || [ -f "$BIOS/.git" ] || fail "not a git checkout: $BIOS"
[ -d "$TARGET" ] || fail "target directory missing: $TARGET"
[ -f "$MASTER_BIB" ] || fail "master bib missing: $MASTER_BIB"

dirty="$(git -C "$SPATIAL_ROOT" status --porcelain -- "$CHAPTER_QMD" "$CHAPTER_DIR/figures" \
         "$TOOLS_DIR" "$CHAPTER_BIB" "$DECLARATIONS" | grep -v 'tools/sync\.log$' || true)"
if [ -n "$dirty" ]; then
  echo "sync_to_prelim: WARNING: SpatialCRT inputs have uncommitted changes; the copy" >&2
  echo "  is generated from the working tree, not exactly from $SHA:" >&2
  echo "$dirty" | sed 's/^/    /' >&2
fi

# 1-2. Bib check, transform, figures ---------------------------------------
RESULT="$(python3 "$TOOLS_DIR/prelim_transform.py" \
  --chapter "$CHAPTER_QMD" --chapter-bib "$CHAPTER_BIB" --master-bib "$MASTER_BIB" \
  --header "$TOOLS_DIR/prelim_header.yml" --out-dir "$DRAFT_DIR" \
  --source-label "SpatialCRT/$CHAPTER_REL")" || fail "bib check / transform (see message above); nothing written"

cp "$DECLARATIONS" "$TARGET/manuscript-declarations.md"
if [ ! -L "$DRAFT_DIR/bios-prelim.cls" ]; then
  ln -s ../../../../templates/skeleton/prelim/bios-prelim.cls "$DRAFT_DIR/bios-prelim.cls"
fi
[ -f "$DRAFT_DIR/bios-prelim.cls" ] || fail "bios-prelim.cls symlink does not resolve"

# 3. Render ----------------------------------------------------------------
PAGES="not rendered"
if [ "$NO_RENDER" -eq 0 ]; then
  [ -x "$QUARTO" ] || fail "RStudio Quarto not found at $QUARTO"
  TINYTEX_BIN="$(ls -d "$HOME"/Library/TinyTeX/bin/* 2>/dev/null | head -1)"
  [ -n "$TINYTEX_BIN" ] || fail "TinyTeX not found under ~/Library/TinyTeX/bin"
  export PATH="$TINYTEX_BIN:$PATH"
  RENDER_LOG="$(mktemp -t sync_to_prelim_render)"
  if ! (cd "$DRAFT_DIR" && "$QUARTO" render project2-incidence-draft.qmd) >"$RENDER_LOG" 2>&1; then
    tail -40 "$RENDER_LOG" >&2
    fail "quarto render failed (full log: $RENDER_LOG); nothing committed"
  fi
  # Unresolved citations / references exit 0, so check the output explicitly.
  if grep -i -E 'citeproc: reference .* not found|unable to resolve|undefined references|Reference .* undefined|Citation .* undefined' "$RENDER_LOG" >&2; then
    fail "render reported unresolved citations or references (log: $RENDER_LOG); nothing committed"
  fi
  PDF="$DRAFT_DIR/project2-incidence-draft.pdf"
  if command -v pdftotext >/dev/null 2>&1; then
    n_qq="$(pdftotext "$PDF" - | grep -c '??' || true)"
    [ "$n_qq" = "0" ] || fail "rendered PDF contains '??' $n_qq time(s) (an unresolved \\ref); nothing committed"
  fi
  if command -v pdfinfo >/dev/null 2>&1; then
    PAGES="$(pdfinfo "$PDF" | awk '/^Pages:/{print $2}') pp"
  else
    PAGES="rendered"
  fi
  rm -f "$RENDER_LOG"
fi

# 4. Commit the Chapter 3 paths only ----------------------------------------
CANDIDATES=(
  "$TARGET_REL/README.md"
  "$TARGET_REL/project2-incidence.qmd"
  "$TARGET_REL/manuscript-declarations.md"
  "$TARGET_REL/.gitignore"
  "$TARGET_REL/draft/project2-incidence-draft.qmd"
  "$TARGET_REL/draft/bios-prelim.cls"
  "$TARGET_REL/draft/figures"
)
PATHS=()
for p in "${CANDIDATES[@]}"; do
  if [ -e "$BIOS/$p" ] || [ -L "$BIOS/$p" ] || [ -n "$(git -C "$BIOS" ls-files -- "$p")" ]; then
    PATHS+=("$p")
  fi
done

git -C "$BIOS" add -A -- "${PATHS[@]}"
# `git add` silently skips ignored files (bios-dissertation ignores *.pdf), so
# a figure that is ignored would vanish from the commit without an error.
ignored="$(git -C "$BIOS" ls-files --others --ignored --exclude-standard -- "$TARGET_REL/draft/figures")"
[ -z "$ignored" ] || fail "figures are gitignored in bios-dissertation and would not be committed: $ignored"
BIOS_COMMIT="no commit (nothing changed)"
if ! git -C "$BIOS" diff --cached --quiet -- "${PATHS[@]}"; then
  git -C "$BIOS" commit -q -m "Sync Chapter 3 from SpatialCRT $SHA" -- "${PATHS[@]}"
  # Belt and braces: the commit must touch nothing outside the target dir.
  outside="$(git -C "$BIOS" show --name-only --format= HEAD | grep -v "^$TARGET_REL/" || true)"
  [ -z "$outside" ] || fail "commit $(git -C "$BIOS" rev-parse --short HEAD) touched paths outside $TARGET_REL: $outside"
  BIOS_COMMIT="bios-dissertation $(git -C "$BIOS" rev-parse --short HEAD)"
fi

# 5. Summary -----------------------------------------------------------------
# RESULT looks like: RESULT qmd=changed figs_copied=N figs_removed=M keys=K
echo "sync_to_prelim: OK: Chapter 3 from SpatialCRT $SHA -> $BIOS_COMMIT; ${RESULT#RESULT }; $PAGES"
