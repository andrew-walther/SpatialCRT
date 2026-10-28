#!/usr/bin/env bash
# ============================================================
# Script: refresh_figures.sh
# Purpose: Copy the current manuscript figures into presentation/figures/ as PNGs
#          for the talk deck, then draw the deck-only NC example map.
# Author: Andrew Walther
# Created: 2026-10-28
# Dependencies: pdftoppm (poppler), Rscript
# ============================================================
#
# Run from the project root:  bash presentation/tools/refresh_figures.sh
#
# The deck reuses the manuscript figures unchanged. Their sources are PDFs, which a
# browser deck cannot show inline, so each is rasterized at 200 dpi under the same
# stem. Fails loud if any source is missing.
set -euo pipefail

FIG=paper/ctj_manuscript/figures
OUT=presentation/figures

sources=(
  "$FIG/compact/nc_planning_map.pdf"
  "$FIG/compact/grid_accuracy_srs.pdf"
  "$FIG/compact/allocation_distribution_by_design.pdf"
  "$FIG/compact/nc_annual_accuracy_six.pdf"
  "$FIG/fig_colleges_needed.pdf"
  "$FIG/fig_power_vs_n.pdf"
  "$FIG/fig_spread_vs_gain.pdf"
  "$FIG/fig_vs_srs_by_regime_rook.pdf"
  "$FIG/real_sud/population_shares.pdf"
)

for src in "${sources[@]}"; do
  [[ -f "$src" ]] || { echo "ERROR: missing $src" >&2; exit 1; }
  stem=$(basename "$src" .pdf)
  pdftoppm -png -r 200 -singlefile "$src" "$OUT/$stem"
  echo "wrote $OUT/$stem.png"
done

# Clean (title-free) county-to-college map is already a PNG; copy it as is
cp paper/ctj_manuscript/figures/community_college_service_area_clusters.png "$OUT/"
echo "wrote $OUT/community_college_service_area_clusters.png"
# Slide copy with the white margins trimmed so the map and legend fill the slide
python3 -c "
from PIL import Image, ImageChops
im = Image.open('$OUT/community_college_service_area_clusters.png').convert('RGB')
b = ImageChops.difference(im, Image.new('RGB', im.size, (255, 255, 255))).getbbox()
im.crop((max(0, b[0] - 30), max(0, b[1] - 30), min(im.width, b[2] + 30), min(im.height, b[3] + 30))).save('$OUT/cc_service_areas_cropped.png')
"

# Wide grid design panel (written directly to presentation/figures by code/25)
[[ -f "$OUT/design_assignments_panels_wide.png" ]] || { echo "ERROR: run code/25_design_assignment_figures.R" >&2; exit 1; }

Rscript presentation/tools/nc_design_examples.R
Rscript presentation/tools/deck_allocation_figure.R
