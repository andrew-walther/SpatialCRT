#!/usr/bin/env bash
# ============================================================
# Script: copy_figures.sh
# Purpose: Copy figures that analysis scripts write under results/ into both
#          manuscript figure folders, under the names the TeX and Quarto sources use.
# Author: Andrew Walther
# Created: 2026-11-02
# Dependencies: bash, cp
# ============================================================
#
# Run from the project root, after re-running the producing scripts:
#   bash paper/tools/copy_figures.sh
#
# Figures written straight into paper/ by their scripts are not listed here:
# code/20 (compact/ grid and NC exhibits, real_sud/regions.png), code/21
# (compact/allocation_distribution_by_design.pdf), code/25 (design panels and
# Table 1 thumbnails) and application/code/manuscript_maps.R (NC maps).
# Fails loud if any source is missing; nothing is copied in that case.
set -euo pipefail

SIX=results/six_design_manuscript
EXH=application/results/real_sud_rev_20261002/exhibits

# source | destination (relative to each figures/ folder) | producing script
map=(
  "$SIX/queen/fig_mse_by_design_6design.pdf|fig_mse_by_design_ranked_6design_queen.pdf|code/14"
  "$SIX/rook/fig_mse_by_design_6design.pdf|fig_mse_by_design_ranked_6design_rook.pdf|code/14"
  "$SIX/queen/fig_biasvar_6design.pdf|fig_biasvar_6design_queen.pdf|code/14"
  "$SIX/rook/fig_biasvar_6design.pdf|fig_biasvar_6design_rook.pdf|code/14"
  "$SIX/queen/si_figures/si_fig_cd_diagram_rank.pdf|app_cd_diagram_rank_queen.pdf|code/14"
  "$SIX/rook/si_figures/si_fig_cd_diagram_rank.pdf|app_cd_diagram_rank_rook.pdf|code/14"
  "$SIX/queen/si_figures/si_fig_tau_coverage.pdf|app_tau_coverage_queen.pdf|code/14"
  "$SIX/rook/si_figures/si_fig_tau_coverage.pdf|app_tau_coverage_rook.pdf|code/14"
  "$SIX/fig_coverage_by_design_6design_queen.pdf|fig_coverage_by_design_6design_queen.pdf|code/12"
  "$SIX/fig_coverage_by_design_6design_rook.pdf|fig_coverage_by_design_6design_rook.pdf|code/12"
  "$SIX/fig_tau_sensitivity_6design_combined.pdf|fig_tau_sensitivity_6design_combined.pdf|code/12"
  "results/srs_benchmark/fig_vs_srs_by_regime_queen.pdf|fig_vs_srs_by_regime_queen.pdf|code/15"
  "results/srs_benchmark/fig_vs_srs_by_regime_rook.pdf|fig_vs_srs_by_regime_rook.pdf|code/15"
  "results/methods_ext_figures/fig_power_vs_n.pdf|fig_power_vs_n.pdf|code/23"
  "results/methods_ext_figures/fig_colleges_needed.pdf|fig_colleges_needed.pdf|code/23"
  "results/methods_ext_figures/fig_spread_vs_gain.pdf|fig_spread_vs_gain.pdf|code/23"
  "$EXH/population_shares_journal.pdf|real_sud/population_shares.pdf|application/code/real_sud_summary.R"
  "$EXH/primary_setting_ratios_journal.pdf|real_sud/primary_setting_ratios.pdf|application/code/real_sud_summary.R"
  "$EXH/primary_yearly_mse_journal.pdf|real_sud/primary_yearly_mse.pdf|application/code/real_sud_summary.R"
  "$EXH/refined_allocation_risk_journal.pdf|real_sud/refined_allocation_risk.pdf|application/code/real_sud_summary.R"
)

# Check every source first so a missing one copies nothing
for row in "${map[@]}"; do
  IFS='|' read -r src dst script <<< "$row"
  [[ -f "$src" ]] || { echo "ERROR: missing $src (run $script)" >&2; exit 1; }
done

for row in "${map[@]}"; do
  IFS='|' read -r src dst script <<< "$row"
  for folder in paper/ctj_manuscript/figures paper/thesis_chapter/figures; do
    cp "$src" "$folder/$dst"
  done
  echo "copied $src -> {ctj_manuscript,thesis_chapter}/figures/$dst"
done
