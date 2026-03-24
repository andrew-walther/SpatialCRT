# Longleaf HPC Setup for IncidenceDesign

Run the MLE simulation (2,560 scenarios) on UNC's Longleaf cluster using SLURM job arrays. Each scenario runs as an independent task (~1-5 min), enabling massive parallelism.

## Directory Structure

```
longleaf_setup/
  simulation.R           # Per-scenario simulation script (called by SLURM)
  submit_array.sl        # SLURM job array submission script
  aggregate_results.R    # Combines per-scenario outputs after completion
  install_packages.R     # First-time R package setup
  submit_single.sl       # Alternative: single long-running job (not used)
  LONGLEAF_CHEATSHEET.md # Quick reference for Longleaf commands
  CLAUDE_CODE_BRIEFING.md # Context for Claude Code sessions
  logs/                  # SLURM stdout/stderr (created at runtime, git-ignored)
  results/rds/           # Per-scenario .rds files (created at runtime, git-ignored)
```

The simulation sources `../code/01-04*.R` modules — it does NOT modify them.

## How It Works

The 2,560 scenarios are a full factorial of all parameter combinations:

| Parameter | Values | Count |
|-----------|--------|:-----:|
| Incidence config | iid, spatial(0.2), spatial(0.5), poisson(0.2), poisson(0.5) | 5 |
| Neighbor type | rook, queen | 2 |
| Rho | 0.00, 0.01, 0.20, 0.50 | 4 |
| Gamma | 0.5, 0.6, 0.7, 0.8 | 4 |
| Spillover type | control_only, both | 2 |
| Design | 1-8 | 8 |
| **Total** | | **2,560** |

`simulation.R` builds this grid via `expand.grid()`. SLURM passes a task ID (1-2560), which indexes into the grid. Each task:
1. Rebuilds the spatial grid and generates incidence/epsilon with deterministic seeding
2. Runs one scenario's estimation (25 design x 10 outcome resamples = 250 iterations)
3. Saves a 1-row .rds file to `results/rds/scenario_<id>.rds`

## One-Time Setup

```bash
ssh awalther@longleaf.unc.edu
cd /work/users/a/w/awalther
git clone https://github.com/<username>/SpatialCRT.git
cd SpatialCRT/projects/IncidenceDesign/longleaf_setup
mkdir -p logs results/rds

# Install R packages (interactive session)
srun -p interact -n 1 --mem=8g -t 1:00:00 --pty /bin/bash
module add r/4.4.0
Rscript install_packages.R
exit
```

## Running the Simulation

```bash
cd /work/users/a/w/awalther/SpatialCRT/projects/IncidenceDesign/longleaf_setup
git pull  # if code was updated

# Step 1: Test with one task
sbatch --array=1-1 submit_array.sl
squeue -u awalther
# After completion:
seff <jobID>               # check memory/time usage
cat logs/slurm-*_1.out     # check output
cat logs/slurm-*_1.err     # check errors

# Step 2: Submit all 2,560 tasks
sbatch submit_array.sl
squeue -u awalther         # monitor progress

# Step 3: Aggregate results
module add r/4.4.0
Rscript aggregate_results.R
# -> saves to ../results/longleaf/sim_data/sim_results_MLE_combined_<timestamp>.rds
```

## Getting Results to Your Mac

```bash
# Run on your Mac — copies the entire longleaf/ results tree:
scp -r awalther@longleaf.unc.edu:/work/users/a/w/awalther/SpatialCRT/projects/IncidenceDesign/results/longleaf/ \
    ~/GithubProjects/SpatialCRT/projects/IncidenceDesign/results/longleaf/
```

## Results Layout

Longleaf results are kept separate from local results under `results/longleaf/`:

```
results/
  sim_data/                         # LOCAL simulation data (.rds files)
  07_results_summary.pdf            # LOCAL reports
  09_MLE_design_recommendation_report.pdf
  mle_per_config/                   # LOCAL per-config figure PDFs
  longleaf/                         # ALL Longleaf outputs
    sim_data/                       #   aggregated .rds files
    reports/                        #   rendered reports (07, 09, etc.)
    mle_per_config/                 #   per-config figure PDFs
```

**Loading Longleaf results in R:**
```r
# IMPORTANT: always use this path when working with Longleaf results
longleaf_results <- load_latest_results(
  results_dir = "results/longleaf/sim_data"
)
```

**Generating reports from Longleaf results:**

The existing report scripts (`06_visualizations.R`, `08_design_recommendations.R`) work
with Longleaf results — just point them at the right data and output directory:

```r
source("code/06_visualizations.R")

# Load Longleaf results
results <- load_latest_results(results_dir = "results/longleaf/sim_data")

# Generate visualizations into the longleaf reports directory
run_all_visualizations(
  results     = results,
  results_dir = "results/longleaf",
  estimation_mode = "MLE_combined",
  output_pdf  = TRUE
)

# Generate design recommendations
source("code/08_design_recommendations.R")
run_recommendation_report(
  results         = results,
  estimation_mode = "MLE_combined",
  results_dir     = "results/longleaf"
)
```

Reports will be saved to `results/longleaf/reports/` and figures to
`results/longleaf/mle_per_config/`, keeping them fully separate from
local-run outputs.

## Resubmitting Failed Tasks

```bash
# Check which tasks completed:
ls results/rds/ | wc -l

# aggregate_results.R reports missing IDs. Resubmit specific ones:
sbatch --array=5,12,47 submit_array.sl
```

## Scaling Up (Phase 2)

After validating the pipeline works, increase precision by editing `simulation.R`:
```r
n_design_resamples  <- 50    # was 25
n_outcome_resamples <- 50    # was 10
```
And increase the time limit in `submit_array.sl`:
```
#SBATCH --time=01:00:00     # was 00:30:00
```

## Resource Sizing

| Resamples | Iterations | Est. time/task | `--time` | `--mem` |
|-----------|:----------:|:--------------:|:--------:|:-------:|
| 25 x 10 (default) | 250 | ~2 min | 00:30:00 | 4g |
| 50 x 25 | 1,250 | ~10 min | 00:30:00 | 4g |
| 50 x 50 | 2,500 | ~25 min | 01:00:00 | 4g |

Always run `seff <jobID>` after your first job to check actual usage and adjust.
