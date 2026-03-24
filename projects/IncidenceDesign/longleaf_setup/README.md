# Longleaf HPC R Simulation Template

A ready-to-use template for running R simulations on UNC's Longleaf HPC cluster using SLURM job arrays. Designed for researchers running many seeds or parameter combinations in parallel.

## File structure

```
.
├── submit_array.sl          # SLURM job array script (primary — use this)
├── submit_single.sl         # SLURM single-run script (for one long job)
├── LONGLEAF_CHEATSHEET.md   # Common commands reference
├── R/
│   ├── simulation.R         # Your simulation script (edit this)
│   └── aggregate_results.R  # Combines output from all tasks after run
├── logs/                    # SLURM output/error logs (auto-created)
└── results/
    ├── rds/                 # Per-task .rds output files
    ├── csv/                 # Per-task .csv output files
    └── figures/             # Per-task plot files
```

## Quickstart

### 1. Copy this template into your repo

Drop these files into your existing Git repo. Commit and push.

### 2. On Longleaf: clone/pull your repo into /work

```bash
ssh onyen@longleaf.unc.edu
cd /work/users/a/n/onyen           # your work directory
git clone https://github.com/you/repo.git
cd repo
mkdir -p logs results/rds results/csv results/figures
```

### 3. Install your R packages (once, interactively)

```bash
srun -p interact -n 1 --cpus-per-task=1 --mem=8g -t 1:00:00 --pty /bin/bash
module add r/4.4.0
R
# install.packages(c("tidyverse", "your_packages_here"))
# q()
exit
```

### 4. Adapt simulation.R

Edit `R/simulation.R`. Key lines to understand:
- `args <- commandArgs(trailingOnly = TRUE)` — receives the task ID from SLURM
- `task_id` is the SLURM array index (1, 2, 3, ..., N) — use it as your seed or parameter index
- Save all outputs with `task_id` in the filename to avoid file collisions across parallel tasks

### 5. Configure and submit

Edit `submit_array.sl`:
- Set `--array=1-N` to the number of runs you want
- Set `--mem=` and `--time=` appropriately (start conservative; tune with `seff` after first run)
- Set `--mail-user=` to your UNC email

Then submit:
```bash
sbatch submit_array.sl
```

### 6. Monitor

```bash
squeue -u onyen         # watch progress
seff <jobID>            # check resource usage after completion
```

### 7. Aggregate results

After all tasks complete, combine outputs:
```bash
Rscript R/aggregate_results.R
```

Or copy results back to your Mac first:
```bash
# Run on your Mac:
scp -r onyen@longleaf.unc.edu:/work/users/a/n/onyen/repo/results/ ~/local/
```

## Resource sizing guide

| Simulation type         | `--mem` | `--cpus-per-task` | `--time`     |
|-------------------------|---------|-------------------|--------------|
| Small MC (< 10k obs)    | 4g      | 1                 | 1:00:00      |
| Medium (10k–1M obs)     | 8g      | 1                 | 4:00:00      |
| Large / complex model   | 16–32g  | 1                 | 12:00:00     |
| Parallel R (future etc) | 16g     | 4–8               | 4:00:00      |

**Always run `seff <jobID>` after your first job to see actual memory and CPU usage, then adjust.**

## Notes

- `results/` and `logs/` are intentionally git-ignored (add them to `.gitignore`)
- The `%10` suffix in `--array=1-100%10` limits concurrent tasks to 10 — useful to be a good cluster citizen on large arrays
- VPN required to SSH into Longleaf from off-campus
- Never run code directly on the login node — always use `sbatch` or `srun`
