# Claude Code Session Briefing: Longleaf HPC Setup

## Who I am & my experience level
I'm a researcher at UNC Chapel Hill. I'm new to HPC and have never submitted a job to a cluster before. I'm comfortable in R and RStudio, I use Git/GitHub for version control, and I've used Claude Code before — but anything involving the Linux command line, SSH, or SLURM is unfamiliar territory. Please explain things clearly and don't assume I know cluster-specific terminology.

## What I'm trying to do
I want to run my R simulations on UNC's Longleaf HPC cluster instead of my Mac, because some simulations take 12+ hours and require my laptop to stay awake the whole time. The HPC can run them faster and unattended.

## What I already have
In a prior Claude.ai session, I got a full orientation to Longleaf and a template set of files. Those files are in this repo. Here's what they do:

- `submit_array.sl` — the main SLURM job submission script. SLURM is the job scheduler Longleaf uses. This script tells the cluster how many CPUs, how much memory, and how long to reserve for my job, then runs my R script. The "array" part means it runs the same R script many times in parallel, each with a different task ID (think of it like running `Rscript simulation.R 1`, `Rscript simulation.R 2`, ... `Rscript simulation.R 100` all at the same time on different machines).
- `submit_single.sl` — same idea but for one single long-running job instead of an array.
- `R/simulation.R` — a template R script showing how to receive the task ID from SLURM and use it as a random seed or parameter index. My actual simulation code needs to replace the placeholder example in the middle.
- `R/aggregate_results.R` — a script to combine all the individual output files from each parallel task into one dataset after the run finishes.
- `LONGLEAF_CHEATSHEET.md` — a reference card of every command I'll need on Longleaf.
- `README.md` — step-by-step instructions for the full workflow.

## The overall workflow (plain English)
1. I write/edit R code on my Mac in RStudio as usual, then push to GitHub.
2. On Longleaf (via SSH in Terminal, or a browser portal), I `git pull` to get my latest code onto the cluster's storage.
3. I edit the SLURM script to request the right resources (memory, time, number of parallel runs).
4. I type `sbatch submit_array.sl` to submit the job. Longleaf queues it and runs it — I can close my laptop.
5. I check progress with `squeue -u onyen`. When done, I use `seff <jobID>` to see how much memory/CPU was actually used (to tune future jobs).
6. I run `aggregate_results.R` to combine outputs, then either `scp` results back to my Mac or push them to Git from Longleaf.

## Key Longleaf facts to know
- Cluster name: Longleaf, at UNC Chapel Hill
- Scheduler: SLURM
- Login: `ssh onyen@longleaf.unc.edu` (requires UNC VPN if off-campus)
- Browser portal: https://ondemand.rc.unc.edu (no SSH needed, can launch RStudio on the cluster)
- R version available: `r/4.4.0` (loaded with `module add r/4.4.0`)
- Where to run jobs: `/work/users/<first-letter>/<second-letter>/<onyen>/` — NOT the home directory (home is slow and small)
- NEVER run code directly on the login node — always submit via `sbatch`
- General partition is the default (no `--partition` flag needed for standard CPU jobs)
- Max runtime: 11 days per job

## What I need help with in this session
Please look at my actual simulation R script(s) in this repo and help me:

1. **Understand what my simulation is doing** — explain it back to me in plain terms so I can confirm you've read it correctly.
2. **Adapt `R/simulation.R`** to wrap my actual simulation code, handling the task ID as seed and/or parameter index as appropriate.
3. **Configure `submit_array.sl`** with appropriate `--mem`, `--time`, and `--array` values based on what my code actually does.
4. **Identify any issues** that might cause problems on a Linux HPC (e.g. hardcoded Mac file paths, packages that need special installation, GUI-dependent code like `View()` or `RStudio::` calls, `setwd()` calls that won't work on the cluster).
5. **Walk me through each step** of getting this running on Longleaf for the first time, assuming I've never SSHed into a cluster before.

Please ask me clarifying questions if anything about my code is ambiguous before making changes.
