# Longleaf Quick Reference

## Connect
ssh onyen@longleaf.unc.edu          # requires VPN if off-campus
# Or: https://ondemand.rc.unc.edu   (browser portal, no SSH needed)

## Navigate to your working directory
cd /work/users/a/n/onyen            # replace a/n/onyen with your path
                                    # (first two chars of onyen form subdirs)

## Get your code
git clone https://github.com/you/repo.git
git pull                            # update existing clone

## Load R
module purge
module add r/4.4.0
module list                         # confirm what's loaded

## Install R packages (first time only — interactive R session)
srun -p interact -n 1 --cpus-per-task=1 --mem=8g -t 1:00:00 --pty /bin/bash
# then inside that session:
module add r/4.4.0
R
# > install.packages("tidyverse")   # installs to ~/R/libs, persists across jobs

## Submit jobs
sbatch submit_array.sl              # submit job array
sbatch submit_single.sl             # submit single job

## Monitor
squeue -u onyen                     # your running/queued jobs
squeue -u onyen --start             # estimated start time for queued jobs
scontrol show jobid <jobID>         # full job details

## After job completes
seff <jobID>                        # CPU + memory utilization — USE THIS to tune resources
sacct -j <jobID> --format=JobID,Elapsed,MaxRSS,State    # timing + memory used

## Cancel
scancel <jobID>                     # cancel one job
scancel -u onyen                    # cancel ALL your jobs (careful!)

## Copy results back to Mac (run on your Mac)
scp onyen@longleaf.unc.edu:/work/users/a/n/onyen/repo/results/all_results.rds ~/local/path/
scp -r onyen@longleaf.unc.edu:/path/results/ ~/local/path/   # copy whole folder

## Storage paths
/nas/longleaf/home/onyen            # home dir, 50 GB, backed up, slow I/O
/work/users/a/n/onyen               # work dir, 10 TB, NOT backed up, fast I/O — use this
/users/a/n/onyen                    # archive, 10 TB, NOT backed up

## Resubmit failed array tasks
# Check which tasks failed:
ls results/rds/ | wc -l             # how many completed
# Resubmit specific task IDs only:
sbatch --array=5,12,47,83 submit_array.sl

## Partitions (usually leave unspecified to use default 'general')
# general    — default, CPU-only, up to 11 days
# interact   — interactive sessions, up to 8 hours
# bigmem     — extreme memory jobs (requires special access)
