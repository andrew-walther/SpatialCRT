#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

APPLICATION_RUN_FULL=TRUE Rscript application/run_application_profiles.R full

