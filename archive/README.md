# Archive

This directory contains legacy and exploratory code preserved for reference.
**None of these files are actively maintained.** For current work, see `projects/`.

---

## Contents

### `PreliminarySpatialSim/`
Early exploratory simulation work predating both main projects.
Contains initial replication attempts and grid/estimation experiments.
- `JarvisReplicate/` — Attempted replication of Jarvis et al. spatial simulation
- `SpatialSimEstimates/` — Early estimator explorations (binary, linear, exponential outcomes)
- `SpatialSimGrid/` — Early spatial grid setups (Sept–Oct 2024 scripts)

### `SpatialSim_Unified/`
An early unification attempt for the NCdocSpatialSim scripts (now SpillSpatialDepSim).
Superseded by the Plan B modular scripts (`code/01_*` through `code/05_*`) in
`projects/SpillSpatialDepSim/`. Not recommended for use.

### `OutcomeIncidenceDesign_Legacy/`
The monolithic predecessor to `projects/IncidenceDesign/`. These Rmds were the original
all-in-one simulation + analysis scripts (~600 lines each) that were refactored into
the modular numbered scripts. Preserved untouched for historical reference.
- `SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd` — Original monolithic simulation
- `SpatialCRT_Incidence_Sim.Rmd` — Earlier variant
- `SpatialCRT_Incidence_Simulation_MLE_Results.RData` — MLE results from original run

---

## Note

All files here have full git history. Use `git log --follow <file>` to trace origins.
