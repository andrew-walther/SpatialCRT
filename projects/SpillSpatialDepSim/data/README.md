# data/

| File | What | Source |
|---|---|---|
| `NC_County_Assignments.csv` | County → Division of Community Supervision district → judicial division, all 100 counties | The county file behind the paper's maps, read off the NC Department of Adult Correction's Community Supervision district map (<https://www.dac.nc.gov/document/community-supervision-district-map/open>; the 2024-12-10 and 2025-01-07 versions carry the same map image). Copied verbatim on 2026-11-02 from the author's OneDrive (`Project 1 - Spillover Effects/NC DOC Training/NC_County_Assignments.csv`, dated 2025-10-16) |
| `NC_Judicial_District_Assignments.csv` | The 31 districts (1–30, with 19A and 19B), their judicial division, and an `Assignment` column | Written by `code/NCTrtMaps.Rmd` (a seed-42 random 50/50 draw within each division). This draw is **not** the allocation in the paper's BSS map, which uses the hard-coded `treatment_districts_list` in `NCTrtMaps.Rmd`, chosen by hand as the most BSS-like allocation on the state's irregular map |
| `NC_LME_MCO_Counties_2025.csv` | County → LME/MCO, all 100 counties, four LME/MCOs | NCDHHS LME/MCO directory as archived on 2025-05-31, <http://web.archive.org/web/20250531105759/https://www.ncdhhs.gov/providers/lme-mco-directory>. The live page (<https://www.ncdhhs.gov/providers/lme-mco-directory>) now lists three, after the Vaya–Partners merger; the paper uses the four |

`presentation/tools/make_figures.R` reads `NC_County_Assignments.csv` and
`NC_LME_MCO_Counties_2025.csv` to draw the North Carolina maps (county shapes from
the US Census via `tigris`, cartographic boundaries, 2023).

## Checks

- Both county files cover the same 100 county names as the Census county shapes
  (`make_figures.R` stops otherwise).
- The redrawn district and MCO maps match the paper's maps district for district and
  county for county. The coastline differs: the paper's maps used tigris's detailed
  county shapes, which count the sounds as land.
- The DAC map is an image, and three counties sit on ambiguous boundaries. A fresh
  reading on 2026-11-02 put Chowan and Dare in district 1 and Onslow in 4; this file
  (the paper's reading) has Chowan in 6, Dare in 2 and Onslow in 3. The paper's reading
  is kept for continuity; on a close look Onslow's is also the better one.

## Rerunning `code/NCTrtMaps.Rmd`

The frozen Rmd reads `../data/NC_Judicial_District_Assignments.csv`, which has no
`County` column, so it fails as committed. It was run with the county file above.
A copy with that path pointed at `NC_County_Assignments.csv` renders end to end
(checked 2026-11-02), given the `janitor` package (`install.packages("janitor")`).
