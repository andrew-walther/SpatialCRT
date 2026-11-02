# data/

| File | What | Source |
|---|---|---|
| `NC_Judicial_District_Assignments.csv` | The 31 Division of Community Supervision districts (1–30, with 19A and 19B), their judicial division, and an `Assignment` column | Written by `code/NCTrtMaps.Rmd` (a seed-42 random 50/50 draw within each division). This draw is **not** the allocation shown in the paper's BSS map; that map uses the hard-coded `treatment_districts_list` in `NCTrtMaps.Rmd` |
| `NC_DCS_District_Counties.csv` | County → district → division, all 100 counties | Transcribed by hand on 2026-11-02 from the NC Department of Adult Correction's Community Supervision district map (updated 2025-01-07), <https://www.dac.nc.gov/document/community-supervision-district-map/open>. The PDF is an image, so there is no machine-readable source |
| `NC_LME_MCO_Counties_2025.csv` | County → LME/MCO, all 100 counties, four LME/MCOs | NCDHHS LME/MCO directory as archived on 2025-05-31, <http://web.archive.org/web/20250531105759/https://www.ncdhhs.gov/providers/lme-mco-directory>. The live page (<https://www.ncdhhs.gov/providers/lme-mco-directory>) now lists three, after the Vaya–Partners merger; the paper uses the four |

`presentation/tools/make_figures.R` reads the two county files to draw the North
Carolina maps (county shapes from the US Census via `tigris`).

## Checks and known differences

- Both county files cover the same 100 county names as the Census 2023 county
  shapes (checked in `make_figures.R`, which stops otherwise).
- The redrawn MCO map matches the paper's MCO map county for county.
- The DCS district map matches the paper's map except two counties, where the
  DAC map was followed:
  - **Chowan** is in district 1 here; the paper's map has it in district 6. On the
    DAC map the label "Chowan" is printed inside district 6, but the county itself
    (east of the Chowan River) lies inside district 1's boundary.
  - **Dare** is in district 1 here; the paper's map has it in district 2. On the
    DAC map a district boundary separates Dare from Tyrrell (district 2).
  Changing either is a one-line edit to `NC_DCS_District_Counties.csv`.
- `code/NCTrtMaps.Rmd` (frozen) expects a `county` column in
  `NC_Judicial_District_Assignments.csv`, which that file does not have, so the
  Rmd does not run as committed. The county-level file it used was not kept.
