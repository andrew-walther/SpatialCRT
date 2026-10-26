Andrew Walther, Tonya Van Deinse & Dr. Feng-Chang Lin\
University of North Carolina at Chapel Hill, Department of Biostatistics\
McGavran-Greenberg Hall, 3101, Chapel Hill, NC 27599

Editorial Office\
*BMC Medical Research Methodology*

October 25, 2026

Re: Correction to "Spatial Cluster Randomized Trials - Sampling Design with Spillover Effects & Spatial Dependence" (doi:10.1186/s12874-026-02968-0)

Dear Editor,

We are writing to request a correction to our article. After it was accepted, while adapting it for the first author's dissertation, we re-checked the simulation code and found two coding errors. Both affect only the results for the 3×3 and 3×4 grids. The 2×4 results are correct, and the article describes the methods correctly. The errors were in how the code carried them out. Correcting them does not change the article's conclusions or its recommendation for the North Carolina probation studies. The corrected results support those conclusions more consistently across the three grids than the published results did.

**What was wrong**

1. *Block stratified allocations in the summary.* The script that summarized the per-allocation results into Tables 2 and 3 and Figures 6 and 7 used the block stratified sampling (BSS) allocations of the 2×4 grid for all three grids. The BSS columns for the 3×3 and 3×4 grids therefore summarized two allocations that do not meet the block stratification condition, instead of the six valid 3×3 allocations and the two 3×4 checkerboards.
2. *Spillover loop in the 3×3 and 3×4 simulations.* Under control-only spillover, the loop that assigns spillover stopped at district 8, a bound carried over from the 2×4 script. Intervention district 9 on the 3×3 grid, and districts 9 to 12 on the 3×4 grid, therefore never passed spillover to neighboring control clusters. This affected the simple random sampling (SRS) results for control-only spillover on those grids. It did not affect any BSS allocation.

**How we checked**

We first reproduced the published per-allocation results from the original code, to within 1e-8, which confirmed both errors. We then reran the 16 affected control-only scenarios with the loop corrected, using the same random seeds, and rebuilt both tables and figures with each grid's correct BSS allocations. We also found that one stored result file (3×3 grid, spillover into both arms, ψ = 0.8, ρ = 0.01) had been generated with a different random number stream from the others. We regenerated it in the same way as the rest, which changes the SRS values in that row of Table 2 only in the third decimal place.

**What changes**

- Tables 2 and 3 (3×3 and 3×4 grids): the BSS columns, and the SRS rows for control-only spillover, are replaced.
- Figures 6 and 7: replaced with figures built from the corrected results.
- Text: the results paragraphs for the 3×3 and 3×4 grids are rewritten from the corrected tables. The Discussion paragraph describing a "checkerboard paradox" is removed. That pattern came from the incorrect allocations. In its place, a short paragraph explains the collinearity between treatment and spillover in checkerboard allocations, which the Limitations section already describes. A few summary sentences in the Abstract, Results, Recommendations, Discussion and Future Considerations are adjusted so that they match the corrected results.

A marked-up copy of the manuscript showing every change is enclosed, together with the corrected manuscript and figures.

**The conclusions do not change**

With the errors corrected, the 3×4 grid shows the same pattern as the 2×4 grid in all 16 scenarios, and the 3×3 grid matches it in 12 of 16. On every grid, SRS has the lower average mean squared error when spillover reaches both treatment arms, and BSS has the lower error when spillover flows into control clusters and spatial dependence is present. The only difference between grids is control-only spillover without spatial dependence, where BSS has the lower error on the 3×3 grid. The article's central message is unchanged: SRS can achieve lower average error, BSS limits the risk of a poor allocation, and the choice depends on the anticipated spillover mechanism, with BSS most advantageous when spillover flows mainly into control clusters.

If the article can still be updated before its final publication, we ask that the enclosed corrected files replace the accepted version. Otherwise, we ask that a Correction be published under the journal's policy on corrections. We have drafted text for that note below.

All authors have read and approved this request. We apologize for the errors and thank you for your help in correcting the record.

Please address correspondence to the corresponding author, Dr. Feng-Chang Lin, at flin@bios.unc.edu, and to me at awalther@unc.edu.

Sincerely,

Andrew Walther, on behalf of all authors

\newpage

**Proposed text for the Correction note**

Correction: Spatial Cluster Randomized Trials - Sampling Design with Spillover Effects & Spatial Dependence

Following publication of the original article [1], the authors identified two coding errors that affected the results for the 3×3 and 3×4 grids. First, the summaries of block stratified sampling (BSS) in Tables 2 and 3 and Figures 6 and 7 used the BSS allocations of the 2×4 grid instead of each grid's own BSS allocations. Second, under control-only spillover, the simulation for the 3×3 and 3×4 grids did not assign spillover from intervention districts numbered above 8, which affected the simple random sampling (SRS) results for that setting. The 2×4 results are unaffected.

Tables 2 and 3 and Figures 6 and 7 have been replaced with corrected versions. The results text for the 3×3 and 3×4 grids has been revised to match them, and a few summary statements in the Abstract, Results, Recommendations, Discussion and Future Considerations have been adjusted accordingly. The Discussion paragraph describing a "checkerboard paradox" has been removed, because that pattern arose from the incorrect allocations. It has been replaced by a short explanation of the collinearity between treatment and spillover in checkerboard allocations.

The conclusions of the article are unchanged. With the corrections, the 3×4 grid shows the same pattern as the 2×4 grid in every scenario. SRS achieves lower average mean squared error when spillover reaches both treatment arms, and BSS achieves lower error when spillover flows into control clusters and spatial dependence is present.

The original article has been corrected.

**Reference**

1. Walther A, Van Deinse T, Lin F-C. Spatial Cluster Randomized Trials - Sampling Design with Spillover Effects & Spatial Dependence. BMC Med Res Methodol. 2026. https://doi.org/10.1186/s12874-026-02968-0
