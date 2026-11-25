# E8 second half: merge Tables S7 and S12 in the curated SI (sign-off)

Author: Andrew Walther. Drafted 2026-11-25. The body rows come from `code/32_si_performance_table.R`
(`results/si_tables/performance_nine.tex`), which first reproduces every printed cell of Tables S7
and S12 (84 + 54 cells, all exact; `performance_nine_check.txt`). The curated SI `\input`s that
file, so the numbers stay tied to the script. The unabridged build keeps Tables S7 and S12 as now.
The chapter's own `tab:performance` and `tab:eight-design` are untouched, so the verifier's
chapter-to-SI check is unaffected.

| ID | Change (curated build only) | Proposed |
|---|---|---|
| M1 | Table S7 (6.1) becomes the nine-strategy table: rows for all eight designs and SRS in queen-MSE order, the same eight columns as S7 (MSE, SE of MSE ×100, bias, coverage, power, average rank, clusters treated). It carries both labels, so every curated reference to S7 or S12 lands on it. | Caption: "Performance of all eight designs and SRS at $\beta = 1$, pooled over the five incidence configurations, under queen and rook contiguity. Each entry averages 160 scenarios of 250 oracle fits each. Rows follow queen-contiguity mean squared error (MSE). SE of MSE is its Monte Carlo standard error, multiplied by 100. Average rank is the within-block MSE rank among the six designs of the detailed comparison and SRS (1 = lowest). The two variants, Saturation Quadrants and Balanced Halves, are outside that comparison and show a dash. Clusters treated is a mean per assignment. Under rook contiguity, $\beta$ is not identifiable for Checkerboard ($\mathbf{W}\mathbf{Z} = \mathbf{1} - \mathbf{Z}$), and its coefficient on $\mathbf{Z}$ estimates $\beta - \psi$." |
| M2 | Table S12 (6.4) | Unabridged only (its MSE and coverage columns are now in M1; its nine-way ranks stay in the unabridged SI). |
| M3 | 6.4, last sentence: "The MSEs and coverages of the six retained designs are those of Table S7." | Omitted (it would point the table at itself). |
| M4 | 6.1: "For reference, SRS had an average rank over the 160 blocks of 2.600 under queen contiguity and 2.450 under rook contiguity." | Omitted (the table now shows both ranks; code/32 reproduces them). |

Expected effect: about half a page shorter (one table instead of two).
