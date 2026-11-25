# Phases 4e and 3c: new curated-SI sections and two CTJ sentences (sign-off)

Author: Andrew Walther. Drafted 2026-11-25. Nothing below is applied yet.

Sources. 4e: `results/incidence_coef_2026-11/full/summary/performance.csv` (numbers recomputed
from it by script in this session, oracle fit, `Intervention_Grid == "All_Beta"`), the findings
(`docs/findings/2026-11-25-incidence-coefficient-findings.md`) and the review. 3c:
`results/igs_sensitivity_2026-11/full/summary/igs_summary.txt`, `igs_burden_share.csv` and the
findings and rerun review. Every number below appears in those files. MSE values have 3 decimals,
ratios to SRS 2 per the rule except where the source table gives 3 (see D3).

Slots (from the 2e outline): 4e is new subsection 3.2 after "Identification diagnostics"; 3c is
Methods 2.2 after "The two related variants" and a new Results section placed after "Robustness,
the SRS benchmark and allocation-specific risk" and before "Case study: accuracy, inference and
reach". All new text prints in both SI builds (no `\iffullsi`).

## Decisions first

| ID | Decision | Recommendation and reason |
|---|---|---|
| D1 | 4e figure: add the queen rejection figure (`incidence_rejection_queen.pdf`, code/31) as a new curated figure, or text only | **Add it.** One figure shows the whole power story by design, τ, configuration and regime; it has no in-figure title and uses τ. It goes into the paper folders through a new line in `copy_figures.sh`. |
| D2 | 3c figure: `fig_igs_re_vs_srs.pdf` (code/30), also via `copy_figures.sh` | **Add it** (planned in 3b). |
| D3 | Ratio decimals in 3c. The IGS summaries print ratios to SRS with 3 decimals (0.551, 8.126). The rule says 2 for ratios to SRS. | **Two decimals** (0.55, 8.13), per the 2026-11-25 decision; the proposed text below already uses two. Ratios to the original design or the twin also get two. |
| D4 | CTJ sentences C1 (4e) and C2 (3c) | **Both**; together they add about 36 words (3,435 → about 3,471 of 3,500). |

## 4e: Estimating the incidence coefficient after a trial

### N1. New subsection 3.2 (after Section 3.1, before `\FloatBarrier` and Section 4)

Current: none (empty slot).

Proposed:

```latex
\subsection{Estimating the incidence coefficient after a trial}\label{incidence-coefficient}

The incidence coefficient \(\tau\) in Equation~\ref{eq:outcome} is a regression coefficient, not a correlation. Incidence \(\mathbf{X}\) lies between 0 and 1 and \(\sigma = 1\), so \(\tau\) is the shift in the outcome, in residual standard deviations, from the lowest to the highest incidence, holding treatment and spillover fixed and before the spatial multiplier spreads it to neighbors. At \(\tau = 1\) that shift equals the primary intervention effect \(\beta = 1\). The value \(\tau = 0\) means that incidence has no association with the outcome, and \(\tau\) can exceed 1 or be negative.

Before a trial, the investigator chooses \(\tau\) for planning, because there are no outcomes yet. After the trial, the investigator has both the prior incidence used for allocation and the measured outcomes, so \(\tau\) can be estimated in the same spatial lag model that estimates \(\beta\). Both grid models already include \(\mathbf{X}_k\) (Section~\ref{estimation}), so every grid fit gives \(\hat\tau\) beside \(\hat\beta\). The estimate describes how prior incidence is associated with the later outcome after adjustment for treatment and spillover. Incidence is not randomized, so \(\hat\tau\) is an adjusted association rather than an effect of incidence.

To see how well \(\tau\) can be estimated, and whether its value changes the design comparison, we repeated the full grid benchmark at \(\tau = 0\), 0.25, 0.5, 1 and 2. The incidence surfaces, assignments and outcome noise were the same at every value of \(\tau\), so only the incidence term of the outcome changed. The 95\% interval for \(\tau\) has the same \(t\)-based, degrees-of-freedom-corrected form as the interval for \(\beta\), and power is the proportion of intervals that exclude zero, a two-sided test at level 0.05. At \(\tau = 1\) the rerun reproduced every stored intervention result of the benchmark exactly.

With the oracle model, the largest absolute bias of \(\hat\tau\) across designs, incidence configurations, neighbor types, regimes and values of \(\tau\) was 0.023. The intervals covered \(\tau\) in 93.8--95.8\% of fits, and at \(\tau = 0\) they excluded zero in 4.3--6.2\% of fits. Without the spillover term, bias reached 0.137, because omitting spillover can shift the fitted association with incidence.

The design affected the precision of \(\hat\tau\) far more than the value of \(\tau\) did. On the queen grid with independent uniform incidence, the mean model standard error of \(\hat\tau\) was 0.667--0.670 for High Incidence Focus and 0.336--0.337 for SRS, and every other design was between 0.334 and 0.341. At \(\tau = 1\), the test rejected in 29.3--30.3\% of High Incidence Focus trials and in 81.1--81.9\% of SRS trials across the two spillover regimes (Figure~\ref{fig:incidence-rejection}). High Incidence Focus reached that level of power only at \(\tau = 2\), where it rejected in 81.3--81.5\% of trials. It had the largest standard error of any design in every incidence configuration, neighbor type and regime. High Incidence Focus treats the clusters above the median incidence, so its mean absolute correlation between \(\mathbf{Z}\) and \(\mathbf{X}\) was 0.866, against 0.079--0.080 for SRS, and little variation in \(\mathbf{X}\) remains once treatment is in the model. For uniform \(\mathbf{X}\), splitting treatment at the median leaves one quarter of the variance of \(\mathbf{X}\) within the arms. In a simple regression with residual standard deviation 1, this doubles the standard error of the incidence coefficient from \(\sqrt{12/100} = 0.35\) to 0.69. Incidence-Guided Saturation also uses incidence, but by region and with every region partly treated, so its correlation was 0.109--0.110 and its standard error 0.339--0.341, close to that of SRS.

Changing \(\tau\) barely moved \(\hat\beta\). Over the five values of \(\tau\), no design's mean MSE for \(\beta\) changed by more than 3.4\% from its value at \(\tau = 1\) in any incidence configuration, neighbor type or regime. Both fitted models include \(\mathbf{X}\), so a change in \(\tau\) is absorbed by the coefficient on \(\mathbf{X}\) and reaches \(\hat\beta\) only through the estimate of \(\rho\). The spillover regime still decided which designs led. Under control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest MSE in every comparison. Under both-arms spillover, the best design's MSE was at most 5\% below that of SRS. The lowest-MSE design changed with \(\tau\) in 6 of the 20 comparisons, and in 5 of these the two saturation designs, whose MSEs were nearly equal, exchanged first place.

The oracle results assume the correct form of spillover, and the non-oracle bias shows that \(\hat\tau\) depends on it. The grid used one number of clusters, ten incidence surfaces per configuration and nonnegative values of \(\tau\), and \(\hat\tau\) is on the scale of incidence after it was scaled to the unit interval. The primary case-study model omits incidence, because the education outcome is not expected to depend on it, so these results give no estimate of \(\tau\) for that trial.

\begin{figure}[!htb]
\centering
\includegraphics[width=0.9\linewidth]{../figures/fig_incidence_rejection_queen.pdf}
\caption{Two-sided rejection rate of the 95\% interval for the incidence coefficient \(\tau\) under queen contiguity with the oracle model, by incidence configuration (rows) and spillover regime (columns), averaged over the five values of \(\beta\). At \(\tau = 0\) the rate is the test size, and the dashed line marks 0.05. High Incidence Focus has the lowest rejection rate at every \(\tau > 0\), and the other designs and SRS overlap.}
\label{fig:incidence-rejection}
\end{figure}
```

Notes: about 1.5 pages with the figure, as planned. "3.4%" is the largest change over the full
β grid (3.7% in the β = 1 slice alone; not printed). "0.137" is the non-oracle maximum over the
full β grid (0.160 in the β = 1 slice). "At most 5% below" is the smallest best-to-SRS ratio,
0.950 (rook; 0.961 queen).

### N2. `paper/tools/copy_figures.sh`: one new map line (not manuscript text)

`"results/incidence_coef_2026-11/full/summary/incidence_rejection_queen.pdf|fig_incidence_rejection_queen.pdf|code/31_incidence_coefficient_summary.R"`

### N3. SI 1.1, pointer where τ is first set (optional)

Current: "\(\mathbf{X}\) is the baseline incidence covariate with coefficient \(\tau = 1\), and \(\sigma = 1\)."

Proposed: "\(\mathbf{X}\) is the baseline incidence covariate with coefficient \(\tau = 1\), and \(\sigma = 1\). Section~\ref{incidence-coefficient} explains how \(\tau\) can be estimated after a trial."

**Recommend accept** (the reader meets τ here first). Rejecting it leaves the sentence as is.

### C1. CTJ, Methods (grid study paragraph), pointer sentence

Current: "The investigator specifies $\tau$ according to whether incidence is expected to shift the outcome. It is not estimated at the design stage, and we set $\tau=1$ in the grid simulation."

Proposed: "The investigator specifies $\tau$ according to whether incidence is expected to shift the outcome. It is not estimated at the design stage, and we set $\tau=1$ in the grid simulation. After a trial, the same model estimates it, less precisely when treatment follows incidence (SI Section~\ref{SI-incidence-coefficient})."

Shorter alternative (C1b, about 12 words): "After a trial, the same model can estimate it (SI Section~\ref{SI-incidence-coefficient})."

**Recommend C1** (+15 words): it carries the one design result, which connects to the
High Incidence Focus sensitivity sentence in Results. Chapter 3 needs the same sentence (2i).

## 3c: Sensitivity of incidence guidance to regions and saturation spread

### N4. New Methods subsection 2.2 (after Table `tab:excluded-designs` and its `\FloatBarrier`)

Current: none (empty slot).

Proposed:

```latex
\subsection{Variants of incidence guidance}\label{igs-variants-methods}

Incidence-Guided Saturation fixes two choices, four quadrants and saturations from 80\% down to 20\%. To see how much each matters, we varied the number of regions and the spread of the saturations on the grid. The grid was divided into 2 halves, 4 quadrants, 10 rectangles or 25 blocks of equal size. With \(R\) regions and spread \(s\), the region ranked \(r\) by mean incidence received the target saturation \(0.5 + s\{1 - 2(r-1)/(R-1)\}\), so the highest-ranked region was treated at \(50\% + s\) and the lowest at \(50\% - s\). We used \(s = 0.1\), 0.2, 0.3, 0.4 and 0.5, which give highest and lowest saturations of 60/40, 70/30, 80/20, 90/10 and 100/0 percent. Four quadrants at 80/20 is the original design.

The 100/0 setting is a boundary case, included to show whether going beyond 90/10 helps or hurts, not as a design to use. It treats one region fully and leaves another untreated. For quadrants its targets are 100\%, 67\%, 33\% and 0\%, or 25, 17, 8 and 0 treated clusters. Targets were rounded by largest remainders so that every assignment treated exactly 50 of the 100 clusters, and ties in regional mean incidence were broken at random. Each guided variant has an unguided twin that gives the same targets to the regions in random order, so the pair isolates the effect of incidence guidance. The unguided twin of the original design is Saturation Quadrants. With 25 blocks at 60/40, every block has two treated clusters, so both versions are \(2\times2\) Blocking and the pair serves as a null control.

The 38 new variants, the original design, Saturation Quadrants and SRS were run on the full grid of Section~\ref{simulation-design-resampling-and-metrics}, with the same incidence surfaces and outcome noise as the other designs, 25 assignments per surface and both fitted models. We report queen contiguity and the oracle model, pooled over the five incidence configurations and the five values of \(\beta\) within each spillover regime. Each comparison uses the 50 paired configuration-by-surface differences in MSE, with a two-sided Wilcoxon signed-rank test and Holm adjustment within each family. The families compare each design with SRS (40 comparisons), each guided variant with the original design (19) and each guided variant with its twin (20). For the figure, the paired 95\% \(t_{49}\) interval for the mean difference is divided by the comparator's MSE and shifted by 1, which puts it on the ratio scale. Every 100/0 assignment gave full-rank design matrices, and no fit of a new variant was aliased or failed.
```

### N5. New Results section (after Section "Robustness, the SRS benchmark and allocation-specific risk", before "Case study: accuracy, inference and reach")

Current: none (empty slot).

Proposed:

```latex
\section{Sensitivity of incidence guidance to regions and saturation spread}\label{igs-sensitivity-results}

Figure~\ref{fig:igs-sensitivity} shows the MSE of each variant relative to SRS (Section~\ref{igs-variants-methods}). SRS had MSE 0.227 under control-only spillover and 0.043 under both-arms spillover. Under control-only spillover, widening the guided saturations from 60/40 to 90/10 lowered MSE for every number of regions, and at each of these spreads fewer, larger regions had lower MSE. At 90/10 the guided ratios to SRS were 0.55 for halves, 0.59 for quadrants, 0.77 for ten regions and 1.05 for 25 blocks, against 0.67 for the original design. Under both-arms spillover, wider spreads raised MSE for halves and quadrants. Guided quadrants rose from 1.00 times SRS at 60/40 to 1.22 at 90/10, and halves from 0.98 to 1.94. With ten regions or 25 blocks, MSE stayed between 1.00 and 1.22 times that of SRS from 70/30 to 100/0. Some variants traded less between the regimes than the original design. Quadrants at 70/30 and ten regions at 80/20 had ratios of 0.84 and 0.89 under control-only spillover and 1.02 under both-arms spillover, where neither differed detectably from SRS.

Moving from 90/10 to the 100/0 boundary had different effects by number of regions under control-only spillover. MSE rose for halves, from 0.55 to 0.88 times that of SRS for the guided version, and the bias of \(\hat\beta\) rose from 0.086 to 0.161. Quadrants changed little (0.59 to 0.59), while MSE fell for ten regions (0.77 to 0.68) and 25 blocks (1.05 to 0.86). At 100/0, both 25-block versions had detectably lower MSE than SRS for the first time. Under both-arms spillover, MSE rose at 100/0 for every number of regions in both versions. The largest rise was for halves, whose MSE reached 8.13 times that of SRS. When one half is fully treated, \(\mathbf{W}\mathbf{Z}\) differs from \(\mathbf{Z}\) only at the 20 clusters along the dividing line, and the two columns have correlation 0.941 under queen contiguity. The model then has full rank but separates treatment from spillover poorly. At 100/0, coverage of the \(t\)-based interval was between 0.932 and 0.953, lowest for halves under control-only spillover. The changes from 90/10 to 100/0 are descriptive, because the tests compare each design with SRS, the original design or its twin.

Incidence guidance changed where treatment went more than how accurate the estimate was. At the original 80/20 quadrants, guided and unguided MSE were close (ratios 1.00 under control-only and 1.04 under both-arms spillover). No guided variant had detectably lower MSE than its unguided twin, and every detectable difference favored the unguided version. At 100/0 under control-only spillover, guided MSE was 1.04 to 1.10 times that of the twin, and all four differences remained after Holm adjustment. In the targeting check (queen contiguity, \(\rho = 0\), \(\psi = 0.5\), control-only spillover, averaged over the five incidence configurations), the original design treated 52.7\% of total incidence. The share grew with the spread and the number of regions, to 57.7\% for 25 blocks at 90/10 and 53.6--60.1\% at 100/0, while the unguided versions and SRS treated about 50\%.

Rook contiguity gave the same directions from 90/10 to 100/0. Its halves at 100/0 had slightly higher MSE than SRS under control-only spillover and about 17.3 times that of SRS under both-arms spillover. With the non-oracle model, MSE fell from 90/10 to 100/0 for quadrants under control-only spillover with both neighbor types, and the both-arms penalty for halves remained, with coverage as low as 0.607 under rook contiguity. Adding the 100/0 variants enlarged the Holm families, and 3 of the 756 earlier comparisons rose above 0.05, all in rook sensitivity analyses. No queen oracle conclusion changed. These results cover equal-size partitions of the grid with linearly spaced targets. The NC map was not part of this study.

\begin{figure}[!htb]
\centering
\includegraphics[width=0.9\linewidth]{../figures/fig_igs_re_vs_srs.pdf}
\caption{MSE of the variants of incidence guidance relative to SRS on the grid under queen contiguity with the oracle model, pooled over the five incidence configurations and five values of \(\beta\). The horizontal axis gives the highest and lowest target saturations in percent, and colors give the number of equal-size regions. The left panels give the targets in order of regional mean incidence, and the right panels give the same targets in random order. Rows are the spillover regimes and have different vertical scales. Bars are paired 95\% \(t_{49}\) intervals on the ratio scale, and the line at 1 marks SRS. Circles mark Incidence-Guided Saturation (left) and Saturation Quadrants (right).}
\label{fig:igs-sensitivity}
\end{figure}
```

Checks on N5 claims (from `igs_summary.txt`, queen oracle):
- "1.00 and 1.22 ... from 70/30 to 100/0" for ten regions and 25 blocks: guided 1.003–1.225,
  unguided 1.001–1.185 (25 blocks at 60/40 is 1.46, the 2×2 Blocking null pair, so it is excluded
  by the range).
- "Quadrants at 70/30 ... 0.84 ... 1.02" and "ten regions at 80/20 ... 0.89 ... 1.02": guided
  0.837/1.018 (Holm 1.000) and 0.887/1.020 (Holm 1.000). This is a descriptive observation, not a
  recommendation (no single recommended design); **D5 below**.
- Original design vs its twin: vsTwin 0.998 (control-only) and 1.039 (both-arms), printed as 1.00
  and 1.04. The summary file gives no twin p-value for the original; the findings say every
  detectable twin difference favored the unguided design, and the text claims no more.

| ID | Decision | Recommendation and reason |
|---|---|---|
| D5 | Keep the "Some variants traded less between the regimes" sentences in N5 | **Keep.** It is what a reviewer asks first after seeing the figure, and the numbers say it; it stays descriptive. Drop it if you would rather not invite "then why not 70/30?". |

### N6. `paper/tools/copy_figures.sh`: one new map line (not manuscript text)

`"results/igs_sensitivity_2026-11/full/summary/fig_igs_re_vs_srs.pdf|fig_igs_re_vs_srs.pdf|code/30_igs_sensitivity_summary.R"`

### C2. CTJ, Discussion, after the High Incidence Focus paragraph's last sentence

Current: "Saturation keeps both arms in every region while varying the local share treated, which helps separate the direct effect from neighboring exposure."

Proposed: "Saturation keeps both arms in every region while varying the local share treated, which helps separate the direct effect from neighboring exposure. On the grid, widening the quadrants' saturation spread beyond 80/20 lowered control-only MSE but raised both-arms MSE (SI Section~\ref{SI-igs-sensitivity-results})."

(+21 words. Quadrant ratios: control-only 0.67 → 0.59 → 0.59, both-arms 1.09 → 1.22 → 1.63.)
Chapter 3 needs the same sentence (2i).

## Chapter 3 (for 2i)

N1 with its figure, N4, N5 with its figure, and C1/C2 go into the Chapter 3 appendix and text.
The chapter's figure copies must meet the thesis rule (8 pt, 6.5 in, embedded fonts): the IGS
figure is already 6.5 × 5.5 in with 8 pt labels; the rejection figure needs checking.

## Decisions (author)

(pending)
