# Incidence-Guided Saturation sensitivity study: independent review

Reviewed 2026-11-24. Scope: `docs/plans/2026-11-24-igs-sensitivity-spec.md`,
`code/29_igs_sensitivity.R`, `code/29_igs_variant_designs.R`,
`code/30_igs_sensitivity_summary.R`, the outputs in `results/igs_sensitivity_2026-11/full/`
and the report `docs/findings/2026-11-24-igs-sensitivity-findings.md`. Nothing was rerun
except small R checks against the saved outputs (scratch scripts, not kept in the repository).

## Verdict

**Sound with fixes.** The code does what the spec says, the seeds and pairing are right, the
reproduction check is valid, and the summary computes the contrasts the report describes. I
re-derived about 40 numbers from the saved outputs and they match. One variant is identical in
distribution to 2×2 Blocking and to its own twin, which the report describes as "almost". Two
claims need fixing before anything reaches the SI or the CTJ: the guided-vs-twin bound ("within
about 8%", "at most 8%", "at most a few percent") is wrong, and the closing sentence about the
North Carolina map is contradicted by stored results. Several smaller overstatements and some
design points a committee will raise follow.

## What was checked and holds

- **Partitions.** `region_index()` gives equal regions for every R (50/50; 4×25; 10×10;
  25×4). Shapes are as stated: R = 2 splits on x only (two 5×10 halves, west and east), R = 10
  is 2 wide by 5 tall, R = 25 matches the 2×2 Blocking blocks of Design 5 exactly.
- **Guided direction.** `rank_r <- R + 1 - random_tie_rank(means)` gives the highest-mean
  region rank 1 and `sats[1] = 0.5 + s`, the highest saturation. This is the same map as
  Design 8's `sat_levels[5 - random_tie_rank(quad_means)]`, with ties at random per draw. For
  R = 4, s = 0.3 the variant reproduces Design 8's distribution (counts 20/15/10/5).
- **Unguided twin.** `sample.int(R)` is a uniform random permutation of the same saturation
  set per draw, as Design 3 does with `sample()`.
- **Counts.** `treated_counts()` always returns 50 (checked in the script and by me over 200
  draws per variant). Ties among equal remainders use `runif()` under the seeded RNG, so they
  are fair and reproducible. Ties occur only at R = 4 (s = 0.2, 0.4) and R = 25 (s = 0.3), where
  the highest- and lowest-ranked regions share a remainder of 0.5; for example R = 4, s = 0.2
  gives 18/14/11/7 or 17/14/11/8.
- **(R = 4, s = 0.3)** is correctly dropped from VARIANTS and covered by Designs 8 and 3.
- **Seeds.** `seed_key()` formats variant IDs as "101.00" etc., so variant Z keys never equal
  production keys. The noise E is keyed by block only and drawn before the design loop, and
  X by (mode, ρ_X, k), so every design in a block shares X and E. The pairing claim holds.
  Z keys do not include the effect size, so the five effect sizes share Z and E (see S6).
- **Overrides.** In 29, `get_designs`, `design_ids` and `summarize_fits` are redefined in the
  global environment after 05 is sourced, and 05's `draw_assignments()`/`run_unit()` look them up
  there at call time, including in forked workers. In 30 the order (10 → 05 → 29 variants →
  override) leaves the wrapper in force; nothing later re-sources 03.
- **t + df.** Same multipliers as `code/26` (df 95 oracle, 96 non-oracle, `qt(0.975, df) ·
  sqrt(N/df)`), and `pick()` gives each estimator its own df. I also checked the rerun's
  `Coverage_Tdf` and `Power_Tdf` for Designs 3, 8 and 9 against
  `results/grid_tdf_2026-10/grid_tdf_scenarios.rds`: 9,600 rows matched, zero difference.
- **Reproduction check.** It is valid: 9,600 rows, zero difference in seven metrics, and it
  stops on mismatch. No fit failed, none was aliased, no warnings, no rank-deficient [1, Z, WZ].
- **Summary logic.** `ratio = 1 + D̄ / MSE_comp` equals MSE_d / MSE_comp exactly here,
  because the slices are balanced and no fit failed. Pooling over effect sizes (mean of
  scenario MSEs) and the 50 (configuration × surface) units are as in `code/15`. Twin matching
  is one-to-one (the `vapply(..., character(1))` would fail otherwise). Burden share is
  Σ X_j z_j / Σ X_j, as in `code/22`.

### Numbers verified (queen, oracle unless stated)

All 32 cells of the main table and their asterisks; SRS MSE 0.227 and 0.043; Saturation
Quadrants 0.668 and 1.047; 16% to 45% (R = 2 and 4, s ≥ 0.2); 0.555, 0.551, 0.588; 2×2
Blocking 2.579 in the production run; 1.940; original 0.667/1.088 and halves at 70/30
0.676/1.051; twin ranges 0.984 to 1.076 and 0.955 to 1.095; the Holm-surviving twin
differences named in the text; burden shares 0.507, 0.577, 0.500 (R = 25, s = 0.1), 0.527,
and 0.522 to 0.530 by configuration; coverage 0.940 to 0.953; rook original 0.810; non-oracle
0.642 and 2.364; 105,600 rows; 78.5 minutes.

## Findings

1. **Error. The R = 25, s = 0.1 variant is identical in distribution to 2×2 Blocking, and the
   guided version is identical to its unguided twin.** Largest-remainder rounding of
   4 × sat_r (2.4 down to 1.6) gives exactly two treated in every one of the 25 blocks, every
   draw. So guidance has nothing to act on, and the design is Design 5. The report says
   "almost everywhere", "close to that of 2×2 Blocking" (2.674 vs 2.579) and "rounding treats two
   in almost every block". The three values (2.674 guided, 2.580 unguided, 2.579 Design 5) are
   the same design under three seeds, so their 4% spread is Monte Carlo noise.
   *Fix:* say the variant is 2×2 Blocking exactly; drop it from the twin family (or keep it as a
   labelled null control, see S2); mark it in the table and figure (its "60/40" label is
   nominal; the realized set is 50/50).

2. **Error (report and CTJ draft). The guided-vs-twin bound is understated.** The summary says
   "within about 8%", the CTJ draft sentence says "changed MSE by at most 8%", and "What this
   changes" says "at most a few percent". Under queen/oracle the largest twin ratio is 1.095
   (R = 25, s = 0.2, both arms), and the body itself says 0.96 to 1.10. Across the rook and
   non-oracle checks it reaches 1.111 (rook, oracle, control-only, R = 25, s = 0.4) and 1.137
   (rook, non-oracle, both arms). Those are Holm-significant penalties, not noise.
   *Fix:* "within about 10% under queen contiguity (up to 14% in the rook and non-oracle
   checks), and every Holm-significant difference was against guidance." Change the CTJ
   sentence to match, or say "at most 10%" and name the slice.

3. **Error (unsupported claim). "The North Carolina map has larger regional differences."** The
   stored NC reference in `results/grid_burden_share/grid_burden_share_log.txt` records an
   Incidence-Guided Saturation rate share of 0.511 against 0.504 for plain saturation on NC, a
   gain of 0.007. On the grid it is 0.527 against 0.500, a gain of 0.027. The sentence implies
   larger gains on NC, which the evidence contradicts, and nothing in this study measured NC.
   AGENTS.md also frames NC as "a small gain".
   *Fix:* delete the sentence, or replace it with the stored NC numbers.

4. **Overstatement. The coverage bullet reads as if it covers the non-oracle fit.** The
   summary says t + df coverage "was 94% to 95% for every version. The pattern held under rook
   contiguity and for the model fit without the spillover term." Non-oracle t + df coverage was
   0.526 to 0.869 under control-only spillover and 0.851 to 0.952 under both-arms spillover.
   Rook/oracle is fine (0.945 to 0.953).
   *Fix:* "Coverage of the t + df interval with the oracle fit was 94% to 95% for every version,
   under queen and rook. The MSE pattern held for the model without the spillover term, whose
   coverage was low, as in the main study."

5. **Overstatement. "Flat sets were close to SRS" under both-arms spillover (summary
   bullet 2).** True for R = 2, 4 and 10 (0.977, 1.005, 1.037). At R = 25 the flat set is 1.463,
   which is just 2×2 Blocking (finding 1). The body sentence carries the qualifier; the summary
   drops it.

6. **Overstatement. "Under rook contiguity the ordering was the same" and "Rook gave the same
   picture."** The rank correlation of the 32 ratios with queen is 0.993 under control-only
   spillover but 0.903 under both-arms spillover. In the twin contrast rook had four
   Holm-significant guided penalties under control-only spillover (queen had one), up to 1.111.
   The direction holds and the magnitude differs.
   *Fix:* "a similar ordering (rank correlation 0.99 and 0.90 by regime)", and give rook's
   largest twin ratio.

7. **Overstatement (interpretation as a finding). "Accuracy came from the saturation
   structure"** (summary bullet 4). The twin ratios support "guidance did not drive accuracy".
   Which part of the structure matters (number of regions, shape, contrast) cannot be separated
   by this design (S1). Suggest "Accuracy followed the regions and saturation set, not the
   guidance."

8. **Minor accuracy points.**
   - "reaching 0.087 for two halves at 90/10": in a paragraph about guided versions, but 0.087 is
     the unguided twin; the guided version is 0.086.
   - "13.2 million fits": 13.2 million simulated datasets, each fit twice (oracle and non-oracle),
     so 26.4 million model fits.
   - "Holm-adjusted across the 32 designs": correct for the code, but the spec says 30 variants
     against SRS and 15 guided against the original. The code uses 32, 15 and a 16-pair twin
     family. Record the deviation in the spec.
   - The pooled-over-regimes slices are computed but never used in the report; harmless.

## Suggestions (design and presentation)

S1. **R is confounded with the size of the saturation contrast, with region shape and with
rounding.** A committee member will ask whether "fewer, larger regions" is a region effect or a
contrast effect.
- For a linear-in-rank set the variance of the saturations is s²(R + 1)/(3(R − 1)). At a fixed s
  that is s² for R = 2, 0.56 s² for R = 4, 0.41 s² for R = 10 and 0.36 s² for R = 25. So
  "spread s" means a larger between-region contrast for two halves than for 25 blocks.
- Shape changes with R: halves and 2×5 rectangles are elongated (and R = 2 splits only
  west/east), while quadrants and 2×2 blocks are square. Under queen contiguity the share of a
  cluster's neighbours that lie in its own region, which is what spillover sees, depends on shape
  as well as size.
- Rounding at R = 25 changes the realized set: s = 0.1 → all 2/4 (finding 1), s = 0.2 → five
  blocks at 3, fifteen at 2, five at 1 (0.75/0.5/0.25, not linear 0.7 to 0.3).

*Fix:* report the realized saturation SD for each variant (or a mechanism summary such as mean
corr(Z, WZ) or the mean share of treated neighbours of controls), and add a panel of the ratio
against it. If the points collapse onto one curve, the "region size" message becomes a
"clustering of treatment" message, which is more defensible. At least state the confounding in
the Limitations section.

S2. **Use the identical pair as an empirical null.** Guided vs unguided at R = 25, s = 0.1 has a
ratio of 1.036 [1.007, 1.065], unadjusted Wilcoxon p = 0.031 (queen, oracle, control-only)
although the two designs are the same. Across the eight slices it ranges from 0.985 to 1.036.
So twin differences of 3 to 4% are within Monte Carlo noise at this precision. That fits the
Holm result, but it should be said, because the report's "more often slightly worse than better"
count (13 of 16) includes this null pair and several ratios of 1.02 to 1.03.

S3. **The guidance penalty is variance, which supports the report's interpretation.** For every
twin pair the squared-SD ratio tracks the MSE ratio (for example 1.091 vs 1.095 for R = 25,
s = 0.2, both arms), and bias differences are small. That fits variance inflation from Z being
correlated with X in the oracle model. The report can say "consistent with" and cite this rather
than "not something this study tested".

S4. **Asterisks and intervals disagree in places.** The table's asterisk is Holm-adjusted
Wilcoxon, while the figure's bars are unadjusted t₄₉ intervals. For example, R = 10, s = 0.2,
control-only has interval [1.006, 1.076], which excludes 1, but carries an asterisk; R = 2,
s = 0.2, both arms is similar. Say in the caption that the bars are unadjusted and the asterisks
are Holm-adjusted Wilcoxon tests.

S5. **The second named primary contrast (each guided variant against the original) is not
reported.** The CSV has it (`vsOrig_*`), but the report gives only the point ratio `RE_orig`
implicitly through the SRS table. Either report it (one column or a sentence) or demote it to
secondary in the spec.

S6. **Pooling over effect sizes changes the size of the gains, not the ordering.** The five
effect sizes share Z and E, so they are not replicates. The 32-design ranking is stable across
them (rank correlation with the pooled ranking 0.992 to 1.000). The control-only gain shrinks as
the effect size grows, though: the original design's ratio runs from 0.618 at β = 0.8 to 0.718
at β = 3, and two halves at 90/10 from 0.417 to 0.699. The "16% to 45%" headline therefore
depends on the pooling. Add one sentence with these ranges to the pooling limitation.

S7. **Reproduction check robustness.** `max(abs(a - b)[...], na.rm = TRUE)` would hide a row
where one run is NA and the other finite. It cannot matter here (no failed fits), but
`identical()` on the matched columns, plus the t + df columns against `results/grid_tdf_2026-10/`
(which I confirmed match), would make the check airtight.

S8. **"What this changes".** Bullet 1 follows from the evidence, subject to the S1 caveat that
"large regions" and "contrasting saturations" are not separated. Bullet 2 ("guidance matches the
accuracy of unguided saturation") holds for the original R = 4 design (twin ratios 0.998 and
1.039) and broadly. It needs the corrected bound from finding 2, plus the note that the
measurable cost of guidance is concentrated in many small regions (R = 10, 25) at wide spreads.
