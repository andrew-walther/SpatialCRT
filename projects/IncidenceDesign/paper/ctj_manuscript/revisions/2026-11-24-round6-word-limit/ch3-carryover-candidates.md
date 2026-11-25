# Round 6: Chapter 3 carry-over candidates (input to Phase 1e)

Author: Andrew Walther. Written 2026-11-25 after 1d. Candidates only; nothing in
`paper/thesis_chapter/Thesis_Chapter.qmd` was edited. Each needs a current → proposed check
against the chapter's own wording in 1e. Length cuts (the 1b T/N/C items) do not carry over:
the chapter has no word limit and keeps the fuller text.

## Likely needed (content or citation)

| CTJ ID | Change | Chapter 3 target |
|---|---|---|
| A1 | `\cite{lee_asymptotic_2004}` after "the Wald interval" | Ch3 line ~148 ("The 95\% Wald interval uses $q = z_{0.975} = 1.96$…"). Needs the Lee entry in `bios-dissertation/prelim/references.bib` first, or the sync hook fails |
| TAU-1 | Fuller first definition of τ (holding treatment and spillover fixed; τ = 0 in the case study, τ > 0 on the grid) | Ch3 model paragraph (~line 78) and the τ paragraph (~line 162) |
| N5 | τ wording: "The investigator specifies τ…", one-SD shift "before spatial spread", equal to β = 1 | Ch3 ~line 162 |
| R2 | Configuration did not change which design had the lowest control-only MSE | Ch3 grid results (check whether already stated; S-03 was marked done for Ch3) |
| T6b | Rounding rule for spatial blocking (nearest even; one-cluster group never treated) | Ch3 Appendix B case-study adaptations (~lines 414–416, 799) |
| HW-03b | Abstract: define Incidence-Guided Saturation at first mention | Ch3 has no abstract; check the chapter's first results mention |

## Style edits on sentences the two documents share (check each)

- "Therefore,"/"However," at sentence start: TH-1 (β in units of the noise SD), TH-2 (grid ranking
  need not carry over), HW-29 (low MSE alone), HW-33 (two changes made together), T7 (worst decile).
- Semicolons and colons removed: HW-05 (Checkerboard definition; Ch3 wording differs, "Chapter 2"),
  HW-10 (Checkerboard identifiability), HW-16, HW-17 (configuration list), HW-19 (NC data only in
  the case study), HW-26 (two-sided test), HW-34 (conditional guidance), C3 (direct effect).
- ", so" joins: HW-07, HW-08, HW-12, HW-21, HW-37.
- Parenthetical results: HW-27 (ratios), HW-28 (β = 2, 3; keep "nearly equal"), HW-30.
- Table 1 / `tab:designs` cells: TC-1, TC-2 (Ch3's table was already reworded without ", so";
  the SI copy still differs, Phase 2d).
- Captions: HW-40 to HW-51 where Ch3 uses the same figure captions.

## Not carried over (by decision)

- C-16: Ch3 keeps the HC1 wording ("the divisor … of the HC1 covariance estimator").
- T11's dropped caveat ("spread … differs somewhat by mechanism") stays in Ch3.
