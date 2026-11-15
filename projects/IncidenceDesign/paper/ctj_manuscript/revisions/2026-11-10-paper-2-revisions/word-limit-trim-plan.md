# CTJ word limit: trim plan (not yet applied)

Author: Andrew Walther. Prepared 2026-11-15, after the first revision round (SpatialCRT `43b6f99`).
Trimming is deferred until after the proposal (Decision C); this file is the starting list.

## Where the text stands

Main text after this round, counted from `CTJ_Manuscript.tex` (Introduction to Conclusions, no
abstract, equations or captions; each inline math expression counts as one word). This count
is about 100 words below a texcount count, which counts inside math. Recount with the journal's
method before submission; the title page now shows a placeholder.

| Section | Words |
|---|---|
| Introduction | 435 |
| Methods: outcome model, estimation and allocation risk | 550 |
| Methods: treatment assignment designs | 169 |
| Methods: simulation study | 489 |
| Methods: case study | 671 |
| Results: simulation study | 371 |
| Results: case study | 578 |
| Discussion (four subsections) | 575 |
| **Total** | **3,838** |

Target: 3,500, and about 3,400 by this count to leave room for the journal's method. That
means cutting about 450 words. The candidates below add up to about 570, so some can be left out.
Moves go to the SI with a pointer. Nothing is deleted outright: a cut sentence that has no SI
home goes to `paper/thesis_chapter/notes/trimmed-material.md` with date and reason.

## Candidates, largest first

| # | Where (C line) | What | Action | Saves |
|---|---|---|---|---|
| T1 | 173 Simulation study | The written-out incidence models (spatial Φ{…}, Poisson log-rate, base-rate sentence, rank-scaling) | Keep one sentence naming the five configurations and the (0, 1) scale; move the formulas to SI Section 2 with a pointer. **Conflicts with FCL-032**, where Dr. Lin asked for the model for X; a compromise keeps only the spatial formula (saves about 70) | 110 |
| T2 | 193, 213 Case-study replication | The precision targets appear in both 193 and 213; the 400-outcome and split-half details | Delete 213's restatement; move the split-half check and 400-outcome detail to SI Section 14, keeping one clause | 90 |
| T3 | 232–238 Implications | Four paragraphs; 236 (High Incidence Focus) and 238 (judging the regime) partly repeat Results and 232 | Merge 236 and 238 into one or two sentences | 60 |
| T4 | 170 Designs | Restates the categories already in Table 1 | Keep the last sentence only ("These restrictions change how closely ...") | 45 |
| T5 | 180 NC data | IRB protocol (repeated in Ethical considerations) and the heart-failure correction detail | Keep the source, years and totals; move the correction detail to SI Section 5 | 45 |
| T6 | 187 NC adaptations | Checkerboard rule and the rounding rule "(1.5 to 2, 0.5 to 0) ..." duplicate Table 1 | Point to Table 1 for the Checkerboard rule; move the rounding rule to the SI | 40 |
| T7 | 222 Case-study allocation risk | Three sentences on worst-decile uncertainty | Shorten to one sentence | 40 |
| T8 | 147 Introduction | "Incidence can direct resources ... None of these goals ensures ..." | Merge into one sentence | 30 |
| T9 | 226 vs 244 | High Incidence Focus 22–24% reach is stated in both | Keep it in Results, refer to it in 244 | 30 |
| T10 | 204 Results | The Isolation Buffer explanation of identical regimes | One clause | 25 |
| T11 | 173 | "All configurations put incidence on the (0, 1) scale ..." repeats the rank-scaling sentence and M-11's "Incidence lies between 0 and 1" | Delete, or merge with M-11 | 25 |
| T12 | 158 | The `lagsarlm` agreement check | Move to the SI (it is already there) | 22 |
| T13 | 224 vs 234 | "22% with 58 colleges to 5% with 29" appears twice | Keep it in 234 (Discussion) only | 20 |
| T14 | 145 | "We do not compare numerical results between the two studies." | Delete or merge | 8 |

## Other things to check when trimming

- The abstract is about 400 words. Check the current Clinical Trials limit for structured
  abstracts before submission (not checked here).
- C 193 ends with two nested parenthetical pointers ("...; the grid power in SI Section 4 ...;
  SI Section 5.5"). Splitting it into its own sentence costs a few words but reads better.
- Captions are not counted, but moving detail into captions should not be used to dodge the limit.

## Suggested order

T2, T4, T6, T9, T11, T12, T13, T14 first: they remove repetition and lose no content (about 280
words). Then T3, T5, T7, T8, T10 (about 200). Leave T1 for last, since it reverses a reviewer
request.
