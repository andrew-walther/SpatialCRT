# Review of the CTJ round 6 word-limit cut list

Written 2026-11-24. Author: Andrew Walther. Read-only review of
`paper/ctj_manuscript/revisions/2026-11-24-round6-word-limit/signoff.md` (from `cuts.py`) against
`CTJ_Manuscript.tex` lines 140-248, the SI source, `proposed-edits.md`, the 2026-11-14 meeting
notes and `docs/reviews/2026-11-24-ctj-revision-audit.md`. Nothing was edited or applied. Word
savings were measured with the same counter `cuts.py` uses (`count()` from `ctj_wordcount.py`);
the measuring script lived in the session scratchpad.

SI section numbers were confirmed from `supplementary_info/.build/Supplementary_Information.aux`:
1.2 heterogeneous incidence; 3 estimation; 4 simulation design (4.1 comparisons); 5.1 observed
inputs; 5.2 geographic adaptations; 5.3 outcome model and Monte Carlo evaluation; 5.4
reproducibility and precision; 5.5 partial participation; 12.3 allocation-specific risk; 14
case-study allocation downside; 15.2 interval calibration; 15.4 case-study sensitivity analyses.

## Summary

- 14 items are fine as written, 9 are fine with a small edit, and 3 should be reconsidered
  (T8, T13, N2). N4 and T1 stay rejected, as the list advises.
- Two notes in the list are wrong about the SI. N6's content is **already in SI 5.3**, so nothing
  needs adding. N12's content is in **SI Section 14, not 15.4**, and it is already there.
- R2 is verified against the grid t + df results. R3 is verified, but SI 15.4 compares **three**
  forms (proportional, count, binary), and the sentence needs a citation to the earlier study.
- With the edits below, the recommended set measures **3,476** words (24 under the limit). The
  margin is thin, since a journal count may treat inline math or citations differently. Adding N8
  and a few of the new candidates C1-C8 gives about **3,380**.

## Verdict per item

| ID | Verdict | Saves (as listed → as edited) | Notes |
|---|---|---|---|
| T14 | fine | 10 | The sentence was part of the approved I-05 paragraph. The AGENTS rule (no numerical cross-project comparison) is met by making none. Low risk. |
| T8 | **reconsider (reject)** | 1 | Saves one word by dropping "historical" (the burden the case study uses is historical, as the previous sentence says) and by joining two sentences into one of about 40 words, against meeting item W-02 (break long sentences). Not worth it. |
| T12 | fine | 18 | SI Section 3 has the validation (9.7×10⁻⁸ in β̂ on 5,120 fits). Optional: add "(SI Section 3)" after the ML citation (+3). |
| N2 | **reconsider: fine with edit** | 10 → 9 | "so that outcome noise did not drive the ranking" overstates. Results say outcome noise accounts for a median 63% of the variance in estimated risk (L210), and that most both-arms variance is Monte Carlo error (L222). Use: "We estimated each $m(a)$ by averaging squared errors over many simulated outcomes, to limit the influence of outcome noise on the ranking of assignments." |
| T4 | fine | 50 | Table 1 gives each rule and category. "Saturation" stays defined at L147 and in Discussion. The result is a one-sentence paragraph; it could join the end of L165's paragraph (no word change). |
| T11 | fine with edit | 36 → 27 | The deleted sentence carries the answer to meeting item S-04 (why incidence is on (0,1): so one τ means the same shift in every configuration). M-11's sentence explains what τ = 1 means, not why the scale is shared. Keep the reason in the rank-scaling sentence: "... to put them on the same $(0,1)$ scale as the other configurations, on which one value of $\tau$ gives the same shift." |
| N9 | fine with edit | 13 → 11 | SI Section 4 (and 4.1) explains the 50 independent units. Keep "between designs" so it is clear what is paired: "MCSEs and paired comparisons between designs used the 50 independent configuration-by-surface units (five configurations by ten surfaces)." |
| N3 | fine | 21 | MCSE stays defined at first use (commit 71d0ce3). No "so" join (FCL-041/M-16 satisfied). |
| T5 | fine | 28 | IRB protocol 14-2036 is in Ethical considerations (L267). The heart-failure correction is in SI 5.1 (line 637). "corrected" in the previous sentence is then unexplained in the body; optional pointer "(SI Section 5.1)" after "prepared by Habib" (+3). |
| T6 | fine | 26 | Table 1's NC cell has the identical rule (E-02), so FCL-039 stays answered in the table. SI 5.2 has the full rule. |
| T6b | fine (SI addition needed) | 18 | SI 5.2 (line 733) says only "Half of each group, rounded". The round-half-to-even rule (1.5 to 2, 0.5 to 0) is not in the SI and should be added there when applied. Pointer to 5.2 is correct. |
| N7 | fine with edit | 11 → 11 | "its share of the working-age population, which is distinct from" lets "which" attach to "share". Keep "a population distinct from" (same count): "We report each design's number of treated clusters and its share of the working-age population, a population distinct from ..." |
| N6 | fine (note is wrong) | 24 | The list says "not there yet; add it when applying". SI 5.3 already says it ("Settings with the same set of possible assignments shared their simulations ... Checkerboard, Saturation Quadrants, Isolation Buffer and SRS reuse the same simulations across years"). No SI change needed. |
| T2 | fine with edit | 24 → 22 | "half of them" is ambiguous (assignments or outcomes). Use "... by selecting it with half of these outcomes and evaluating it with the other half." 400 outcomes are in SI 5.3 and the split-half check is in SI 5.4. |
| T10 | fine with edit | 19 → 9 | Table 2 shows Isolation Buffer at 0.161 and 0.159 in the two regimes; without the Monte Carlo clause a reader may take the difference as real. Use: "For Isolation Buffer the two regimes define the same model, because no treated cluster has a treated neighbor, and its small differences between regimes are Monte Carlo error." |
| T2b | fine | 19 | The targets stay in Methods L193, which no item removes. |
| T7 | fine (optional edit) | 27 → 26 | Retains a ", so" join (W-03). Optional: "(SI Section 14), and small both-arms differences therefore do not show better downside protection." SI 14 pointer is correct. |
| T13 | **reconsider (swap)** | 20 → 18 | Cutting the Results clause leaves the 22% → 5% numbers first appearing in the Discussion (L234), and L232's guidance "when ... most colleges take part" would then rest on no Results statement. Keep the Results clause and instead delete the Discussion sentence "In the case study, its advantage over SRS fell from 22\% with 58 colleges to 5\% with 29." (saves 18). R-09's point (FCL-059, "not best in every setting") still stands on the two grid instances. If you would rather keep R-09 verbatim, reject T13 outright. |
| N12 | fine with edit (pointer is wrong) | 14 → 11 | The content is in SI **Section 14** (Supplementary line ~2121: "Alternative regional summaries change the regional ordering in some years ... We retain the prespecified mean-rank rule"), with the orderings in SI 5.2's Table. SI 15.4 does not discuss it. No move needed. Because Methods L191 still lists this sensitivity analysis and Results would no longer report it, add a pointer there: "... or pooled regional rates instead of mean rank (SI Section 14)." (+3). |
| T3 | fine (optional edit) | 26 | Retains ", so it does not support" (W-03); optional "and therefore does not support" (same count). The FCL-058 reply (reply only, pointing to C 238) should now point to the merged paragraph. Drops "is consistent with its control-only advantage", a small loss. |
| N13 | fine with edit | 36 → 35 | After the deletion, "Before using them to choose an allocation" has no antecedent (the previous sentence is about consequences). Change to "Before using these results to choose ...". Limitations (L241) keeps the equal-cost point, and Results (L226) keeps the reach figures. |
| R2 | fine (verified) | -13 | See the R2 check below. |
| R3 | fine with edit | -27 → -31 | Verified (below). Edit: cite the earlier study (CTJ must cite the BMC paper, not say "Project 1") and name the three forms the SI compares: "Spillover here is proportional to the share of treated neighbors, whereas our earlier study used binary (on-or-off) spillover;\cite{walther_spatialcrt_2026} SI Section 15.4 compares proportional, count and binary forms in the case study." |
| N5 | fine with edit (optional) | 27 → about 24 | Keeps D-2, M-11 and Q-02 content; "The investigator sets τ by whether incidence is expected to shift the outcome" covers D-2's "investigators choose it for their own outcome". Two fixes: "On the grid $\tau=1$" reads clipped, and the "because" clause needs a closing comma: "On the grid we set $\tau=1$, and because incidence lies between 0 and 1, this raises ...". If T11 is applied as edited, N5 still works with it. |
| N8 | fine (could be recommended) | 5 | Keeps FCL-002/003/010 content; "and vice versa" is clear. |
| N4 | reject (agree) | 14 | Author's own M-14 wording answering meeting item S-01. |
| T1 | reject (agree) | 37 | Removes the Poisson model, the 100,000 population and the rank-scaling formula that answer FCL-032, -033, -034 and -036 (S-02). |

Style check on all proposed text: no new colons, no "fitted" used as a verb, American spelling
throughout. The only ", so" joins are ones carried over from the current text (T3, T7).

## R2 verification

`results/grid_tdf_2026-10/grid_tdf_scenarios.rds`, oracle, queen, `True_Tau == 1`, Designs 1, 2, 4,
5, 6, 8, 9; mean MSE over the 16 (ρ, ψ) settings in each cell (1,120 rows).

Lowest three by configuration:

| Configuration | Control only | Both arms |
|---|---|---|
| iid | IGS 0.139, IB 0.157, SRS 0.203 | SRS 0.0419, BQ 0.0423, IGS 0.0456 |
| Poisson ρ_X = 0.2 | IGS 0.134, IB 0.162, BQ 0.216 | BQ 0.0411, SRS 0.0427, IGS 0.0465 |
| spatial ρ_X = 0.2 | IGS 0.129, IB 0.159, SRS 0.217 | BQ 0.0419, SRS 0.0420, IGS 0.0453 |
| Poisson ρ_X = 0.5 | IGS 0.137, IB 0.160, SRS 0.217 | SRS 0.0430, BQ 0.0440, IGS 0.0472 |
| spatial ρ_X = 0.5 | IGS 0.138, IB 0.160, SRS 0.216 | BQ 0.0414, SRS 0.0422, IGS 0.0448 |

(IGS Incidence-Guided Saturation, IB Isolation Buffer, BQ Balanced Quartiles.)

- The claim holds. Under control-only spillover Incidence-Guided Saturation had the lowest MSE in
  all five configurations (0.129-0.139), with Isolation Buffer second in all five. Under both-arms
  spillover SRS and Balanced Quartiles traded first place within 0.0016, so R2 correctly limits
  the claim to control-only spillover.
- The existing "<10% except High Incidence Focus" sentence also holds: the largest-to-smallest
  configuration ratio minus 1 was 19.2% for High Incidence Focus under control-only spillover and
  at most 8.9% for every other design and regime (SRS 8.9%, Incidence-Guided Saturation 8.2%).

## R3 verification

SI 15.4 ("Sensitivity analyses" of Section 15, `\label{nc-sensitivity-results}`, aux 15.4) has a
"Spillover form" paragraph and `tab:nc-forms` comparing proportional, count and binary forms on the
58 NC colleges. SI 5.6 defines the binary form as an indicator of any treated queen neighbor (times
1 − Z under control-only spillover). The grid has only the proportional form, as the meeting
decision (Q-01) notes. The pointer is correct; the edit above adds the count form and the
citation.

## New candidates not on the list

Exact current text → proposed, words saved (measured).

| ID | Where | Current | Proposed | Saves | Notes |
|---|---|---|---|---|---|
| C1 | L185 | "The fit included an intercept, treatment and the true spillover term but not incidence." | *(delete)* | 14 | L158 already says the fit includes incidence "in the grid study". The τ = 1 sensitivity (L191) still says incidence was added. |
| C2 | L156 | ", which multiplies the two vectors cluster by cluster" | *(delete)* | 8 | Restates "elementwise product". Keep if it was added for clinical readers. |
| C3 | L156 | "The intervention effect $\beta$ is a direct effect: the effect of a cluster's own treatment" | "It is a direct effect, the effect of a cluster's own treatment" | 2 | Also removes a prose colon. |
| C4 | L165 | "SRS randomly selects 50 of the 100 clusters for treatment, and on the grid every design treated 50 clusters except Isolation Buffer, whose number treated depended on the neighbor graph." | "On the grid every design treated 50 clusters except Isolation Buffer." | 19 | SRS is defined in the Introduction, Table 1 and the Table 2 caption; Table 1 gives Isolation Buffer's count. |
| C5 | L222 | "Checkerboard and High Incidence Focus each have a single possible assignment for these inputs, so they show no variation across assignments but can still have substantial error." | *(delete)* | 27 | Figure 3's caption says a one-assignment design "appears as a circle only", and Discussion L232 makes the "high error for every assignment" point. Dr. Lin asked to delete the Methods version of this sentence (FCL-029, M-09). Also removes a ", so". |
| C6 | L241 | "The primary model used the true form of spillover, and omitting the spillover term produced biased estimates and poor coverage." | "The primary model used the true form of spillover." | 11 | The omission result is reported in Results L208. |
| C7 | L244 | "Under control-only spillover even all 58 colleges gave less than 80\% power at $\beta=1$, so the expected effect size and outcome variability must be established first." | "The expected effect size and outcome variability must also be established first, because power was low under control-only spillover." | 6 | Repeats L224; removes a ", so". |
| C8 | L226 | "Treatment budgets differed by design and year (Table~\ref{tab:allocation-rules} and SI Table S5). High Incidence Focus treated 29 clusters" | "High Incidence Focus treated 29 clusters (SI Table S5)" | 9 | Methods L187 already says budgets differed (Table 1). |
| C9 | L147 | "In this paper, we also use a case study" | "We also use a case study" | 3 | Repeats "In this paper" from L145; touches your I-06 edit, so your call. |
| C10 | L158 | "approximate and was chosen by calibration (SI Section 15.2)." | "approximate and was chosen by calibration." | 3 | The next sentence but one points to SI 15.2 again. |

Lower-priority repetition, not proposed: the "incidence guides allocation only" point appears at
L147, L185 (M-21, FCL-047) and in the Figure 1 caption; each answers a different reader question.

## On the recommended/optional split

- Drop T8 from the recommended set; swap T13 for the Discussion deletion described above.
- N8 is safe enough to be recommended.
- With the edits above (T8 out, T13 swapped, N12 pointer added, R3 lengthened), the recommended
  set measures 3,476. That leaves 24 words, which a reading-flow pass (W-03) can easily use up.
  Adding N8 and C1, C2, C4, C6, C7, C8 and C5 brings it to about 3,375 (3,382 measured before N8);
  N5 would add about 24 more.
