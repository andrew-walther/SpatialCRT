# CTJ V5 final proofread (read-only)

Date: 2026-10-27. Scope: `paper/ctj_manuscript/manuscript_ctj/CTJ_Manuscript.tex`,
`CTJ_Exhibits.tex`, the built `CTJ_Reading.pdf` (9 pp.) and `CTJ_Manuscript.pdf` (23 pp.),
`figures/Figure_Legends.md`, `submission_checklist.md`. Nothing was edited except this file.

Checks that passed:
- Every prose sentence in both sources appears in both PDFs, so the PDFs match the sources. `CTJ_Exhibits.tex` (23:37) is newer than the PDFs (21:11), so rebuild before committing, or the sync hook will warn.
- No `??` or undefined references. The four legends in `Figure_Legends.md` match the `\caption` text exactly.
- Exhibits are numbered in order of first citation (Table 1, Figure 1, Table 2, Figure 2, Figure 3, Figure 4). Panel cross-references (Figure 1C) are correct.
- Every abbreviation is defined at first body use: SRS, MSE, SUD, NC, SI, ML, SAR, MCSE. SEER is defined but used once.
- Design names in text, tables and figures follow the approved set. The figure PDFs use a true × in "2×2 Blocking", and "spatial blocking" is used for the NC adaptation.
- Citations follow punctuation throughout.
- The body is 3,499 words (texcount: 3,444 text + 55 headers; 96 math units excluded) and the abstract is 404 words, so the author's counts are confirmed.
- The only highlighted placeholders are the two known ones. See item 10 for unhighlighted draft notes.

## Proofreading and formatting findings

1. **Methods, Outcome model ¶3: "With queen contiguity Checkerboard is identified"**. Problem: there is no comma after the introductory phrase, while the parallel sentence that opens the same paragraph has one ("With rook contiguity, Checkerboard gives…"). Fix: "With queen contiguity, Checkerboard is identified".
2. **Serial comma (the manuscript otherwise leaves it out of simple lists, e.g. "clinics, schools or geographic areas", "population, clusters and counties")**. Problem: the following lists use it:
   - Introduction ¶3: "additional strategies, and an irregular observed geography". Fix: delete the comma before "and".
   - Methods, Simulation study ¶1: "($\rho_X=0.2,0.5$), and Poisson incidence". Fix: delete the comma before "and".
   - Abstract, Methods: "added an incidence-related baseline to the outcome, and ordered regions". Fix: delete the comma before "and".
   - Discussion, Limitations: "when treatment effects vary between clusters, or when outcomes are measured". Fix: delete the comma before "or".
   (Leave the commas in lists whose items contain "and", such as Methods ¶ "Primary queen settings…", Continued work "…deliver the education, and set…", and Supplementary material. The Sage Funding wording "research, authorship, and/or publication" is standard text, so keep it.)
3. **Title page (submission): "Main-text word count: 3499."** Problem: the text elsewhere uses a thousands separator (1,000; 23,523; 25,594,321). Fix: "3,499".
4. **Keywords: "…incidence-guided design; sudden unexpected death."** Problem: the keyword list ends with a period. Fix: delete the final period in `\CTJKeywords`.
5. **Table 1, 2×2 Blocking NC cell: "k-means clustering on location, incidence rank and population"**. Problem: Methods says the clustering uses "standardized location, incidence rank and log-population". Fix: change the cell to "location, incidence rank and log-population".
6. **Table 1, High Incidence Focus NC cell: "Treat the 29 clusters with the highest observed rate."** Problem: every other row ends its cells with a "(n treated)" count. Fix (optional for consistency): "Treat the clusters with the highest observed rate (29 treated)."
7. **Table 2 column headers: "Ratio to SRS (95% CI)"**. Problem: "CI" is never defined in the caption or the text, which say "confidence intervals" / "95% interval". Fix: add "(CI)" after "95% $t$-based confidence intervals" in the caption.
8. **Ranges written with "to" in place of an en dash**. Problem: ranges elsewhere use `--` (22--24%, 27--31, 21.1--29.9%). "Results, SUD application ¶6: "needed 35 to 37 colleges"" does not. Fix: "needed 35--37 colleges". "ranged from … to …" constructions and the negative bias range "$-0.12$ to $-0.40$" are correct as written.
9. **Value-list notation for ρ and γ**. Problem: one sentence mixes notations. Methods ¶ "Primary queen settings crossed $\rho=0,0.2,0.5$ with $\gamma=0.5,0.8$ … rook contiguity at $\rho\in\{0,0.5\}$ and $\gamma\in\{0.5,0.8\}$". Simulation study ¶1 uses the "=" list form. Fix: use one form within the sentence, e.g. "$\rho\in\{0,0.2,0.5\}$ with $\gamma\in\{0.5,0.8\}$".
10. **Unhighlighted draft notes that must be removed before submission**:
    - Submission title page: "Clinical Trials submission draft. Final submission approval remains pending."
    - Submission running header on every page: "Clinical Trials Submission draft (review)", from `\def\volumenumber{Submission draft\ }` / `\issuenumber{review}`.
    - Reading title block: "Reading preview only. Final author approval remains pending." This is acceptable in the preview but listed for completeness.
    Fix: delete the title-page line and blank `\volumenumber`/`\issuenumber` (or `\journalname`) for the submitted file.
11. **Submission layout: Tables 1 and 2 overflow the text block** (`.build/CTJ_Manuscript.log`: "Overfull \hbox (16.23997pt too wide) … lines 286" and "(16.42804pt too wide) … lines 288"). The reading layout's Table 1 is 8.7pt too wide (line 174). Fix: reduce the Table 1 column fractions (e.g. 0.37→0.35 and 0.28→0.27\linewidth). For Table 2, set `\tabcolsep` to 4pt or use `\footnotesize` in the submission branch.
12. **Figure 1 panel A title: "A. Observed 2018-2021 SUD rates (pooled)"**. Problem: the year range prints with a hyphen (soft hyphen `­` in `code/20_manuscript_figure_revision.R` line 248), while the caption and text use an en dash (2018–2021). Fix: use `–` in that label and re-export `nc_planning_map.pdf`.
13. **Regime labels across figures**. Problem: Figures 2 and 3 strips read "Spillover reaches both arms" / "Spillover reaches controls only", while Figure 4 reads "Both arms" / "Control only". Fix: use the same pair in all three figures (and, if the author wishes, match the text's "both-arms"/"control-only").
14. **Figure 2 legend: "Independent uniform"**. Problem: the text calls this configuration "independent Uniform$(0,1)$ incidence". Fix: change the legend to "Independent Uniform" (or "Independent").
15. **MSE decimal places**. Problem: grid MSEs have 4 decimals (0.1355, 0.2146) and application MSEs have 5 (0.30761, 0.07638, 0.07428, 0.07682), in the text and in the Figure 4 caption. Fix: if not deliberate, round the application values to 4 decimals (0.3076, 0.0764, 0.0743, 0.0768).

### References (as rendered by SageV)

16. **Ref 16 (`seer_population_data`): "U.s. county population data 1969–2023"**. Problem: BibTeX lowercased the acronym. Fix: in the bib, `title = {{U.S.} County Population Data 1969--2023}`.
17. **Refs 1 and 10 (`donner_design_2000`, `bivand_applied_2008`): "ISBN 978-0-470-71100-2 978-0-340-69153-3" and "ISBN 978-0-387-78170-9 978-0-387-78171-6"**. Problem: each prints two ISBNs. Fix: keep one ISBN per entry or delete the `isbn` fields; Vancouver does not need them.
18. **Article-title capitalization is mixed**. Title case: refs 3 ("Toward Causal Inference With Interference"), 11, 12, 13, 17. Sentence case: refs 4, 5, 6, 7, 8, 14. Fix: convert 3, 11, 12, 13 and 17 to sentence case, protecting proper nouns with braces (e.g. `{North Carolina}`, `{United States}`).
19. **Ref 2 (`hayes_cluster_2017`): "Cluster Randomised Trials, Second Edition"**. Problem: the edition is part of the title. Fix: `title = {Cluster Randomised Trials}` with `edition = {2nd}`.
20. **Ref 9 (`walther_spatialcrt_2026`): "Spatial cluster randomized trials - sampling design with spillover effects & spatial dependence"**. Problem: a spaced hyphen is used as a dash. Fix: match the published title exactly; if it has a dash, use `--`, otherwise a colon.
21. **DOI and URL redundancy is inconsistent**. Refs 2, 3, 4, 5, 9, 11, 12, 13 and 17 print both a DOI and a publisher URL, while refs 6, 8 and 14 print the DOI only. Fix: remove the `url` fields from entries that have a DOI, keeping URLs only where there is no DOI (ref 16).

## Journal requirements

**Source note:** WebFetch of https://journals.sagepub.com/author-instructions/ctj returned **HTTP 403 Forbidden**. A curl retry was blocked, and web search did not turn up the CTJ page. These findings therefore rely on `submission_checklist.md` and the guideline text recorded when the live page was checked on 2026-10-02:
- `docs/findings/2026-10-02-manuscript-completion-findings.md` (lines 146–154): "original articles permit 3,500 body words and six exhibits; structured abstracts permit 425 words (Background/Aims, Methods, Results, Conclusions)".
- `docs/findings/2026-10-02-ctj-layout-review.md` (lines 17–20): "request double-spaced, minimally formatted review text; tables on separate sheets at the end; separate figure legends; no more than six main exhibits".

Re-check the live page before submission.

J1. **Body length, compliant with no margin.** It is 3,499 of 3,500 (texcount confirmed). Any net addition of 2 or more words exceeds the limit.

J2. **Abstract, compliant.** 404 of 425 words, under the four required headings "Background/Aims", "Methods", "Results" and "Conclusions".

J3. **Exhibit count, compliant at six (2 tables + 4 figures), with one risk.** Table 1 embeds eight graphics (`figures/compact/design_thumbs/thumb_*.pdf`). Sage production generally expects tables as editable text, so a table with images may be treated as a figure, which would make seven exhibits. This is unverified because the live guidelines could not be fetched. Confirm with the editorial office or the live instructions.

J4. **Ethics statement missing from the title page.** `submission_checklist.md` records: "(Sage also asks for the statement on the title page)". The title page has no ethics statement; it appears only in the "Ethical considerations" declaration and the Methods. Fix: once the determination wording is final, add the ethics statement (IRB protocol 14-2036 and the determination) to the title page.

J5. **Required declarations are present, but two still hold the known placeholders.**
   - Ethical considerations: "[Placeholder for the IRB determination and consent statement, pending confirmation from Dr. Simpson.]"
   - Author contributions: "All authors critically revised the manuscript and approved the final version. [Placeholder: pending coauthor approval.]"
   The other declarations are present: conflicting interests, funding, data availability, AI disclosure (in Acknowledgements, which the checklist says Sage allows in "methods or acknowledgements") and supplementary material.

J6. **Conflicting-interests wording is not the Sage standard sentence, while Funding is.** Current text: "The authors declare that there is no conflict of interest." Sage's standard declaration (from Sage's general template, not the fetched page) is "The author(s) declared no potential conflicts of interest with respect to the research, authorship, and/or publication of this article." Fix: use the Sage sentence so it matches the Funding statement.

J7. **Funding not confirmed for every author.** The tex comment above Funding reads "confirm with Habib and Simpson", but the checklist ticks "[x] Confirmed current author list, corresponding author, no funding/conflicts." Fix: confirm with both coauthors, then remove the comment or untick the item.

J8. **Coauthor ORCIDs outstanding.** The title page lists Walther's ORCID only. The checklist item is still open.

J9. **Draft labels on the submission file.** See item 10. The title-page footer and running header must not say "submission draft" in the submitted PDF.

J10. **Checklist internal inconsistencies (`submission_checklist.md`)**:
   - "[x] Generative-assistance acknowledgement removed from the draft (author, 2026-10-12)." This contradicts the later "[x] Generative-AI disclosure … Draft added as Acknowledgements (2026-10-14)" and the current manuscript. Fix: mark the 2026-10-12 item as superseded.
   - "Checked 2026-10-02 against journal instructions": the live page has not been re-verified since then, and today's fetch failed. Fix: re-check and update the date.
   - "Declarations not recounted since 2026-10-14 (266 words then)": stale now that Ethical considerations has been added. Fix: recount or drop the figure.
   - Files list "the four referenced vector PDFs in `figures/compact/`" omits the eight Table 1 thumbnails in `figures/compact/design_thumbs/`. Fix: add them to the files list, because they must be uploaded or the table supplied as rendered.

J11. **References follow Sage Vancouver (SageV numeric superscripts after punctuation), as the recorded check expects.** There are 17 references. The recorded notes give no reference limit, so this is unverified. Fix the reference-formatting errors in items 16–21.

J12. **Keywords: six are present, as the checklist records.** Remove the trailing period (item 4).

## Live guideline check (2026-10-28)

The page loaded in the app's browser on 2026-10-28, so the items above that relied on recorded notes were re-checked against the live text.

- **Confirmed:** 3,500 body words excluding abstract, references, tables and figures; six tables or figures; abstract ≤ 425 words under Background/Aims, Methods, Results and Conclusions; 3–10 keywords; Sage Vancouver references; double spacing, left-ranged unjustified text without hyphenation; tables at the end on separate sheets; legends on a separate sheet and not on the figures; running head ≤ 40 characters (ours has 32); title page with word count, first names, affiliations with city/state/country, and the corresponding author's address, email and telephone.
- **Corrected (J6 was wrong):** section 2.5 gives the conflict sentence "The Author(s) declare(s) that there is no conflict of interest", which the original text matched. Restored as "The authors declare that there is no conflict of interest."
- **Corrected:** section 2.4 gives the no-funding sentence "This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors." This replaces the generic Sage sentence.
- **J4 withdrawn:** no title-page ethics statement is required. Section 2.6 asks Methods to state that the IRB approved or waived approval, give its full name, institution and approval number, and say whether consent was obtained (written or verbal). Name, institution and number are in Methods; the determination and consent wording wait on Simpson.
- **J8 softened:** the system collects the corresponding author's ORCID only.
- **New:** section 4.1 strongly discourages acronyms except universally recognized ones. Author decision.
- **Unchanged:** J3 (Table 1 with images) is not addressed by the guidelines; still a question for the editorial office.
