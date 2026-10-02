# Completed long-form chapter and derived CTJ drafts — 2026-10-02

The observed-incidence application, full thesis chapter/appendix and full CTJ
manuscript/supplement drafts are complete. This means analyses and reviewable
manuscripts, not journal submission or author approval of the final text.
Author-only ethics/data-use documentation remains an explicit placeholder.
The author superseded older calendar deadlines: complete promptly, with CTJ
remaining a required deliverable.

## Scientific interpretation and independent review

- Incidence guides primary NC education-outcome allocation; β = 0 in generation
  and no X in the primary fit. The grid study retains β = 1. The β = 1/X-adjusted
  application sensitivity changes baseline and adjustment jointly and assumes no
  empirically established incidence–education relationship. τ is the constant
  structural intervention coefficient before spatial propagation.
- SRS and all eight strategies appear in the main full-grid benchmark and NC
  application. Detailed grid rankings remain a labeled six-design subset. Plain
  saturation and Balanced Halves were not silently omitted from recommendations.
- Average MSE and the distribution of conditional MSE both matter. Saturation
  has favorable mean and descriptive allocation-tail results under control-only
  spillover. SRS remains competitive under both-arms spillover. Small differences
  establish neither equivalence nor uniquely optimal assignment.
- The grid allocation-risk pilot is supporting selected-surface/corner evidence;
  it never replaces full-study numbers. The NC risk extension adds outcomes for
  the same allocations. Tail membership remains uncertain; singleton allocation
  support does not imply negligible estimation error.
- The NC Design 8 mean-rank summary is an author-approved adaptation. Grid code
  averages X directly; only rank-scaled Poisson X gives the mean-rank analogue.
  Continuous grid incidence uses its values, not their ranks. Alternatives were
  assessed using frozen regions without choosing the primary by performance.
- The fresh independent agent reviewed study alignment, application verification
  and the completed documents. Three mandatory CTJ scope corrections were applied:
  zero aliases/precision success refer to the application; failure to recover τ
  refers to rook Checkerboard, not identified queen Checkerboard; mean rank does
  not exactly reproduce the continuous grid mean-X rule. Final source review
  reported no scientifically mandatory issue. This was an independent agent
  check, not a journal peer review or a Claude review performed by the author.
- 936 denotes computational sources, not necessarily distinct statistical
  distributions: buffer's regimes can share a DGP while retaining separate
  approved RNG streams. Reused yearly rows are not independent confirmation.

Project 1 remains final/read-only. The accepted work's practical argument about
reasonably good estimation and limiting poor-allocation downside is represented
within its studied conditions. The historical numerical discrepancy remains
unresolved. No direct numerical cross-project comparison or renewed audit was
performed. The chapter uses the approved conceptual Chapter 2 link; CTJ/SI
instead cite the verified published paper and stand alone.

## Long-form source and derivation

Canonical source: `paper/dissertation_chapter/Dissertation_Chapter.qmd`.
It retains detailed September-revision grid methods/results, adds a full
nine-design/SRS benchmark and properly separated grid pilot, replaces the old
application with completed education methods/results, and adds Appendix A7
for reproducibility, geography, annual tables and matched sensitivities.
Technical identification details moved to Appendix A6, with a short main-text
explanation and explicit reference.

`CTJ_Manuscript.tex` trims and reorganizes that source into Introduction,
Methods, Results and Discussion plus declarations. It is a full article, not
an outline. `Supplementary_Information.tex` preserves the source's full methods,
detailed grid results, completed application and seven appendix components as
S1–S10. All 22 source table bodies are retained. CTJ omits dissertation-only
references and uses a self-contained literature link.

Earlier live chapter/CTJ/SI sources and the old declaration abstract were saved
in `paper/archive_manuscript/pre_completion_20261002/`. Previously generated
synthetic/application results remain preserved and unused in the revised drafts.
The current declaration file now carries the current abstract rather than
copying its superseded April results into bios-dissertation.

## Resolution of former NEEDS-AUTHOR-CONFIRMATION comments

The old application background/proposal was rewritten for the actual study.
Removing an unsupported claim is not author confirmation of that claim.

| Former item | Resolution |
|---|---|
| Uncited broad SUD burden claim | Removed; focus is the specific education-pilot planning question. |
| SUDDEN validation / expert panel claim | Removed; no upstream validation claim inferred from a bibliography entry. |
| Mirzaei-to-SUDDEN methodological connection | Removed. Grid base-rate example remains explicitly general, not calibration to observed NC rates. |
| Unverified household-income gradient detail | Removed from the application rationale. |
| Satellite campus example / comprehensive service commitments | Removed. The implemented 100-county/58-cluster analytical partition is described as a study convention. |
| Intercollege programmatic communication | Replaced by a geographic exposure proxy; no observed communication network is claimed. |
| Lenoir caption assertion | Removed from the map caption. |
| Unverified upstream case exclusions | No undocumented exclusion reconstructed; input reconciliation is documented, including corrected heart-failure numerator provenance. |
| Planned endpoint/protocol specifics | Education outcome clarified by author; participant eligibility, implemented endpoint and trial protocol remain future decisions, not fabricated details. |
| Geographic allocation adaptations | Author approved retaining mechanisms with frozen four-region/yearly-block adaptations, differing budgets and 29-treatment Balanced Halves correction; implemented and verified. |
| Application results pending | Replaced with complete verified production and separate tail refinement. |
| IRB / restricted-source permissions | Explicit author-approved placeholder retained; no approval, waiver or consent determination asserted. |

## Number and exhibit provenance

| Content | Authoritative source |
|---|---|
| Full nine-design grid, τ sweep | `results/sim_data/sim_results_MLE_tau_sweep_combined_20260927_000502.rds`; 14,400 scenarios per estimator. Original eight-design rows reproduce September 24. |
| Main regime-specific grid benchmark | `results/srs_benchmark/manuscript_regime_means.csv`, generated by `code/19_manuscript_exhibits.R`; τ=1, queen, 80 scenarios per design/regime. |
| Existing detailed grid means/ranks/tests/configuration splits | `results/six_design_manuscript/{six_design_comparison_report.rds,six_design_summary.txt,dissertation_results_extract.txt}`, queen/rook assets and non-oracle counterpart. |
| Grid technical decomposition/rare aliases | Chapter `notes/a5_technical_notes.md` and its recorded rerun; A6 explicitly distinguishes rerun diagnostics from stored full-run metrics. |
| Selected-surface grid risk | `results/allocation_risk/summary/` and `docs/plans/allocation_risk_findings_2026-09-27.md`; pilot and refinement kept separate. |
| Production counts/precision/diagnostics | `application/results/real_sud_rev_20261002/production/{performance.csv,results.rds,manifest.rds,verification.txt}`. |
| Annual primary MSE, bias, coverage, budgets | `exhibits/yearly_primary_design_means.csv`; six equally weighted ρ/γ pairs, each year/regime separate. |
| Within-setting SRS differences | `exhibits/primary_setting_srs_comparison.csv`; approximate setting-specific MC intervals, no unsupported pooled intervals. |
| Hypothetical baseline, matched rook and regional alternatives | `exhibits/yearly_sensitivity_means.csv` and `sensitivity_setting_matched_comparison.csv`; denominator matches year/design/parameters exactly. |
| Same-allocation refined downside | `tail_confirmation/` and `exhibits/refined_tail_comparison.csv`; copied prefixes counted once. |
| Input maps, regions and assignment examples | Frozen `setup.rds`, named cached weights, authorized observed-cluster inputs and `application/report/real_sud_assets/`. Examples illustrate rules, not performance-selected masks. |

Paths in the final six rows are below `application/results/real_sud_rev_20261002/`.
`paper/tools/verify_manuscripts.py` compares 36 grid cells and 20 primary annual
MSE cells in all three source documents; 360 full annual performance cells,
72 budget cells and 32 matched-sensitivity cells to current CSVs. It also checks
all 22 supplement tables against the chapter, source references, bib keys and
figure paths. A deliberately erroneous-number fixture is rejected. This is a
numerical/source check; independent review assesses interpretation and scope.

## Figure/table selection

| Placement | Material retained |
|---|---|
| Chapter body | Detailed grid design rules, main grid accuracy/coverage/ranks/configuration/effect-size findings, full nine-design/SRS benchmark, observed geography/rate maps, education model, budgets, annual primary results and refined risk. |
| Chapter appendix | Detailed eight-design subset justification, random-number/MC mechanics, statistical diagnostics, rook figures, non-oracle results, identification/error decomposition, frozen NC regions/examples, setting-specific SRS ratios, all annual nine-design metrics and sensitivity ratios. |
| CTJ Table 1 | Full nine-design queen grid benchmark separated by spillover regime. |
| CTJ Table 2 | Four separate yearly primary NC MSE rows: SRS, plain and guided saturation as relevant to each regime. |
| CTJ Figure 1 | Current detailed six-design grid MSE figure; main table supplies all nine designs. |
| CTJ Figure 2 | All nine NC primary yearly MSE comparisons. |
| CTJ Figure 3 | Refined conditional mean/q90/worst-decile error, both regimes, same allocations. |
| CTJ Figure 4 | Annual treated-population shares/ranges; cluster counts remain in text/SI. |
| CTJ supplement | Full source detail, observed maps, all annual nine-design bias/coverage/budgets, matched sensitivities and explicit diagnostics; no substantive finding discarded to meet the main cap. |

Journal versions omit embedded titles; captions explain their scope. The revised
risk plot includes its stated worst-decile series as well as mean and q90.
Separate legends are saved in `paper/ctj_manuscript/Figure_Legends.md`.

## Live publication/reference checks

Checked October 2 against [Clinical Trials instructions](https://journals.sagepub.com/author-instructions/ctj):
original articles permit 3,500 body words and six exhibits; structured abstracts
permit 425 words (Background/Aims, Methods, Results, Conclusions). The older
plan's 250-word unstructured abstract expectation is superseded. The draft has
2,217 main-text words and a 283-word abstract (texcount fragments; body includes
headings, excludes abstract/declarations/references/exhibits). It has six keywords,
a separate author/title page, a 33-character running head, double-spaced review
text, tables/figures after references, Vancouver references and declarations.

Author confirmed Walther/Simpson/Habib/Lin, Lin corresponding, no current funding
or conflicts. Telephone/ORCID, contribution wording, final coauthor consent and
ethics/data-use documentation remain author-side submission tasks. No actual
trial was conducted; no trial registration or CONSORT result is invented.
See `paper/ctj_manuscript/submission_checklist.md` and the cover-letter draft.

[Publisher page for the preceding work](https://link.springer.com/article/10.1186/s12874-026-02968-0)
verifies Walther, Van Deinse and Lin, the title, BMC Medical Research Methodology,
DOI 10.1186/s12874-026-02968-0 and online publication August 26, 2026.
The bib entry `walther_spatialcrt_2026` uses these details, without guessing a
volume/page/article number. No Project 1 source file was changed.

The [Sage AI policy](https://www.sagepub.com/journals/publication-ethics-policies/artificial-intelligence-policy)
requires disclosure when generative assistance affects code/text/visuals rather
than language polishing alone. An acknowledgement disclosure draft states the
Codex assistance and author responsibility, without claiming that all authors
have already completed review. There are no AI co-author lines in Git commits.

## Verification and rendering

- Application data, allocation/model behavior, tail diagnostics and summary tests
  pass. Production and separate tail verification independently reconciled every
  cached source; all 1,248/144 reporting settings are complete and satisfy approved
  mean/coverage precision gates. Full earlier ML/grid validation remains documented;
  actual NC engine comparisons were checked against lagsarlm by the independent reviewer.
- The companion's actual JavaScript passes a DOM fixture covering nine-design
  selection, planned cells, main/refined filters, value scaling and withholding
  false tail intervals. Its maps/plots were inspected as images. In-app file-URL
  browsing is blocked; no browser visual review is claimed.
- Fresh chapter, CTJ and supplement builds succeed. A visual pass inspected all
  chapter/SI pages and CTJ pages; a final pass checks corrected title/abstract,
  legends and updated plots. No unresolved citations/references or clipped
  exhibits were found. Final CTJ/SI logs contain no overfull boxes or unresolved references;
  expected float-only-page/spacing notices are not model warnings.
- Local Quarto uses a temporary Lua font cache and disables automatic package
  installation. CTJ uses the existing Sage Review class in a single-column,
  double-spaced A4 layout with a separate breakable abstract. A custom author
  title sheet avoids the class's unbreakable production-layout abstract.
- Chapter checkpoint `58c1c82` invoked the approved hook successfully: bios-dissertation
  commit `fc63326`, 72-page generated Chapter 3/Appendix B, seven figure assets
  copied and 31 cited keys checked; no push. Unrelated bios-dissertation changes
  were not modified or staged by this work.
  No direct edits to bios-dissertation were made; all generated changes used the hook.

## Teach-back for changes in this completion pass

**`real_sud_summary.R` export block:** Purpose: make the completed results usable
in both the companion and journal. Logic: retain the original titled plot, save
an extra PDF with title/subtitle/caption removed, and include worst-decile means
in the risk plot. Inputs are verified yearly/tail CSVs; outputs are four regular
and four journal PDFs plus existing PNGs. ggplot2 changes presentation only.
There is no recomputation of τ estimates, fit caches or model manifests. Risk
quantiles remain descriptive noisy conditional estimates, without invented CIs.

**`render_real_sud_companion.R` status/link block:** Purpose: let the author review
completed analyses and drafts from one offline page. It reads current verified
aggregate outputs, reports computational-source counts and links the completed
PDFs/review prompt. Input/outputs and DOM behavior are unchanged. Relative local
links require keeping the project layout; the page needs no external scripts.

**`code/19_manuscript_exhibits.R`:** Purpose: trace the current full-grid benchmark
and export its six-design journal figure. It reads the named verified nine-design
oracle RDS, asserts its scenario/fit counts, selects τ=1/queen and averages 80
scenarios per design/regime. Then it plots the detailed six-design subset in
pooled-MSE order, preserving the chapter's overlaid-configuration interpretation.
Inputs are the current 14,400-row RDS; outputs are an 18-row aggregate CSV, an
MD5 provenance note and one vector PDF. Existing visualization helpers and
`ggplot2` avoid a new dependency. It runs serially and performs no simulation.
If the study/version changes, its explicit source/count assertions require a
reviewed update rather than silently selecting whichever result file is newest.

**`paper/tools/verify_manuscripts.py`:** Purpose: catch transcription drift and
broken derivation. `csv_rows` reads aggregate dictionaries; `table` selects an
explicit labeled tabular body; `numerical_rows` extracts printed cells while
excluding headers; `require_equal` fails with actual/expected values; `verify`
assembles independently rounded CSV expectations and compares all target cells,
then checks source table derivation, assets, refs, citations and exhibit count.
Its local `normalized` helper ignores comments/whitespace and the approved
Chapter-2-to-standalone phrase change, retaining table cells and commands.
Inputs are three source documents, aggregate CSVs and the shared bib; output is
PASS text or a nonzero assertion. Python stdlib avoids new dependencies. Exact
format comparisons are intentional: changed rounding/order requires review.
The negative fixture changes a number and proves that equality verification
rejects it. This does not replace prose review or a PDF layout inspection.


Final standalone PDFs: chapter 65 pp, CTJ 18 pp, supplement 43 pp. The generated
prelim copy has 72 pp under its own chapter/appendix layout. Final verification
text is `manuscript_verification_2026-10-02.txt`; a SHA-256 source/output manifest
is `paper/manuscript_provenance_20261002.json`. LaTeX intermediates are scoped
ignored; frozen archived TeX sources are explicitly retained.

Final geographic wording uses “projected cluster coordinates,” not geometric
centroids: the frozen setup uses polygon representative interior points. This
clarifies the compactness/distance description without changing any allocation,
result, seed or manifest. The final chapter revision is synced by the same hook.

Final manuscript checkpoint: `a47b0ae`; its approved hook resynced the coordinate
wording to bios-dissertation `6b07b95` (72 pp, no figure changes, 31 citekeys).
Generated chapter/appendix numbering was checked in extracted text and sampled
rendered pages. Project Git status was clean after the manuscript checkpoint.
The final metadata checkpoint records this result and refreshes the source hashes;
no data, allocation, estimator or manuscript finding changed.
