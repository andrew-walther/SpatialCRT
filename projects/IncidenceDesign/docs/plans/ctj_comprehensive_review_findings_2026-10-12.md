# CTJ comprehensive review — consolidated findings (2026-10-12)

Four read-only tracks run against the state after commit `2b142f9` (brief:
`claude_comprehensive_review_prompt_2026-10-12.md`). No edits were made by reviewers.

## Headline

- **Numbers: clean.** `verify_manuscripts.py` and `test_manuscript_exhibits.R` pass; every
  abstract/body/table/caption number and ~60 SI cells were recomputed from source CSVs/RDS and
  match. Scratch rebuild: 0 undefined refs/citations. Chapter ↔ SI numeric diff: no conflicts;
  all 22 shared figures byte-identical.
- **Code: clean.** No issue changes a manuscript number. All cheap suites pass
  (application summary/tail/verify production+tail, allocation-risk summary, simulation
  revision with `SKIP_EQUIVALENCE=1`); code/20 and code/21 re-runs reproduce tracked CSVs byte for byte.
- **Journal limits met:** body 3,165 words (limit 3,500), abstract 279 (≤425, required
  headings), 6 exhibits (≤6), running head 32 characters. Live page is behind Cloudflare; rules
  taken from Sage's 2026-04-30 PDF capture of the CTJ guidelines + live Sage AI policy.

## Tier A — recommended fixes before sharing with the advisor

| ID | Sev | Where | Issue | Fix |
|---|---|---|---|---|
| A1 | ERROR | bib (Cai, ref 7) | Internal note printed: "key kept for existing citations" | Clean `note` to a proper ICLR 2024 citation |
| A2 | WARN | CTJ abstract, :164, :202, :225 | No caveat that incidence ordering adds nothing over plain saturation (SI/chapter say so; plain ≤ guided on every headline metric) | One sentence in Results/Discussion; decide abstract wording (Q1) |
| A3 | WARN | CTJ :176 | Habib citation attached to the *corrected* counts, but the cited report used the earlier numerator; coauthor cited as "personal communication" | Cite report for the original dataset; state the correction plainly |
| A4 | WARN | CTJ :164 | "Quartile and half stratification" still mentions dropped Balanced Halves | "Quartile stratification distributes…" |
| A5 | WARN | Table 1 / Figs 3–4 / :216,:220 | NC design names differ (2×2 Blocking vs Spatial blocking; Checkerboard vs Graph Checkerboard; "2 x 2") | Add NC names in Table 1 NC column; one spelling of 2×2 |
| A6 | WARN | body | SRS/SUD/MSE defined only in abstract; SAR, MCSE, MC, SEER, q90 never defined | Define at first body use; drop MC/MCSE where possible |
| A7 | WARN | :171,:181,:187,:189,:209,:233 | Pipeline jargon ("named cached W", "computational sources", "supports", "masks", "copied prefixes", "opposite outcome halves", manifest/checkpoint sentence) | Plain-language rewrites (Track 2 wording) |
| A8 | WARN | :211 | Garbled "including setting-specific approximate MC intervals" | "…and approximate Monte Carlo intervals for each setting-specific difference excluded zero" (verified true: 48/48) |
| A9 | WARN | :202 | Pooled six-design ranking claim reads against regime-specific framing | Say "pooled across regimes", note SRS not in that ranking, both-arms leader differs |
| A10 | WARN | :157 / Fig 2 | ρ_X undefined; X not named as incidence vector | "(ρ_X = 0.2, 0.5)" |
| A11 | WARN | :169 | Unclear subject ("Under queen control-only…") | "For queen Checkerboard, under…" |
| A12 | WARN | refs 6, 9 | Baird DOI loses underscore; BMC ref lacks volume/article no., " - " in title | Escape DOI; add volume/article number; fix dash |
| A13 | WARN | Table 1 PDF | First column fully justified ("High      Incidence Focus") | Ragged-right column |
| A14 | WARN | contributions :244 | No "all authors critically revised and approved" | Append sentence |
| A15 | WARN | Cover_Letter_Draft.md | Says eight strategies, a disclosure, and aggregate outputs accompany — all now false | Align with current manuscript |
| A16 | WARN | SI :1084–1087, :2035 | Still says "public repository" and "ethics statement pending" — contradicts Oct-12 CTJ data statement | Align SI with CTJ |
| A17 | NOTE | :193 | "approximately 94%" grid coverage; table spans 0.932–0.942 | "approximately 93–94%" |
| A18 | NOTE | :193 | Incidence config changes MSE "only slightly" (HIF varies 19%) | "modestly (under 10% except High Incidence Focus)" |
| A19 | NOTE | :222 | 3.68–4.44 / 1.42–2.52 are setting-level ranges | add "across settings" |
| A20 | NOTE | Table 1 | Buffer "about 21 treated" is queen (rook 38) | add "(queen)" |
| A21 | NOTE | PDF | Keywords run into abstract; header "Submission draft(review)" | `\par`; add space |
| A22 | NOTE | Fig 3 caption | Colors not explained | "SRS orange; designs blue" |
| A23 | NOTE | SI | "the the" (:227); "detailed grid analysis of the detailed analysis" (:1465,:1518); "primary analysis in the detailed analysis" (:33); dated Oct 2; downside caption should say ratio of averages (:1343) | Fix wording/date |
| A24 | WARN | chapter :742 | "spread small for every randomized design" — false for Isolation Buffer | "for SRS, balancing and saturation designs" |
| A25 | WARN | chapter :351, :738 | "first or second under both regimes" true only among six designs (5th of 9 both-arms) | add "among the six designs" |
| A26 | NOTE | chapter :626/SI :1256; :655/SI :1311; :732 | "all nine designs" vs eight-row heatmap; "Random sampling" header; Discussion "six designs" | Reword |

## Tier B — author/coauthor decisions (not fixable by editing alone)

| ID | Sev | Issue |
|---|---|---|
| B1 | ERROR (submission) | No ethics/IRB statement; CTJ requires one in Methods for human data. Needs Habib/Simpson IRB determination. |
| B2 | ERROR (submission) | Sage AI policy requires disclosure of AI-generated/checked code and visualization (methods or acknowledgements). Author removed it on Oct 12. |
| B3 | WARN | Title page requires corresponding author's telephone. |
| B4 | WARN | Public repo tracks 58-cluster × year death counts (`observed_cluster_inputs.csv`, `setup.rds`; min cell 5). Confirm the data-use agreement's small-cell rule permits this. |

## Tier C — deferred / hygiene (no reader impact)

- Script headers credit AI tools (code/20 "Codex", code/21 "drafted with Claude Code",
  CTJ_Exhibits/CTJ_Reading "Author: Codex"); standing rule is "Author: Andrew Walther" only.
- code/21 `aggregate()` would silently drop NA (none present) → add finiteness asserts;
  code/20 docs say "four figures"; Balanced Halves `sample()` latent single-element hazard;
  killed-worker NULL not caught early in application runners; quantile-type inconsistency.
- Figure ticks (Fig 2 left limit, Fig 3 past last tick, Fig 4 color bar to 2); reference ISBN
  clutter; unreferenced figure files; companion `regions.png` still pre-"fixed" label;
  `sources.txt`, `CTJ_Exhibits.tex` header date, AGENTS page counts stale.
- Optional: ADEMP reporting clause; cite SEER and Mirzaei base rate; name data custodian.

## Author decisions and resolution (2026-10-12)

- Applied: A1, A4–A13, A15–A26 and the chapter items; Lin telephone (919) 966-9431 from
  his UNC Gillings profile (B3 resolved); Codex/Claude mentions removed from code headers,
  TeX comments, cover letter, declarations draft and companion (Tier C headers).
- A2 reframed per author: no "no benefit" caveat in the CTJ. SI/chapter/report now say plain
  and incidence-guided saturation performed similarly, so incidence guidance targets
  high-burden regions without loss of accuracy; the advantage over SRS comes from the
  shared saturation structure.
- A3 not changed: Habib produced the corrected counts after amending the adjudication method.
- A14 (all-authors-approved sentence) deferred to coauthor approval.
- B1 (IRB), B2 (AI disclosure) moved to TODO; B4 presumed acceptable pending IRB permission
  to publish aggregates.
- BMC reference (A12): no volume/article number yet (online-first); title's " - " matches
  the publisher record, so it is unchanged.
- Deferred hygiene: code/21 NA guards, figure axis ticks, reference ISBN clutter, unused
  figure files, Balanced Halves `sample()` hazard, killed-worker handling.
