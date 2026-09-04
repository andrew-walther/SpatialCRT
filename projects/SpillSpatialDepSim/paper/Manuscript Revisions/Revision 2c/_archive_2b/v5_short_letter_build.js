const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType, VerticalAlign,
} = require("/tmp/v5build/node_modules/docx");

// ---- layout constants (US Letter, 1" margins -> content width 9360 DXA) ----
const CONTENT = 9360;
const COL = [2680, 4360, 2320]; // comment | response | location  (sums to 9360)
const HEADER_FILL = "D9D9D9"; // light gray (B&W safe)
const gborder = { style: BorderStyle.SINGLE, size: 4, color: "808080" };
const borders = { top: gborder, bottom: gborder, left: gborder, right: gborder,
  insideHorizontal: gborder, insideVertical: gborder };
const cellMargins = { top: 80, bottom: 80, left: 120, right: 120 };

// ---- text helpers ----------------------------------------------------------
function para(text, opts = {}) {
  const runs = Array.isArray(text) ? text : [new TextRun(text)];
  return new Paragraph({
    children: runs,
    spacing: { after: opts.after ?? 120, before: opts.before ?? 0, line: 276 },
    alignment: opts.align,
  });
}
// build paragraphs from text where **bold** and *italic* are inline-marked
function rich(s) {
  const parts = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0, m;
  while ((m = re.exec(s)) !== null) {
    if (m.index > last) parts.push(new TextRun(s.slice(last, m.index)));
    const tok = m[0];
    if (tok.startsWith("**")) parts.push(new TextRun({ text: tok.slice(2, -2), bold: true }));
    else parts.push(new TextRun({ text: tok.slice(1, -1), italics: true }));
    last = re.lastIndex;
  }
  if (last < s.length) parts.push(new TextRun(s.slice(last)));
  return parts;
}
// multi-paragraph cell content: array of strings -> array of Paragraphs
function cellParas(arr, opts = {}) {
  return arr.map((s, i) =>
    new Paragraph({ children: rich(s), spacing: { after: i === arr.length - 1 ? 0 : 100, line: 264 } }));
}
function headerCell(txt, w) {
  return new TableCell({
    borders, width: { size: w, type: WidthType.DXA }, margins: cellMargins,
    shading: { fill: HEADER_FILL, type: ShadingType.CLEAR },
    verticalAlign: VerticalAlign.CENTER,
    children: [new Paragraph({ children: [new TextRun({ text: txt, bold: true })] })],
  });
}
function bodyCell(parasArr, w) {
  return new TableCell({
    borders, width: { size: w, type: WidthType.DXA }, margins: cellMargins,
    children: parasArr,
  });
}
// three-column item table: comment / response / location  (each = array of strings)
function itemTable(comment, response, location) {
  return new Table({
    width: { size: CONTENT, type: WidthType.DXA }, columnWidths: COL,
    rows: [
      new TableRow({ tableHeader: true, children: [
        headerCell("Reviewer Comment", COL[0]),
        headerCell("Author Response", COL[1]),
        headerCell("Location & Evidence", COL[2]),
      ]}),
      new TableRow({ children: [
        bodyCell(cellParas(comment), COL[0]),
        bodyCell(cellParas(response), COL[1]),
        bodyCell(cellParas(location), COL[2]),
      ]}),
    ],
  });
}
function itemHeading(t) {
  return new Paragraph({
    spacing: { before: 260, after: 100 },
    children: [new TextRun({ text: t, bold: true, size: 24 })],
  });
}
function sectionHeading(t) {
  return new Paragraph({
    spacing: { before: 320, after: 140 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "000000", space: 2 } },
    children: [new TextRun({ text: t, bold: true, size: 26 })],
  });
}

// ============================================================================
// CONTENT
// ============================================================================
const children = [];

// ---- Title block (centered) ----
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 60 },
  children: [new TextRun({ text: "Response to Reviewers", bold: true, size: 32 })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
  children: [new TextRun({ text: "Spatial Cluster Randomized Trials — Sampling Design with Spillover Effects & Spatial Dependence", italics: true, size: 22 })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
  children: [new TextRun({ text: "BMC Medical Research Methodology — Second Revision (Condensed Version)", size: 22 })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
  children: [new TextRun({ text: "Andrew Walther, Tonya Van Deinse, Feng-Chang Lin", size: 22 })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 },
  children: [new TextRun({ text: "June 2026", size: 22 })] }));

// ---- Salutation + opening ----
children.push(para([new TextRun("Dear Editors (Michael Grayling & Vrushali Mule) and Reviewer,")], { after: 160 }));

children.push(para(rich(
  "We thank the Editor and the Reviewer for their continued, careful engagement with our manuscript. The second-round comments identified genuine weaknesses in how the model’s estimands, spillover specification, notation, estimated model, and inferential scope were presented, and addressing them has materially strengthened the work. We have carefully revised the manuscript in response to each of the six major comments and to the Editor’s overarching requests, **including the request to reduce the manuscript’s length**. This is the condensed version of the revision: it is identical in scientific content to the full revision and additionally implements the structural length reductions described in the **Manuscript Length** section below. All changes are implemented in the revised file **Revisions_V2b_short.tex**; the original submission has been left unmodified for reference."
), { after: 140 }));

children.push(para(rich(
  "For every comment below we restate each of the Reviewer’s comments, describe the change, and give its location in the revised manuscript — section number and page in **Revisions_V2b_short.pdf**, with equation references where relevant — to facilitate cross-referencing with the annotated manuscript (**Revisions_Annotated_V2b_short.pdf**), in which all added text is underlined and all deleted text is struck through. A consolidated “Summary of All Changes” is provided at the end of this letter. One notational note: the Reviewer writes the spillover parameter as φ; in our manuscript this parameter is denoted ψ, and we use ψ throughout our responses."
), { after: 60 }));

// ============================================================================
// EDITOR COMMENTS
// ============================================================================
children.push(sectionHeading("EDITOR COMMENTS"));
children.push(para(rich(
  "The Editor (Dr. Michael Grayling) concurred with the Reviewer’s report, noted that the previous revision had not yet addressed the deeper issues, and asked us to think carefully about each problem and revise thoroughly. The accompanying editorial instruction set three specific requirements: that results be **accurately reported**, that **overstated conclusions be rewritten**, and that the **limitations of the work be fully explained**. The Editor additionally asked that the manuscript’s length be reduced; this condensed version implements that reduction, summarized in the final row of the table below and detailed in the **Manuscript Length** section. The table below maps each requirement to the corresponding revisions; Dr. Grayling’s request to address each previously-raised problem is answered point-by-point in the Reviewer Comments section that follows."
), { after: 140 }));

children.push(new Table({
  width: { size: CONTENT, type: WidthType.DXA }, columnWidths: [2680, 4360, 2320],
  rows: [
    new TableRow({ tableHeader: true, children: [
      headerCell("Editorial Requirement", 2680),
      headerCell("How the Revision Addresses It", 4360),
      headerCell("Location", 2320),
    ]}),
    new TableRow({ children: [
      bodyCell(cellParas(["Results accurately reported"]), 2680),
      bodyCell(cellParas([
        "The simulation’s estimation behavior is now stated explicitly and the bias reported candidly. §3.2 specifies the full model fitted by SAR maximum likelihood (both x_i and z_i included) and explains the ≈−ψ bias that arises under x_i–z_i collinearity (see C4); §5.2 now quantifies the magnitude of that bias (50–80% of β) rather than leaving it implicit (see C6).",
      ]), 4360),
      bodyCell(cellParas(["§3.2 (p. 17); §5.2 (p. 31)"]), 2320),
    ]}),
    new TableRow({ children: [
      bodyCell(cellParas(["Overstated conclusions rewritten"]), 2680),
      bodyCell(cellParas([
        "The blanket endorsement of BSS has been replaced with a conditional recommendation. §4.4 now recommends BSS only when spillover is expected to flow primarily into control clusters, and explicitly names the regime (bidirectional spillover with spatial dependence) in which SRS or a buffer-zone design is preferable; the §5.1 “checkerboard paradox” discussion reinforces this (see C6).",
      ]), 4360),
      bodyCell(cellParas(["§4.4 (p. 26); §5.1 (p. 29)"]), 2320),
    ]}),
    new TableRow({ children: [
      bodyCell(cellParas(["Limitations fully explained"]), 2680),
      bodyCell(cellParas([
        "§5.2 now fully develops the identifiability limitation (C2/C4), the inferential-scope limitation (C5), and the severity of the bias (C6); §5.3 adds buffer-zone (“fried-egg”) designs as the analysis-stage remedy that addresses the limitation at its root.",
      ]), 4360),
      bodyCell(cellParas(["§5.2 (p. 31); §5.3 (p. 34)"]), 2320),
    ]}),
    new TableRow({ children: [
      bodyCell(cellParas(["Manuscript length reduced"]), 2680),
      bodyCell(cellParas([
        "Two structural condensations shorten the article with no loss of results, figures, or references: the three full result tables (Tables 1\u20133) are relocated to a new Appendix E, and the Introduction\u2019s method-catalog is condensed into a concise two-category summary. The Results-section narrative is correspondingly tighter (see the Manuscript Length section).",
      ]), 4360),
      bodyCell(cellParas(["\u00a73.4\u20133.6 & App. E (p. 45); \u00a71 (p. 3)"]), 2320),
    ]}),
  ],
}));

// ============================================================================
// REVIEWER COMMENTS
// ============================================================================
children.push(sectionHeading("REVIEWER COMMENTS"));

// ---- C1 ----
children.push(itemHeading("C1 — Estimands and parameter interpretation (Major)"));
children.push(itemTable(
  [
    "“The estimands are not clear in the article. The authors have identified ρ as a nuisance parameter and φ as the causal design level quantity, which is referred to as the indirect effect. However, in the autoregressive model I do not believe that φ represents the indirect effect. For example, if z_j is cluster j’s allocation and z_j′ the allocation of the other clusters with potential outcome Y_i(z_j, z_j′) then we might specify an indirect effect as something like E[Y(0,1)] − E[Y(0,0)]. If the autoregressive model is the correct model, then this quantity would include a contribution of ρ through the spatial autoregression as well as φ. Indeed, with ρ ≠ 0 ‘control-only’ spillover still includes a spillover mechanism from treated to treated clusters. It would help to be explicit about what the treatment effects are, what the parameters mean, and if the model identifies them. In other places β is just referred to as ‘the treatment effect’. In Section 2.3.1 it notes ‘… β̂ is more precisely interpreted as a total treatment effect.’ This doesn’t seem to be correct in the context of the model where β would be only the direct effect. Indeed (6) provides the ‘equilibrium’ state of the spatially autoregressive model and it’s clear that the expectation of any individual y_ik is a function of the whole mean vector and parameter ρ and that E[Y(0,1)] − E[Y(0,0)] does not equal any single parameter in the model.”",
  ],
  [
    "We agree on both points and have substantially revised the treatment of estimands. The misleading sentence describing β̂ as a “total treatment effect” has been removed; §2.3.1 now defines β as the **direct structural coefficient** of the intervention — the change in a cluster’s own outcome from its own receipt of treatment, holding z_i and the spatial lag fixed.",
    "More fundamentally, a new paragraph in §2.5 states explicitly that β, ψ, and ρ are **structural coefficients** of the SAR model, not marginal potential-outcome contrasts. We give the model’s reduced (equilibrium) form, y = (I − ρW)^(−1)(α + Xβ + Zψ + ε) [new Eq. (4)], and use it to make precisely the Reviewer’s point: because a change in any cluster’s status propagates through the spatial multiplier (I − ρW)^(−1), a potential-outcome indirect effect such as E[Y(0,1)] − E[Y(0,0)] is a *derived* function of β, ψ, ρ, and W, and equals no single coefficient in the model.",
    "We further note that when ρ ≠ 0, even “control-only” spillover induces propagation among treated clusters through this multiplier, so ψ must be read as the structural spillover coefficient on z_i rather than as the marginal indirect effect. We have accordingly removed the word “causal” from the description of ψ. The simulation is now framed unambiguously as evaluating how well each design recovers the structural coefficient β.",
  ],
  [
    "§2.3.1 (p. 7); β/ψ/ρ estimands paragraph and reduced-form Eq. (4) in §2.5 (pp. 12–14); scope note in §5.2 (p. 31).",
    "*Starts with:* “The quantities β, ψ, and ρ are structural coefficients of the spatial autoregressive model rather than marginal potential-outcome contrasts …”",
  ],
));

// ---- C2 ----
children.push(itemHeading("C2 — Spillover specification and identifiability (Major)"));
children.push(itemTable(
  [
    "“The spillover effect is unclear – the authors reference both distance and neighbour functions. The spillover effect is written as φ(d_i) where d_i is ‘distance from the nearest intervention source’, but an intervention source is not defined. Is it the boundary of the cluster or another treated unit? If it’s at individual level, then the spillover is heterogeneous across individuals. But then later it’s assumed to be neighbour 1/0 only, but this would be at the level of the cluster only. Then in the model (6) φ is now a parameter and not a function of anything. If so, then for almost all allocations on the small grids the treatment indicator and spillover indicator will be exactly collinear and the effect not identifiable. It’s also unclear because in the notation they vary at the individual level. Similarly, the study observation locations (e.g. Fig 4) are simulated to have random locations but if everything is at cluster level and not distance based then this location information is discarded, right?”",
  ],
  [
    "We agree the φ(d_i) notation was inconsistent and have removed it entirely. §2.3.2 now states that, although a spillover effect can in general be specified as a function of distance to the nearest intervention cluster or of the number of treated neighbors, this study uses a single **binary cluster-level exposure**: z_i in {0,1} equals 1 when cluster i shares a Rook edge with at least one intervention cluster (and the spillover type permits it). Because z_i is a 0/1 indicator, ψ is a scalar multiplier, not a function; and because z_i is cluster-level, spillover does not vary across individuals within a cluster. The undefined “intervention source” and “distance from the nearest intervention source” language has been removed throughout.",
    "Crucially, we now **concede the identifiability point directly**: §2.3.2 states that on these small grids z_i can be exactly collinear with x_i for some allocations, in which case β and ψ are not separately identifiable from the model fit alone; the consequences are examined in §3 and §5.2.",
    "Finally, we confirm the Reviewer’s inference about the simulated coordinates. §3.2 now states that subject coordinates are used only to assign cluster membership and to render the grid figure and **do not enter the SAR fit** — the fitted model depends on each subject only through which cluster it belongs to, and within-cluster variation contributes only to the error variance, carrying no information about β, ψ, or ρ. Subject-level/location-based spillover is identified as future work.",
  ],
  [
    "Spillover/identifiability paragraphs in §2.3.2 (p. 8); coordinate-role clarification in §3.2 (p. 17).",
    "*Starts with:* “A consequence of representing spillover through a binary cluster-level indicator is that … z_i can be strongly correlated with the treatment indicator x_i …”",
  ],
));

// ---- C3 ----
children.push(itemHeading("C3 — Notation (Major)"));
children.push(itemTable(
  [
    "“There are several issues with the notation.",
    "a. (1) w_ij is not defined.",
    "b. (1) why is y indexed by both i and k but not x or ε? (3) x and z have i,k indexes.",
    "c. (4)–(6) the cluster index is dropped altogether now, why?",
    "d. The matrix W – how is ‘neighbour’ defined here? The entries appear to be at individual level but the neighbour status is at cluster level, and if the rows are row-normalised then really this is a spillover of the cluster mean from cluster to cluster, not individual level observations.",
    "e. φ is used both as a function and as a scalar.",
    "f. (4) the indexing over N seems to be redundant since w_ij is zero if it’s not in the set of neighbours.",
    "g. (5) contains a different expression than (4).”",
  ],
  [
    "All seven items have been corrected.",
    "**(a)** w_ij is now defined inline immediately after Eq. (1): the elements of the row-normalized spatial weights matrix W, equal to 1/|N(i)| for each Rook neighbor and 0 otherwise (full construction in Appendix A).",
    "**(b)** Indexing is harmonized: in Eq. (1) ε_i → ε_ik to match y_ik; in Eq. (3) x_ik, z_ik → x_i, z_i, since treatment and spillover status are cluster-level quantities shared by all subjects in a cluster, while outcomes and errors carry the subject index ik.",
    "**(c)** The cluster index is no longer dropped; the ik (and j,k′) indexing is carried through Eqs. (1)–(6) and the appendix likelihood. Within a cluster, because x_i and z_i are constant, individual outcomes differ only through ε_ik.",
    "**(d)** W is now stated explicitly to be a **cluster-level n×n adjacency matrix**, with i, j in w_ij referring to clusters, not individuals. We acknowledge the Reviewer’s reading is correct: with row-normalized cluster-level weights the spatial lag is effectively the spillover of cluster-level mean outcomes from cluster to cluster, and the text now says so.",
    "**(e)** ψ is a scalar throughout; the φ(d_i) function notation has been removed (see C2).",
    "**(f)** The redundant set notation is removed: the spatial sum is now written Σ_{j≠i} w_ij y_jk′, since w_ij = 0 for non-neighbors already restricts the sum to neighbors.",
    "**(g)** The spatial-dependence term and the estimation equation now use the identical form Σ_{j≠i} w_ij y_jk′, resolving the discrepancy between the two expressions.",
  ],
  [
    "SAR model and notation block, §2.5, Eqs. (1)–(6) (pp. 12–15); spatial-dependence term §2.6 (p. 14); weights construction Appendix A (p. 42).",
    "*Starts with:* “where N(i) denotes the set of clusters that share a Rook edge with cluster i …”",
  ],
));

// ---- C4 ----
children.push(itemHeading("C4 — The estimated model (Major)"));
children.push(itemTable(
  [
    "“The simulation study remains a little unclear to me. What model is being estimated for the simulation study? While (6) gives the DGP, the results suggest that the models do not include the spillover indicator, but potentially do include the autoregressive component. As per point 2, the control-only spillover bias for block stratified sampling is effectively the value of φ, which is what you’d expect with perfect collinearity and omitted variable bias. The estimated model should be explained.”",
  ],
  [
    "We have made the estimated model explicit. A new paragraph in §3.2 states that all four structural parameters (α, β, ψ, ρ) are estimated jointly by SAR maximum likelihood (lagsarlm()), with **both** the intervention indicator x_i **and** the spillover indicator z_i included alongside the spatial lag — Eq. (6). The DGP and the estimation model therefore share the same structure; the bias does not arise from omitting z_i.",
    "We then explain the mechanism the Reviewer correctly identifies: in scenarios where spillover propagates only to adjacent control clusters, z_i is determined by adjacency to treated clusters and can become exactly collinear with x_i; when that collinearity is present, the jointly estimated β̂ absorbs part of the spillover signal and exhibits a systematic negative bias of approximately −ψ. This is exactly the omitted-variable-like pattern the Reviewer describes — here arising from collinearity *within* the fitted model rather than from omission — and it matches the bias visible in the control-clusters-only columns of Tables 1–3.",
    "We frame this as a structural property of the small-grid geometry (a limited feasible-allocation set), not a defect of the estimator (§5.2).",
  ],
  [
    "Estimation-model paragraph in §3.2 (p. 17); identifiability framing in §5.2 (p. 31); Eq. (6); Tables 1–3 (now Appendix E, p. 45).",
    "*Starts with:* “The structural parameters (α, β, ψ, ρ) were estimated jointly for each replication via SAR maximum likelihood, with both the intervention indicator x_i and the spillover indicator z_i included …”",
  ],
));

// ---- C5 ----
children.push(itemHeading("C5 — Frequentist inference under near-deterministic BSS (Major)"));
children.push(itemTable(
  [
    "“In my first round comments I queried whether the lack of permutations of the treatment allocation made sense in the context of a randomised controlled trial. The authors’ response was that ‘BSS is a design-stage tool evaluated through exhaustive enumeration, not a randomization distribution for permutation inference.’ But this does not address the issue for Frequentist inference which relies on interpretation around the sampling distribution under the null. The sampling distribution for BSS seems to have only 1 or 2 values and the treatment allocation would be almost deterministic. If it’s only a design-stage tool, what purpose does it serve if it’s not the actual allocation mechanism?”",
  ],
  [
    "We take this point seriously and address it in two parts — conceding the inference limitation and clarifying the design rationale.",
    "*(Concede.)* We now state plainly (§3.1; §5.2) that because the block-stratification condition admits only 2–6 valid allocations, BSS does **not** generate a usable randomization distribution and does **not** support permutation- or randomization-based inference; we make no such claim.",
    "*(Clarify.)* Frequentist inference for the intervention effect does not, however, require a randomization distribution over allocations. For any single realized allocation — SRS or BSS — β̂ has a well-defined sampling distribution arising from the stochastic outcome process, and **model-based SAR maximum likelihood (lagsarlm())** conditions on the realized allocation and derives that distribution from the error term. This is the appropriate inferential framework for a practitioner implementing either design, and the small size of the BSS allocation set does not invalidate it; the replication-to-replication variation in our simulation is correspondingly outcome-level Monte Carlo variation for a fixed allocation, not re-randomization of the allocation.",
    "On the Reviewer’s final question — what purpose BSS serves if it is not the allocation mechanism — we clarify that **BSS is the allocation mechanism**: an investigator adopting BSS deliberately selects one of the few spatially isolated allocations rather than drawing at random. Its near-determinism is the point: by excluding the high-error allocations SRS can draw by chance, BSS bounds the worst-case estimation error across the feasible allocations. The relevant design question is therefore not how many allocations a method admits, but whether the principled fixed allocation it selects yields a low-error estimator — which the exhaustive MSE comparison answers directly.",
  ],
  [
    "Inferential-scope block in §3.1 (p. 16); inference-constraint paragraph in §5.2 (p. 31).",
    "*Starts with:* “Block stratified sampling is itself the treatment allocation mechanism …”",
  ],
));

// ---- C6 ----
children.push(itemHeading("C6 — Severity of bias and buffer-zone designs (Major)"));
children.push(itemTable(
  [
    "“Interpretation. The results suggest that with the models described here the bias might be potentially catastrophic and so none of these designs should be used. BSS might be slightly better in some scenarios but it’s so severe that ‘fried egg’ or ‘buffer zone’ designs should really be the preferred approach.”",
    "[Editor’s annotation: “‘It’s so severe’??? – what is severe”]",
  ],
  [
    "We agree the magnitude of the bias must be stated candidly, and we have done so while also clarifying that it is not unique to the designs studied. §5.2 now **quantifies the severity directly** (this also answers the Editor’s annotation): with β = 1.0 and ψ in [0.5, 0.8], an exactly collinear allocation yields β̂ ≈ β − ψ, a bias of approximately −ψ — i.e., 50–80% of the true effect, large enough to materially understate or reverse a practical conclusion.",
    "We place this in context with two points. **First**, this error is almost entirely *bias* arising from lost identifiability (the exact x_i–z_i confounding), not *variance*; it afflicts any allocation that surrounds every control cluster with treated neighbors, and the equally collinear SRS allocations exhibit comparable bias — so it is not a property of block stratification specifically. **Second**, where the allocation instead permits z_i to be estimated separately (as in many SRS allocations under bidirectional spillover), the bias of β̂ falls to roughly 1–2% of β.",
    "Accordingly: **(i)** the BSS recommendation is now *conditional*, not universal (§4.4) — recommending BSS only when spillover is expected to flow primarily into control clusters, and pointing to SRS or buffer-zone designs when bidirectional spillover with dependence is anticipated; and **(ii)** we have added buffer-zone (“fried-egg”) designs to Future Considerations (§5.3, citing Jarvis et al. 2017) as the analysis-stage remedy that addresses identifiability at its root — excluding contaminated boundary controls restores z_i = 0 and breaks the x_i–z_i collinearity by construction.",
    "We note candidly that buffer-zone exclusion is infeasible on the compact grids studied here (every control borders an intervention cluster, so exclusion would remove the entire control arm) and therefore requires larger layouts; comparing SRS, BSS, and buffer-zone designs on such layouts is identified as the key next step. We respectfully retain the paper’s core contribution — characterizing the comparative bias–variance behavior of SRS vs. BSS across grid geometries and spillover regimes — because that comparison is itself the actionable design guidance, and it points investigators toward buffer-zone designs precisely where they are needed.",
  ],
  [
    "Severity quantification in §5.2 (p. 31); conditional recommendation in §4.4 (p. 26); buffer-zone paragraph in §5.3 (p. 34).",
    "*Starts with:* “The magnitude of this bias warrants candid discussion. With a true intervention effect of β = 1.0 and spillover magnitudes ψ in [0.5, 0.8] …”",
  ],
));

// ============================================================================
// SUMMARY TABLE
// ============================================================================
children.push(sectionHeading("MANUSCRIPT LENGTH"));
children.push(para(rich(
  "The Editor requested that the manuscript\u2019s length be reduced. We have done so through two structural condensations that shorten the article without removing any results, figures, or references. The condensed manuscript is **Revisions_V2b_short.pdf**; a cumulative annotated version showing all changes relative to the original submission is provided in **Revisions_Annotated_V2b_short.pdf**."
), { after: 120 }));
children.push(para(rich(
  "**1. Result tables relocated to an appendix.** The three full result tables (Tables 1\u20133), which previously occupied substantial space within the Results section (\u00a73.4\u20133.6), have been moved to a new **Appendix E** (p. 45). The three grouped-bar MSE figures remain inline in \u00a73.4\u20133.6, so the Results narrative retains its visual evidence while the dense numeric tables are consolidated in the appendix. We verified that every table data row is identical to the full revision \u2014 no results were altered."
), { after: 120 }));
children.push(para(rich(
  "**2. Introduction method-catalog condensed.** The review of prior spatial-adjustment methods in the Introduction (\u00a71, p. 3) has been condensed from a study-by-study enumeration into a concise two-category summary (spatial-variable vs. spatial-modeling approaches). All cited works are retained \u2014 the bibliography is unchanged \u2014 and the downstream discussion of specific studies (e.g., Silcocks & Kendrick, Chao et al.) is preserved."
), { after: 120 }));
children.push(para(rich(
  "These changes reduce the length of the Results section body and tighten the Introduction. The condensed manuscript is otherwise identical in scientific content to the full revision: every response to comments C1\u2013C6 above applies verbatim, and the section and equation structure is unchanged apart from the addition of Appendix E."
), { after: 60 }));

children.push(sectionHeading("SUMMARY OF ALL CHANGES"));
children.push(para(rich("The table below consolidates all revisions made for this round, with locations in Revisions_V2b_short.pdf."), { after: 120 }));

const SUMC = [1100, 6060, 2200];
function sumRow(id, summary, loc) {
  return new TableRow({ children: [
    bodyCell(cellParas([`**${id}**`]), SUMC[0]),
    bodyCell(cellParas([summary]), SUMC[1]),
    bodyCell(cellParas([loc]), SUMC[2]),
  ]});
}
children.push(new Table({
  width: { size: CONTENT, type: WidthType.DXA }, columnWidths: SUMC,
  rows: [
    new TableRow({ tableHeader: true, children: [
      headerCell("ID", SUMC[0]), headerCell("Change Summary", SUMC[1]), headerCell("Location (Revisions_V2b_short)", SUMC[2]),
    ]}),
    sumRow("Editor", "Results reported accurately, overstated conclusions rewritten, limitations fully explained — realized through the C4, C6, C5, and C2 revisions below.", "§3.2, §4.4, §5.1–5.3"),
    sumRow("Length", "Reduced manuscript length per the Editor’s request: relocated result Tables 1–3 to new Appendix E (MSE figures retained inline; all data rows verified unchanged) and condensed the Introduction method-catalog (all citations retained, bibliography unchanged).", "§3.4–3.6, App. E (p. 45); §1 (p. 3)"),
    sumRow("C1", "Removed “total treatment effect” claim; β redefined as the direct structural coefficient; added estimands paragraph and reduced-form Eq. (4) showing E[Y(0,1)]−E[Y(0,0)] is a derived quantity, not any single parameter; “causal” removed from ψ.", "§2.3.1 (p. 7); §2.5 + Eq. (4) (pp. 12–14); §5.2 (p. 31)"),
    sumRow("C2", "Removed φ(d_i)/distance notation; spillover defined as binary cluster-level indicator z_i with scalar ψ; conceded exact x_i–z_i collinearity / non-identifiability; clarified subject coordinates do not enter the SAR fit.", "§2.3.2 (p. 8); §3.2 (p. 17)"),
    sumRow("C3", "Fixed all seven notation items (a–g): defined w_ij; harmonized ik / cluster indexing; carried cluster index through Eqs. (1)–(6); declared W cluster-level n×n; scalar ψ; removed redundant set indexing; matched Eqs. (5)–(6) summation form.", "§2.5–2.6, Eqs. (1)–(6) (pp. 12–15); App. A (p. 42)"),
    sumRow("C4", "Stated the estimated model explicitly: α, β, ψ, ρ jointly estimated via SAR MLE with both x_i and z_i included; explained the ≈−ψ bias under collinearity, matching Tables 1–3.", "§3.2 (p. 17); §5.2 (p. 31)"),
    sumRow("C5", "Conceded that BSS (2–6 allocations) supports no randomization/permutation inference; clarified that model-based SAR MLE conditions on the realized allocation and is the valid framework; reframed near-determinism as the deliberate allocation mechanism.", "§3.1 (p. 16); §5.2 (p. 31)"),
    sumRow("C6", "Quantified bias severity (50–80% of β; answers Editor’s annotation); contextualized as identifiability bias common to equally collinear SRS; made BSS recommendation conditional; added buffer-zone designs to Future Considerations.", "§5.2 (p. 31); §4.4 (p. 26); §5.3 (p. 34)"),
  ],
}));

// ============================================================================
// CLOSING
// ============================================================================
children.push(new Paragraph({ spacing: { before: 280, after: 140 }, children: rich(
  "We are grateful to the Reviewer and the Editor for pressing on these points. Addressing them has materially sharpened the manuscript’s treatment of estimands, identifiability, and inferential scope, and has produced a more honest and more useful set of design recommendations. This condensed version preserves that content in full while meeting the Editor’s request for a shorter manuscript. We hope the revised manuscript now meets the standard for publication and look forward to your assessment."
)}));
children.push(para([new TextRun("Sincerely,")], { before: 120, after: 60 }));
children.push(para([new TextRun({ text: "Andrew Walther, Tonya Van Deinse, and Feng-Chang Lin", bold: true })], { after: 20 }));
children.push(para([new TextRun("University of North Carolina at Chapel Hill")], { after: 20 }));
children.push(para([new TextRun("Department of Biostatistics / School of Social Work")], { after: 20 }));
children.push(para([new TextRun("awalther@unc.edu | flin@bios.unc.edu")], { after: 0 }));

// ============================================================================
// DOCUMENT
// ============================================================================
const doc = new Document({
  styles: { default: { document: { run: { font: "Calibri", size: 22 } } } }, // 11pt
  sections: [{
    properties: { page: {
      size: { width: 12240, height: 15840 },
      margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
    }},
    children,
  }],
});

const OUT = process.argv[2];
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(OUT, buf); console.log("WROTE", OUT, buf.length, "bytes"); });
