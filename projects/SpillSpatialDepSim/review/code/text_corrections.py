"""
Script: review/code/text_corrections.py
Purpose: Apply the post-publication text corrections (3x3/3x4 results, retired
         "checkerboard paradox", and the summary statements they invalidate) to
         the accepted LaTeX manuscript and to the Chapter 2 .qmd, from one list of
         edits so the two documents stay word-for-word in step.
Author: Andrew Walther (review written with Claude Code)
Created: 2026-10-25
Usage (from review/):
  python3 code/text_corrections.py tex <in.tex> <out.tex>
  python3 code/text_corrections.py qmd <in.qmd> <out.qmd>
  python3 code/text_corrections.py md  <manuscript-declarations.md> <out.md>   (abstract only)
Each edit's "old" text must match exactly once (whitespace-insensitive), or the
script stops. Edit IDs (T1-T13) follow notes/2026-10-25-correction-plan.md.
Every number in the new text comes from tables/beta_mse_table_<grid>_corrected.csv,
tables/optimal_combos_<grid>_corrected.csv or COMPARISON.md.
"""
import re
import sys
import textwrap

# Edits, written in the manuscript's LaTeX. Cross-references use {SEC:...},
# {TAB:...} and {FIG:...} tokens so they can be rendered for either format.
EDITS = [
    # N1: Chapter 2 only, the post-publication note (decision D4: body text, no footnote)
    ("N1", "qmd",
     r"spatial dependence levels, and spillover application types.",
     r"spatial dependence levels, and spillover application types. The 3$\times$3 and "
     r"3$\times$4 results reported here were recomputed after the article was accepted, to "
     r"correct two coding errors that affected only those grids. The 2$\times$4 results "
     r"are unchanged."),

    # T13: abstract. Chapter 2 has no abstract; its verbatim copy of the paper's abstract
    # in manuscript-declarations.md gets these edits (format "md").
    ("T13a", "abstract",
     r"but block stratified treatment assignments lead to minimal variation in MSE",
     r"but block stratified treatment assignments generally lead to less variation in MSE"),
    ("T13b", "abstract",
     r"but block stratified assignment minimizes variance among combinations",
     r"but block stratified assignment generally reduces variance among combinations"),
    ("T13c", "abstract",
     r"BSS bounds worst-case allocation risk through near-zero variation in estimation error",
     r"BSS bounds worst-case allocation risk through lower variation in estimation error"),

    # C1, C2: source comments that named the retired "checkerboard paradox" (manuscript only;
    # comments do not print, but the source should not point readers to it)
    ("C1", "tex",
     r"% R1.5: Risk-management framing ($3\times4$ grid); R2.M4 checkerboard paradox noted",
     r"% R1.5: Risk-management framing ($3\times4$ grid); rewritten in the post-publication correction"),
    ("C2", "tex",
     r"% R2.M4: Chessboard paradox note",
     r"% Collinearity of checkerboard allocations (replaces the R2.M4 note, retired in the correction)"),

    # T1: 2x4 section, effect of rho
    ("T1", "both",
     r"The presence of spatial dependence ($\rho = 0.01$) consistently increased MSE for "
     r"both methods relative to the no-dependence setting, a pattern that persists across "
     r"all grid sizes.",
     r"The presence of spatial dependence ($\rho = 0.01$) increased the MSE of SRS relative "
     r"to the no-dependence setting in every scenario, a pattern that persists across all "
     r"grid sizes. For BSS on this grid, spatial dependence lowered MSE (e.g., $0.165$ vs. "
     r"$0.262$ at $\psi = 0.5$)."),

    # T2: 3x3 results, both paragraphs
    ("T2a", "both",
     r"For the Control Clusters Only (TrtNoSpill) spillover application, BSS consistently "
     r"yielded significantly lower MSE, SD, and bias compared to SRS. For instance, with "
     r"$\rho=0$, BSS MSE was minimal (e.g., $0.001$ for $\psi=0.5$), while SRS MSE ranged "
     r"from $0.092$ to $0.235$. This pattern held true with spatial dependence "
     r"($\rho=0.01$), where BSS MSE remained very low (e.g., $0.004$ for $\psi=0.5$), "
     r"contrasting with SRS MSE values ranging from $0.259$ to $0.372$. SRS estimates "
     r"consistently showed negative bias, increasing with $\psi$, while BSS bias was "
     r"negligible.",
     r"For the Control Clusters Only (TrtNoSpill) spillover application, BSS yielded lower "
     r"MSE, SD, and bias than SRS in every scenario. With $\rho=0$, BSS MSE ranged from "
     r"$0.042$ to $0.106$, while SRS MSE ranged from $0.114$ to $0.291$. With spatial "
     r"dependence ($\rho=0.01$), BSS MSE ranged from $0.147$ to $0.240$, against $0.274$ to "
     r"$0.455$ for SRS. Both designs showed negative bias that grew in magnitude with "
     r"$\psi$, and the bias was smaller for BSS (e.g., $-0.078$ vs. $-0.226$ at $\rho=0$, "
     r"$\psi=0.5$)."),
    ("T2b", "both",
     r"In the Intervention \& Control Clusters (TrtSpill) spillover scenario, BSS continued "
     r"to demonstrate superior performance with substantially lower MSE, SD, and bias "
     r"compared to SRS. For $\rho=0$, BSS MSE was consistently around $0.0004$, whereas SRS "
     r"MSE ranged from $0.002$ to $0.005$. With $\rho=0.01$, BSS MSE was also consistently "
     r"very low (e.g., $0.039$ for $\psi=0.5$), while SRS MSE ranged from $0.103$ to "
     r"$0.168$.",
     r"In the Intervention \& Control Clusters (TrtSpill) spillover scenario, SRS achieved "
     r"lower average MSE than BSS in every scenario. For $\rho=0$, SRS MSE ranged from "
     r"$0.002$ to $0.005$, whereas BSS MSE ranged from $0.042$ to $0.106$. With "
     r"$\rho=0.01$, SRS MSE ranged from $0.103$ to $0.167$ and BSS MSE from $0.147$ to "
     r"$0.240$. Unlike the 2$\times$4 and 3$\times$4 grids, where the BSS allocations are "
     r"the two checkerboards, the six BSS allocations on this grid differ in structure, so "
     r"their MSE also varies across allocations. Under TrtSpill, the SD of BSS was larger "
     r"than that of SRS."),

    # T3: 3x4 results, all four paragraphs
    ("T3a", "both",
     r"These results reveal both the typical strengths of BSS and a notable exception that "
     r"illuminates the limits of the checkerboard design.",
     r"These results follow the same pattern as the 2$\times$4 grid in every scenario."),
    ("T3b", "both",
     r"Under the Control Clusters Only (TrtNoSpill) spillover condition, BSS dramatically "
     r"outperforms SRS across all combinations of $\rho$ and $\psi$. With no spatial "
     r"dependence ($\rho = 0$), BSS MSE is an order of magnitude smaller than SRS (e.g., "
     r"BSS MSE $\approx 0.0004$ vs. SRS MSE $= 0.034$ at $\psi = 0.5$). With spatial "
     r"dependence ($\rho = 0.01$), BSS retains a substantial advantage (e.g., BSS MSE "
     r"$\approx 0.035$ vs. SRS MSE $= 0.120$ at $\psi = 0.5$). BSS SD is near zero "
     r"throughout, while SRS SD is substantially larger, confirming that the two valid BSS "
     r"checkerboard allocations for this grid provide stable and consistently low-error "
     r"estimates whenever spillover is confined to control clusters.",
     r"Under the Control Clusters Only (TrtNoSpill) spillover condition, performance again "
     r"depends on the level of spatial dependence. With spatial dependence ($\rho = 0.01$), "
     r"BSS achieves lower MSE than SRS (e.g., BSS MSE $= 0.133$ vs. SRS MSE $= 0.244$ at "
     r"$\psi = 0.5$). With no spatial dependence ($\rho = 0$), SRS achieves lower average "
     r"MSE (e.g., SRS MSE $= 0.133$ vs. BSS MSE $= 0.226$ at $\psi = 0.5$). The SD of BSS "
     r"is smaller than that of SRS throughout, because the two valid BSS checkerboard "
     r"allocations for this grid give similar estimation error."),
    ("T3c", "both",
     r"Under the Intervention \& Control Clusters (TrtSpill) spillover condition, results "
     r"are more nuanced. With no spatial dependence ($\rho = 0$), BSS achieves lower MSE "
     r"than SRS (e.g., BSS MSE $\approx 0.0003$ vs. SRS MSE $= 0.001$ at $\psi = 0.5$). "
     r"However, when spatial dependence is present ($\rho = 0.01$), BSS MSE ($\approx "
     r"0.223$ at $\psi = 0.5$) substantially exceeds SRS MSE ($= 0.112$). This inversion "
     r"occurs because the strict checkerboard structure of the 3$\times$4 BSS allocation "
     r"ensures that every intervention cluster is surrounded entirely by control neighbors "
     r"eligible for spillover; under TrtSpill, this configuration amplifies the systematic "
     r"bias in $\hat{\beta}$, a ``checkerboard paradox\" discussed further in {SEC:Discussion}.",
     r"Under the Intervention \& Control Clusters (TrtSpill) spillover condition, SRS "
     r"achieves lower average MSE in every scenario (e.g., SRS MSE $= 0.112$ vs. BSS MSE "
     r"$= 0.133$ at $\psi = 0.5$, $\rho = 0.01$, and $0.001$ vs. $0.226$ at $\rho = 0$). In "
     r"both checkerboard allocations, every control cluster borders an intervention "
     r"cluster and no intervention cluster borders another, so the spillover indicator "
     r"equals $1 - x_i$ and the two spillover types produce identical BSS results. "
     r"$\beta$ and $\psi$ are then not separately identifiable, and $\hat{\beta}$ absorbs "
     r"the spillover effect, with bias near $-\psi$ when $\rho = 0$ (e.g., $-0.463$ at "
     r"$\psi = 0.5$; {SEC:Limitations})."),
    ("T3d", "both",
     r"The 3$\times$4 results reinforce the broader conclusion: BSS reliably bounds "
     r"worst-case estimation error under TrtNoSpill, while the TrtSpill regime with spatial "
     r"dependence represents the primary scenario in which a practitioner should carefully "
     r"weigh the risk-management benefits of BSS against its potential for amplified bias.",
     r"The 3$\times$4 results reinforce the broader conclusion: BSS bounds worst-case "
     r"estimation error and has its advantage under TrtNoSpill with spatial dependence, "
     r"while SRS offers lower average MSE when spillover reaches both treatment arms."),

    # T4: Recommendations, qualification built on the paradox
    ("T4", "both",
     r"If substantial bidirectional spillover among treated units is instead anticipated "
     r"together with background spatial dependence, the spatial separation enforced by BSS "
     r"can amplify rather than reduce bias (the ``checkerboard paradox,\" {SEC:Discussion}); "
     r"in that regime",
     r"If substantial spillover among treated units is instead anticipated, SRS achieved "
     r"lower average MSE in every simulated scenario of that type ({SEC:Discussion}). In "
     r"that regime"),

    # T5: Discussion, bias/variance trade-off
    ("T5", "both",
     r"The results indicate that randomly selected treatment assignments lead to reduced "
     r"average error (lower bias) when estimating intervention effects, while block "
     r"stratified treatment assignments consistently result in reduced variation (lower "
     r"variance) among a series of estimates.",
     r"On the 2$\times$4 and 3$\times$4 grids, randomly selected treatment assignments led "
     r"to lower bias when estimating intervention effects, while block stratified treatment "
     r"assignments resulted in less variation in MSE across allocations in most scenarios. "
     r"On the 3$\times$3 grid, whose six BSS allocations are not all checkerboards, BSS had "
     r"the smaller bias in most scenarios, and under TrtSpill its MSE varied more across "
     r"allocations than that of SRS."),

    # T6: Discussion, effect of rho
    ("T6", "both",
     r"The presence of even a small degree of spatial dependence (e.g., $\rho=0.01$) "
     r"generally increased MSE for both sampling methods.",
     r"The presence of even a small degree of spatial dependence (e.g., $\rho=0.01$) "
     r"increased MSE for SRS in every scenario. For BSS, the effect depended on the grid: "
     r"MSE rose on the 3$\times$3 grid but fell for the checkerboard allocations of the "
     r"2$\times$4 and 3$\times$4 grids."),

    # T7: Discussion, spillover types
    ("T7", "both",
     r"When spillover was restricted to control clusters only, the MSE for the intervention "
     r"effect was higher for SRS compared to BSS, particularly with increasing spillover "
     r"magnitude. Conversely, when spillover was permitted for both control and "
     r"intervention clusters, SRS often exhibited lower MSE, especially when spatial "
     r"dependence was absent.",
     r"When spillover was restricted to control clusters only and spatial dependence was "
     r"present, the MSE for the intervention effect was higher for SRS than for BSS on "
     r"every grid. Without spatial dependence, SRS had the lower MSE on the 2$\times$4 and "
     r"3$\times$4 grids and BSS on the 3$\times$3 grid. Conversely, when spillover was "
     r"permitted for both control and intervention clusters, SRS exhibited lower MSE in "
     r"every scenario, especially when spatial dependence was absent."),

    # T8: Discussion, optimal allocations and "minimal MSE"
    ("T8a", "both",
     r"there are 7 different treatment assignment combinations, each of which",
     r"there are 7, 4, and 6 different treatment assignment combinations, respectively, "
     r"each of which"),
    ("T8b", "both",
     r"This serves as further evidence to suggest that it is prudent to recommend the Block "
     r"Stratified Sampling strategy because it provides results that consistently exhibit "
     r"minimal MSE of the intervention effect estimate although the BSS treatment "
     r"combinations are never the absolute best case in terms of MSE for any given "
     r"simulation scenario. By implementing a BSS strategy, investigators will consistently "
     r"obtain strong estimates of the intervention effect (although not the best for any "
     r"specific single case) under the conditions of the spatial autoregressive (SAR) model.",
     r"The BSS treatment combinations are never the absolute best case in terms of MSE for "
     r"any given simulation scenario, and BSS does not give the lowest average MSE in "
     r"general. Its benefit is that it excludes the poorly balanced allocations that SRS "
     r"can draw, which limits the worst case, most clearly when spillover flows into "
     r"control clusters and spatial dependence is present."),

    # T9: Discussion, the "checkerboard paradox" paragraph (retired)
    ("T9", "both",
     r"An important nuance observed in the 3$\times$4 simulation results warrants "
     r"attention in the context of practical design recommendations. When BSS forces an "
     r"exact checkerboard arrangement, maximally separating intervention and control "
     r"clusters, every intervention cluster is surrounded exclusively by control neighbors. "
     r"Under TrtSpill, this means spillover flows from each intervention cluster into its "
     r"control neighbors \emph{and} back into neighboring intervention clusters (which are "
     r"diagonally adjacent but not Rook neighbors). This structural property creates a "
     r"systematic pattern of amplified spillover-induced bias specifically for BSS in the "
     r"TrtSpill regime with spatial dependence, producing estimation error that exceeds "
     r"that of SRS in that scenario. We refer to this as the ``checkerboard paradox\": the "
     r"very design feature that makes BSS effective at minimizing spillover contamination "
     r"under TrtNoSpill, complete spatial separation of treatment arms, becomes a "
     r"liability when bidirectional spillover is expected.",
     r"Strict checkerboard allocations, which are the BSS allocations on the 2$\times$4 and "
     r"3$\times$4 grids, warrant a further note in the context of practical design "
     r"recommendations. In a checkerboard, every control cluster borders an intervention "
     r"cluster and no intervention cluster borders another, so the spillover indicator "
     r"equals $1 - x_i$ under both spillover types. The two spillover types then produce "
     r"identical data, $\beta$ and $\psi$ are not separately identifiable, and "
     r"$\hat{\beta}$ absorbs the spillover effect, with bias near $-\psi$ when $\rho = 0$ "
     r"({SEC:Limitations}). This is why the BSS results on these grids are identical under "
     r"the two spillover types."),

    # T10: Discussion summary
    ("T10", "both",
     r"an unbiased average estimate (favored by SRS) versus the need for consistent, "
     r"low-variance estimates (achieved with BSS)",
     r"an unbiased average estimate (generally favored by SRS) versus the need for "
     r"consistent, low-variance estimates (generally achieved with BSS)"),

    # T12: Limitations, why the collinear error is bias
    ("T12", "both",
     r"because the standard deviation of the block-stratified estimator is near zero, this "
     r"error is almost entirely",
     r"because the variance of these collinear estimates is small relative to their squared "
     r"bias, this error is almost entirely"),

    # T11: Future Considerations, sentence built on the paradox
    ("T11", "both",
     r"especially in the bidirectional-spillover-with-dependence scenarios where BSS "
     r"exhibited amplified bias,",
     r"especially for allocations in which treatment and spillover exposure are collinear,"),
]

# Cross-reference tokens -> (LaTeX, Quarto)
REFS = {
    "{SEC:Discussion}": (r"Section~\ref{section:Discussion}", "@sec-discussion"),
    "{SEC:Limitations}": (r"Section~\ref{section:Limitations}", "@sec-limitations"),
}


def unescape(s):
    """The raw strings above write a closing straight quote as \\" ; drop the backslash."""
    return s.replace('\\"', '"')


def to_tex(s):
    """Render tokens for the LaTeX manuscript."""
    s = unescape(s)
    for k, (tex, _) in REFS.items():
        s = s.replace(k, tex)
    return s


def to_qmd(s):
    """Render an edit string in the Chapter 2 .qmd's markup."""
    s = unescape(s)
    for k, (_, qmd) in REFS.items():
        s = s.replace(k, qmd)
    s = re.sub(r"(\d)\$\\times\$(\d)", "\\1×\\2", s)   # 3$\times$4 -> 3×4
    s = s.replace(r"\&", "&").replace("``", "“")
    s = re.sub(r"\\emph\{([^}]*)\}", r"*\1*", s)
    return s


def pattern(old, fmt):
    """Whitespace-insensitive regex for an old passage; Quarto refs may be bracketed."""
    parts = [re.escape(w) for w in old.split()]
    rx = r"\s+".join(parts)
    if fmt == "qmd":
        rx = re.sub(r"@sec\\-([a-z]+(?:\\\\-[a-z]+)*)", r"\\[?@sec\\-\1\\]?", rx)
    return re.compile(rx)


def rewrap(text, pos):
    """Re-wrap the blank-line-delimited prose paragraph containing position pos."""
    a = text.rfind("\n\n", 0, pos)
    a = 0 if a < 0 else a + 2
    b = text.find("\n\n", pos)
    b = len(text) if b < 0 else b
    para = text[a:b]
    if re.match(r"[#!|<\\:-]", para):   # not plain prose: leave as is
        return text
    wrapped = "\n".join(textwrap.wrap(" ".join(para.split()), 76,
                                      break_long_words=False, break_on_hyphens=False))
    return text[:a] + wrapped + text[b:]


def apply(text, fmt):
    for eid, scope, old, new in EDITS:
        applies = {"both": ("tex", "qmd"), "qmd": ("qmd",), "abstract": ("tex", "md"),
                   "tex": ("tex",)}[scope]
        if fmt not in applies:
            continue
        render = to_tex if fmt == "tex" else to_qmd   # "md" uses the Quarto rendering
        old_r, new_r = render(old), render(new)
        hits = pattern(old_r, "qmd" if fmt == "md" else fmt).findall(text)
        if len(hits) != 1:
            sys.exit(f"{eid}: old text matched {len(hits)} times in the {fmt} file")
        m = pattern(old_r, "qmd" if fmt == "md" else fmt).search(text)
        text = text[:m.start()] + new_r + text[m.end():]
        if fmt in ("qmd", "md"):   # re-wrap the edited paragraph to the chapter's ~76 columns
            text = rewrap(text, m.start())
        print(f"{eid}: applied")
    return text


if __name__ == "__main__":
    fmt, src, dst = sys.argv[1:4]
    out = apply(open(src).read(), fmt)
    open(dst, "w").write(out)
