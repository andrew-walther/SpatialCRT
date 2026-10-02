# ============================================================
# Script: render_real_sud_companion.R
# Purpose: Build an offline visual companion from frozen inputs and current results.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud modules, ggplot2, sf
# ============================================================
real_sud_define_only <- TRUE
.rs_report_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.rs_report_dir, "run_real_sud.R"))

#' Encode a reporting data frame as a JSON array without an additional dependency
#' @param d Data frame containing atomic character, logical or numeric columns.
#' @return One JSON string; non-finite numbers are null.
#' @family real_sud_report
#' @seealso rs_companion
#' @examples
#' rs_json(data.frame(x = 1, label = "SRS"))
rs_json <- function(d) {
  #' Encode one atomic value with JSON-compatible missingness
  #' @param x One atomic scalar.
  #' @return JSON scalar string.
  #' @examples
  #' # scalar(NA_real_) returns "null"
  scalar <- function(x) {
    if (is.na(x) || (is.numeric(x) && !is.finite(x))) "null" else
      if (is.logical(x)) tolower(as.character(x)) else if (is.numeric(x)) format(x, digits = 16, scientific = FALSE, trim = TRUE) else encodeString(x, quote = '"')
  }
  rows <- vapply(seq_len(nrow(d)), function(i) {
    values <- vapply(d, function(x) scalar(x[i]), character(1))
    paste0("{", paste(paste0(encodeString(names(d), quote = '"'), ":", values), collapse = ","), "}")
  }, character(1))
  paste0("[", paste(rows, collapse = ","), "]")
}

#' Create maps and an interactive, locally readable study review document
#' @param root Study directory containing frozen setup and optional profile results.
#' @param profile Results profile displayed; smoke and pilot are explicitly preliminary.
#' @return Path to saved HTML; maps use authorized cluster aggregates only.
#' @family real_sud_report
#' @seealso rs_reporting_rows, rs_json
#' @examples
#' # rs_companion(profile = "production")
rs_companion <- function(root = file.path(.rs_app, "results", "real_sud_rev_20261002"), profile = "pilot") {
  s <- readRDS(file.path(root, "setup.rds"))
  p <- file.path(root, profile, "performance.csv")
  results <- if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else data.frame()
  tp <- file.path(root, "tail_confirmation", "performance.csv")
  tail_results <- if (file.exists(tp) && file.exists(file.path(dirname(tp), "verification.txt")))
    read.csv(tp, stringsAsFactors = FALSE) else data.frame()
  dest <- file.path(.rs_app, "report"); assets <- file.path(dest, "real_sud_assets")
  dir.create(assets, recursive = TRUE, showWarnings = FALSE)
  # Display simplification changes only map rendering, never the cached legal W.
  g <- sf::st_simplify(s$clusters, dTolerance = 250)
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12))
  annual_maps <- do.call(rbind, lapply(as.character(2018:2021), function(y) {
    a <- g; a$Year <- y; a$Rate <- s$annual[[y]]$rate_per_100k; a
  }))
  map <- ggplot2::ggplot(annual_maps) + ggplot2::geom_sf(ggplot2::aes(fill = Rate), color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_viridis_c(name = "SUD / 100k") + ggplot2::facet_wrap(~Year, ncol = 2) +
    ggplot2::labs(title = "Observed SUD incidence: four fixed yearly planning settings", subtitle = "58 community-college clusters; denominator is population aged 18–64") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(file.path(assets, "observed_years.png"), map, width = 11, height = 7, dpi = 150)
  g$Region <- factor(s$regions)
  map <- ggplot2::ggplot(g) + ggplot2::geom_sf(ggplot2::aes(fill = Region), color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_brewer(palette = "Set2") + ggplot2::labs(title = "One frozen region partition for all four years", subtitle = "Selected by geography and balance before outcome performance is evaluated") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(file.path(assets, "regions.png"), map, width = 10, height = 4, dpi = 150)
  c <- rs_config("smoke"); units <- rs_units(c)
  examples <- do.call(rbind, lapply(1:9, function(id) {
    u <- units[which(units$Design_ID == id)[1], ]
    a <- g; z <- rs_draws(u, s, 1)[, 1]
    a$Design <- paste0(id, ". ", u$Design, " (", sum(z), "/58)")
    a$Treatment <- factor(z, levels = 0:1, labels = c("Control", "Treatment")); a
  }))
  map <- ggplot2::ggplot(examples) + ggplot2::geom_sf(ggplot2::aes(fill = Treatment), color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_manual(values = c(Control = "#e0e7ef", Treatment = "#167c80")) + ggplot2::facet_wrap(~Design, ncol = 3) +
    ggplot2::labs(title = "Reproducible example allocations on the NC geography", subtitle = "2018, queen graph, draw 1. These maps illustrate rules; they do not compare performance.") +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 9))
  ggplot2::ggsave(file.path(assets, "design_examples.png"), map, width = 14, height = 8, dpi = 150)
  status <- if (!nrow(results)) "Implementation and verification in progress; no performance results loaded." else
    paste0(if (profile != "production") "PRELIMINARY " else "PRODUCTION ", toupper(profile), ": ", nrow(results),
      " reporting rows; ", length(unique(results$Source_ID)), " distinct simulated distributions; ",
      sum(!results$Complete), " incomplete rows; ", sum(!results$Precision_OK), " rows outside precision targets.")
  evidence <- character()
  ex <- file.path(root, "exhibits")
  if (profile == "production" && file.exists(file.path(ex, "yearly_primary_design_means.csv"))) {
    a <- read.csv(file.path(ex, "yearly_primary_design_means.csv"))
    plain <- subset(a, Design_ID == 3 & Regime == "control_only")
    guided <- subset(a, Design_ID == 8 & Regime == "control_only")
    decrease <- range(100 * (1 - guided$Ratio_of_Mean_MSE_SRS))
    evidence <- c('<section><h2>What the completed primary analysis shows</h2>',
      paste0('<p>Under control-only spillover, plain saturation reduces descriptive mean MSE by ',
        formatC(100 * (1 - plain$Ratio_of_Mean_MSE_SRS[1]), digits = 1, format = "f"),
        '% relative to SRS. Incidence-guided saturation reductions range from ',
        formatC(decrease[1], digits = 1, format = "f"), ' to ', formatC(decrease[2], digits = 1, format = "f"),
        '% across the four yearly settings. These averages give equal weight to six rho/gamma pairs. Plain saturation and SRS reuse the same performance sources across years.</p>'),
      '<p>Under both-arms spillover, SRS remains competitive with balanced and saturation designs. Small differences do not establish equivalence or a uniquely optimal design. Fixed High Incidence Focus can perform well in some settings but varies across years and is sensitive to the matched baseline/adjustment specification.</p>',
      paste0('<p class="note">Primary yearly mean coverage ranges from ', formatC(100 * min(a$Coverage), digits = 1, format = "f"),
        ' to ', formatC(100 * max(a$Coverage), digits = 1, format = "f"),
        '%, below nominal 95%. Meeting Monte Carlo precision targets does not establish nominal interval coverage.</p>'),
      '<img src="../results/real_sud_rev_20261002/exhibits/primary_yearly_mse.png" alt="Yearly mean MSE for all nine designs under both spillover regimes">',
      '<img src="../results/real_sud_rev_20261002/exhibits/primary_setting_ratios.png" alt="MSE ratios against SRS in every primary parameter setting">',
      '<img src="../results/real_sud_rev_20261002/exhibits/population_shares.png" alt="Yearly treated population shares and sampled ranges">',
      '<p>High Incidence Focus treats 29 clusters but about 22–24% of the working-age population. Buffer allocations treat about 16 clusters and 27% of that population. These population shares describe geographic reach, not the proportion of healthcare professionals trained. Each yearly Spatial Blocking partition contains a singleton that remains control under retained rounding; its total budgets are 28/31/27/29.</p>',
      '<p>Queen-versus-rook sensitivity contrasts match the same rho/gamma corners. The matched beta=1 sensitivity changes both the outcome baseline and covariate adjustment, so its differences cannot be attributed to incidence causing poorer education.</p></section>')
    if (file.exists(file.path(ex, "refined_allocation_risk.png"))) evidence <- c(evidence,
      '<section><h2>Allocation downside after focused outcome refinement</h2><img src="../results/real_sud_rev_20261002/exhibits/refined_allocation_risk.png" alt="Mean and estimated 90th-percentile allocation MSE for every design and year"><p>The refined corners support lower estimated upper-tail allocation MSE for saturation designs under control-only spillover. Cross-selected outcome halves retain that descriptive pattern. Under both arms, small tail differences and unstable tail membership do not establish superior protection.</p><p class="note">Every non-singleton reporting row flags uncertain membership around the worst-decile boundary. These are the same 100 sampled allocations with more independent outcomes, not new allocation confirmation. For singleton graph and high-incidence assignments, allocation q90 equals conditional mean risk; zero variation and Jaccard=1 do not mean the estimation error or its outcome uncertainty is zero.</p></section>')
  }
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>NC SUD application companion</title>',
    '<style>body{font:17px/1.55 system-ui,sans-serif;background:#f4f6f8;color:#193247;margin:0}main{max-width:1160px;margin:auto;padding:34px 24px}h1{font-size:38px;line-height:1.15}h2{margin-top:0}section{background:white;padding:26px;border-radius:12px;margin:24px 0;border:1px solid #dce4e9}.note{border-left:5px solid #167c80;padding:14px;background:#edf6f5}.flow{display:flex;flex-wrap:wrap;gap:14px;align-items:center}.box{flex:1;min-width:200px;padding:18px;background:#edf6f5;border-radius:8px}.eq{font:22px Georgia,serif;text-align:center;padding:18px}img{width:100%;height:auto}label{display:inline-block;margin:8px 16px 8px 0}select{font:inherit;padding:5px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;border-bottom:1px solid #dde5ea;padding:8px}.scroll{overflow:auto}svg{width:100%;min-height:270px}a{color:#106b79}small{color:#536c7b}.two{display:grid;grid-template-columns:1fr 1fr;gap:20px}@media(max-width:750px){.two{grid-template-columns:1fr}h1{font-size:29px}}.badge{background:#f8edcb;padding:8px 14px;border-radius:6px}</style></head><body><main>',
    '<h1>Designing the NC SUD education study</h1><p>A visual companion to the observed-incidence application, its development and its evidence.</p>',
    paste0('<p class="badge">', status, '</p>'),
    '<section><h2>The question we agreed to study</h2><p>Which allocation designs estimate a known education intervention effect accurately on this geography, when some designs use observed SUD incidence to direct treatment?</p><div class="flow"><div class="box"><b>Fixed observed inputs</b><br>2018–2021 cluster SUD rates<br>Named geography and neighbors</div><span aria-hidden="true">→</span><div class="box"><b>Allocation design</b><br>Draw treatment/control assignments<br>Keep duplicate draw frequencies</div><span aria-hidden="true">→</span><div class="box"><b>Simulated education responses</b><br>Known τ = 1, spatial dependence<br>Spillover and independent noise</div><span aria-hidden="true">→</span><div class="box"><b>Estimate and evaluate</b><br>SAR ML estimate of τ<br>MSE, bias, coverage and allocation risk</div></div></section>',
    '<section><h2>Model and interpretation</h2><div class="eq">Y = (I − ρW)<sup>−1</sup>[τZ + S(Z) + ε]</div><div class="two"><div><b>Primary education outcome</b><p>Observed SUD incidence guides assignment. It does not enter the outcome baseline. The arbitrary known effect τ is constant across clusters. Fit an intercept, treatment Z and the true spillover covariate using the validated SAR ML estimator.</p></div><div><b>Explicitly hypothetical sensitivity</b><p>Add βX with β = 1, X = (average rank − 0.5)/58, to the simulated baseline and adjust for X in the fit. This is a matched baseline-and-adjustment sensitivity. It does not assert that historical SUD incidence predicts education ability.</p></div></div><p>Control-only: S(Z) = γ(1 − Z)WZ. Both arms: S(Z) = γWZ. Queen weights are primary; rook weights are a focused geographic sensitivity. Noise has SD 1.</p><p class="note">This study evaluates estimation of a simulated effect. It does not measure an actual education effect or establish that education reduces SUD.</p></section>',
    '<section><h2>Observed inputs stay fixed</h2><img src="real_sud_assets/observed_years.png" alt="Maps of observed per-capita cluster SUD rates in 2018, 2019, 2020 and 2021"><p>Every yearly rate is deaths divided by the corresponding working-age population, multiplied by 100,000. No synthetic incidence surface is generated. Each year is a separate planning setting; pooled incidence is descriptive and supports the reference partition.</p></section>',
    '<section><h2>Adapting the grid designs</h2><img src="real_sud_assets/regions.png" alt="Map of four frozen geographic regions"><p>Incidence-guided saturation uses per-capita cluster rates → cluster ranks → mean rank within each frozen region → 80/60/40/20% saturation. Mean raw cluster rates and pooled regional deaths/population are separate summary sensitivities. Rules with identical regional ordering and ties reuse the same source explicitly.</p><img src="real_sud_assets/design_examples.png" alt="Example treatment assignments for all nine designs"><p>Graph Checkerboard is a greedy interspersion adaptation, not a literal rectangular checkerboard. Spatial blocks are frozen once per year. Balanced Halves allocates 14 and 15 treatments across its two 29-cluster halves. Design budgets otherwise retain their rule-specific differences, and population shares are reported separately.</p></section>',
    '<section><h2>Verification and Monte Carlo precision</h2><p>Stochastic designs start with 100 allocation draws and 100 outcomes per unique allocation. Proven singleton assignments start with 1,000 independent outcomes. Refinement preserves all prefixes and responds to the dominant uncertainty component. Targets are relative MC SE of mean MSE ≤ 5% and coverage MC SE ≤ 0.01.</p><p>Content manifests check inputs, frozen partitions, computation code, package versions, RNG and BLAS before checkpoint reuse. Boundary estimates, non-identification, aliases, warnings and failures are recorded explicitly. Duplicate assignments retain their draw frequencies without being counted as extra independent outcome fits.</p><p class="note">Allocation-tail summaries describe estimated conditional MSEs. Upper-tail membership can be uncertain with finite outcomes and draws. Split-half diagnostics and focused outcome refinement qualify tail comparisons. Reused yearly rows are shared evidence, not independent yearly confirmation. Error bars below concern each setting; they are not uncertainty for pooled averages across settings.</p></section>',
    evidence,
    '<section><h2>Explore current performance</h2><label>Evidence set <select id="dataset"><option value="main">Full study: main production estimates</option><option value="tail">Refined risk corners: same allocations, R ≥ 400</option></select></label><p>For allocation-downside comparisons select the refined risk corners. Full-study estimates cover the broader grid; the refined set adds outcomes at fixed allocations.</p><p id="resultStatus"></p><div id="filters"></div><label>Chart metric <select id="metric"><option value="Mean_MSE">Mean MSE (with MC uncertainty)</option><option value="Q90_Estimated">90th percentile of estimated allocation MSE</option><option value="Worst10_Mean_Estimated">Mean estimated MSE in worst 10% of draws</option><option value="Mean_Population_Share">Mean treated population share</option></select></label><svg id="chart" role="img" aria-label="Current performance chart"></svg><div class="scroll"><table><thead><tr><th>Design</th><th>MSE ± MC SE</th><th>Bias</th><th>Coverage ± MC SE</th><th>Treated</th><th>Population share</th><th>J / R</th><th>Complete / precision</th></tr></thead><tbody id="rows"></tbody></table></div></section>',
    '<section><h2>Development record and next deliverables</h2><ol><li>Study decisions and independent alignment review: complete; concrete implementation approved.</li><li>Input freeze, safeguards and behavioral verification: complete; smoke/pilot remain engineering checks.</li><li>Production and focused tail refinement: completed and separately verified; earlier outputs preserved.</li><li>Thesis chapter and appendix: integrate verified results, trace numbers and exhibits, render and review.</li><li>Standalone CTJ manuscript and supplement: derive from the chapter, verify current journal requirements and Project 1 citation.</li><li>Comprehensive Claude review prompt: prepare after the completed analyses and drafts.</li></ol><p><a href="../../docs/plans/nc_sud_implementation_plan_2026-10-02.md">Approved implementation plan</a> · <a href="../../docs/plans/nc_sud_independent_review_2026-10-02.md">Independent alignment review</a> · <a href="../results/real_sud_rev_20261002/">Current output directory</a></p><small>Only authorized 58-cluster aggregates are displayed. Restricted county source files are excluded. This document is an offline project artifact; it requires no external scripts or network requests.</small></section>',
    paste0('<script>const mainData=', rs_json(results), '; const tailData=', rs_json(tail_results), '; let data=mainData;</script>'),
    '<script>const fields=["Year","Model","Neighbor","Summary","Regime","Rho","Gamma"]; const fmt=x=>x==null?"—":Number(x).toFixed(3); const filters=document.getElementById("filters"); /* Rebuild controls for the chosen evidence set. */ function configure(){filters.replaceChildren();if(data.length){for(const f of fields){const label=document.createElement("label");label.textContent=f+" ";const sel=document.createElement("select");sel.id=f;for(const v of [...new Set(data.map(r=>String(r[f])))]){const o=document.createElement("option");o.value=v;o.textContent=v;sel.append(o);}label.append(sel);filters.append(label);sel.onchange=draw;}}draw();}document.getElementById("dataset").onchange=()=>{data=document.getElementById("dataset").value==="tail"?tailData:mainData;configure();};document.getElementById("metric").onchange=draw;
/* Render matched settings; intervals concern mean MSE Monte Carlo error only. */ function draw(){const rows=data.filter(r=>fields.every(f=>String(r[f])===document.getElementById(f).value)).sort((a,b)=>a.Design_ID-b.Design_ID);document.getElementById("resultStatus").textContent=data.length?`${rows.length} designs in this selected setting. Empty combinations were not part of the approved grid.`:"No performance file is loaded yet.";const tbody=document.getElementById("rows");tbody.replaceChildren();for(const r of rows){const tr=document.createElement("tr");const values=[r.Design,fmt(r.Mean_MSE)+" ± "+fmt(r.SE_Mean_MSE_Joint),fmt(r.Bias),fmt(r.Coverage)+" ± "+fmt(r.MCSE_Coverage),fmt(r.Mean_Treated),fmt(r.Mean_Population_Share),r.J+" / "+r.R,String(r.Complete)+" / "+String(r.Precision_OK)];for(const v of values){const td=document.createElement("td");td.textContent=v;tr.append(td);}tbody.append(tr);}const metric=document.getElementById("metric").value;const svg=document.getElementById("chart");svg.replaceChildren();const H=Math.max(200,rows.length*38+45),W=1060,L=345;svg.setAttribute("viewBox",`0 0 ${W} ${H}`);const max=Math.max(0.01,...rows.map(r=>(r[metric]||0)+(metric==="Mean_MSE"?1.96*(r.SE_Mean_MSE_Joint||0):0)));const scale=(W-L-45)/max;const el=(tag,attrs,text)=>{const e=document.createElementNS("http://www.w3.org/2000/svg",tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text)e.textContent=text;svg.append(e);};rows.forEach((r,i)=>{const y=20+i*38;el("text",{x:0,y:y+16,"font-size":14,fill:"#193247"},r.Design);if(r[metric]!=null){el("rect",{x:L,y,width:r[metric]*scale,height:24,rx:3,fill:r.Design_ID===9?"#59718c":"#167c80"});if(metric==="Mean_MSE"){const lo=Math.max(0,r.Mean_MSE-1.96*r.SE_Mean_MSE_Joint),hi=r.Mean_MSE+1.96*r.SE_Mean_MSE_Joint;el("line",{x1:L+lo*scale,x2:L+hi*scale,y1:y+12,y2:y+12,stroke:"#172d3c","stroke-width":2});}el("text",{x:L+r[metric]*scale+8,y:y+17,"font-size":13},fmt(r[metric]));}else el("text",{x:L,y:y+17,"font-size":13},"Incomplete: performance withheld");});}configure();</script></main></body></html>')
  path <- file.path(dest, "real_sud_companion.html")
  writeLines(html, path, useBytes = TRUE)
  cat("Saved companion: ", path, "\n", sep = "")
  invisible(path)
}
args <- commandArgs(TRUE)
rs_companion(profile = if (length(args)) args[1] else "pilot")
