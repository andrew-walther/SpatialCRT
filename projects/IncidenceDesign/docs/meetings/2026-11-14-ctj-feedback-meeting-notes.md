# CTJ manuscript feedback meeting: decisions and revision list

Written 2026-11-14. Author: Andrew Walther.

Sources: Andrew's notes from the manuscript review meeting with Feng-Chang Lin on 2026-11-10,
and the Zoom AI summary of the same meeting. Both are preserved at the end of this file.
Feng-Chang's written annotations (the FCL register in
`paper/ctj_manuscript/revisions/2026-11-10-paper-2-revisions/`) followed the meeting.

**Priority.** Chapter 3 ready for the prelim is the deliverable. The targeted revisions are made
first in the CTJ manuscript and then carried into Chapter 3. The CTJ 3,500-word limit is not a
priority now.

**Reading the Zoom summary.** It contains transcription errors. "Sunderland data" is the SUDDEN
data. "gamma XZ" should read γ S(Z), with γ a scalar multiplier outside S(Z). "trickle ball"
is probably "Checkerboard", and "trackable" probably "tractable". Where the summary and Andrew's
notes disagree, Andrew's notes govern.

## Decisions from the meeting

| ID | Decision | Applies to | Status |
|---|---|---|---|
| D-1 | Reparameterize the outcome model. β multiplies the treatment indicator Z and is the intervention effect, τ multiplies the incidence covariate X, and spillover is written γ S(Z) with γ a scalar outside S(Z) = WZ. | Chapter 3, CTJ, SI, exhibits, code comments; check Chapter 2 notation | Already recorded in DECISIONS 2026-11-11, not yet applied |
| D-2 | The coefficient on incidence is pre-specified in the simulations, not estimated. Its value is given in the simulation settings, not in the model definition. Investigators choose it for their own outcome, by whether they believe incidence is associated with that outcome. | Chapter 3, CTJ | To do |
| D-3 | Rename "Application" to "Case Study". The North Carolina section uses observed incidence and real county or college geography as fixed inputs, with simulated outcomes, so it is a simulation-based case study, not an analysis of trial data. Working title: "Case Study: Simulation Based on Sudden Unexpected Death in North Carolina". Its description can stay in the Methods. | Chapter 3, CTJ, SI | To do |
| D-4 | Report confidence intervals with a citation for the t interval with the degrees-of-freedom correction and for the normal (Wald) interval. A general quantile notation can cover both. Both stay reported, Wald in the supplement. | Chapter 3, CTJ, SI | To do (citations needed) |
| D-5 | Estimating the incidence coefficient after a real trial, to test whether incidence is associated with the outcome, is future work. | Chapter 3 Future Considerations | To do |

## Revision list

Each item has an ID for the triage table in the Chapter 3 review session. "Where" names the
target documents. Cross-reference each item with the written FCL annotations, since many will
overlap.

### Notation and model (Methods)

| ID | Revision | Where |
|---|---|---|
| M-01 | Apply D-1 throughout (β on Z, τ on X, γ S(Z)). Update every derived statement, for example the Checkerboard estimand: Checkerboard estimates the intervention effect minus the spillover effect, τ − γ in the current notation and β − γ after D-1, not a contrast. | All |
| M-02 | Write the model equation on one line. Drop ε ~ N(0, I) from the equation and state it in the prose: ε is independent standard normal noise, the intervention effect is in standard-deviation units, and ρ is the spatial dependence of the outcomes. No parentheses. | Ch3, CTJ |
| M-03 | Define the spatial weight matrix W clearly. Do not write "the ith element of WZ is the proportion of cluster i's ...". Define the indexing first or describe the quantity in words. | Ch3, CTJ |
| M-04 | "X is cluster incidence" becomes "X is cluster-level incidence". | Ch3, CTJ |
| M-05 | Replace "baseline outcome" for the incidence term. The study is not longitudinal, so say the incidence term shifts the outcome (an intercept shift). Rewrite "the coefficient sets how strongly incidence shifts the baseline outcome" accordingly. | Ch3, CTJ |
| M-06 | Move the coefficient's value (currently "we set it to 1 in the grid study") to the simulation settings, and explain why it is pre-specified (D-2). | Ch3, CTJ |
| M-07 | Move "queen contiguity for the primary analysis, rook contiguity as a sensitivity analysis" from the model definition to the simulation design subsection. | Ch3, CTJ |
| M-08 | "Checkerboard is identified" becomes "identifiable". Under rook contiguity the effects are not identifiable, and under queen contiguity they are identifiable but imprecise. Drop "(cannot be estimated)", which the preceding words already say. | Ch3, CTJ, SI |
| M-09 | Delete the sentence saying a design with only one possible assignment has no spread. It is obvious. | Ch3, CTJ |
| M-10 | Cite the t interval, the normal interval, and the standard-error correction SE·√(N/(N−p)). The sources need not be spatial papers (D-4). | Ch3, CTJ, SI |
| M-11 | Define m(a), the conditional MSE of allocation a that underlies the downside allocation risk, clearly at first use, in words as well as the formula. | Ch3, CTJ |
| M-12 | Define MSE at first use. | Ch3, CTJ |
| M-13 | Split "We estimated τ ... incidence" into separate sentences, one idea each. | Ch3, CTJ |

### Simulation design

| ID | Revision | Where |
|---|---|---|
| S-01 | State that the SUDDEN data are used only in the case study, not in the grid simulation, which uses a 10 × 10 grid. | Ch3, CTJ |
| S-02 | Write out the full generating model for incidence X, including where ρ_X comes from. | Ch3, CTJ (full detail in SI) |
| S-03 | Justify the five incidence configurations or reduce them. They did not change which design performed best. Say that, and consider moving the detail to the appendix or SI. | Ch3, CTJ |
| S-04 | Explain why incidence is on the 0 to 1 scale, or relax it. It does not have to be. | Ch3, CTJ |

### Case study (currently "Application")

| ID | Revision | Where |
|---|---|---|
| A-01 | Rename the section per D-3 and update every reference to "application". | All |
| A-02 | Say plainly what is observed (incidence, geography) and what is simulated (outcome Y). Write "The observed rates were used as fixed inputs" and drop "rather than simulated". | Ch3, CTJ |
| A-03 | Explain how blocks were built on irregular geography: k-means clustering on the pairwise distances between clusters, to form compact groups that stand in for the 2 × 2 blocks. Check the code first. Chapter 3 does not mention k-means now, and the only mention found is in `code/14_manuscript_supplement_figures.R`. | Ch3, CTJ, SI |

### Introduction

| ID | Revision | Where |
|---|---|---|
| I-01 | Cite sources (two or more) for "Two allocations with the same number of treated clusters can provide very different information about the treatment effect". Add references to the first paragraph generally. | Ch3, CTJ |
| I-02 | Define Checkerboard when the prior work is first discussed, before any other use, and state there that it is the prior study's block stratified sampling design. Then use "Checkerboard" only, capitalized consistently. | Ch3 (largely done 2026-11-11), CTJ, SI |
| I-03 | Rewrite "in a strict checkerboard with edge neighbors, a cluster's exposure to treated neighbors is fully determined by its own treatment", which is confusing. | Ch3, CTJ |
| I-04 | Replace "that work" with "the prior project" or name the citation. In CTJ, keep the self-contained wording ("In previous work ...") per the standing rule. | Ch3, CTJ |
| I-05 | Revise the sentence "The comparison includes ... here called Checkerboard". Move the Checkerboard definition earlier (I-02), drop "meaning complete randomization of a fixed number of clusters", make clear that the comparison is this paper's, and do not single out SRS and Checkerboard without saying why. The comparison is not limited to two designs. | Ch3, CTJ |
| I-06 | "Our application supports planning" becomes "In this paper, we plan ..." (with "case study" per D-3). | Ch3, CTJ |

### Writing throughout

| ID | Revision | Where |
|---|---|---|
| W-01 | Avoid parentheses unless they are necessary, and keep the content in the main sentence. | All |
| W-02 | Break up long sentences with many commas. | All |
| W-03 | Use formal transitions. Avoid "so" as a connective. | All |
| W-04 | Define every term at first use, and avoid undefined terms such as "tractable". | All |

### Open questions (not decided in the meeting)

| ID | Question | Suggested handling |
|---|---|---|
| Q-01 | Spillover is proportional to the share of treated neighbors here and binary in Project 1. Should binary-spillover results be included? | **Settled (Andrew, 2026-11-14):** include binary-spillover results only where they already exist, as supplementary results, and run no new analysis. They exist for the case study: the Chapter 3 appendix compares binary, count and proportional spillover forms on the 58 North Carolina colleges (`Thesis_Chapter.qmd` around lines 1656–1662, in "Participation, Power and Sensitivity Analyses for the NC Application"). The grid study has only the proportional form. State the Project 1 / Project 2 difference in Methods and point to that appendix result. |
| Q-02 | What does the incidence coefficient's value represent, and how should it be set? | In practice it was specified as 0 or 1: 1 in the grid study, 0 for the primary case-study outcome, and 1 as a hypothetical sensitivity analysis. Describe it that way, as a choice of whether incidence shifts the outcome, per D-2. Work out and state what a value of 1 means on the scale used, with incidence on [0, 1] and unit noise variance, and check this against the code. Alternative setups raised in the meeting (an on/off indicator, a value chosen to suit the investigator's outcome) go in the Discussion or Future Considerations with D-5. |

## Original notes

### Andrew's notes

Project 2 Manuscript Feedback

* Using "m(a)" for downside allocation risk - make sure this is clearly explained/defined
* Current spillover: proportional to neighbors that are treated/not (update to binary spillover as well?)
   * Project 1 (binary) vs. Project 2 (proportional)
* Introduction
   * "Two allocations with the same number of treated clusters can provide very different information about the treatment effect" ^{3,4} (cite something)
   * Make sure to be consistent with "Checkerboard" vs "checkerboard"
      * Introduce the design of checkboard before referencing (reference it without defining)
   * Mention
   * Clean this up "in a strict checkerboard (define this) with edge neighbors, for example, a cluster's exposure to treated neighbors is fully determined by its own treatment" (this is confusing)
   * Clean up/unify usage of "block stratified sampling" vs "checkerboard
   * Instead of "that work" say "the prior project"
   * "The comparison includes… here called Checkerboard" (make the definition earlier when discussing prior work)
      * Drop "meaning complete randomization of a fixed number of clusters"
      * Make clear that the comparison is for this paper/project
      * Don't specifically call out SRS & Checkerboard w/o making clear
   * "Our application supports planning" --> "In this paper, we plan"
* Method

Model

* "Queen contiguity was used for the primary analysis and rook contiguity for a sensitivity analysis" --> move to the actual simulation subsection (not in model definition)
* Make "ith element" of WZ more clear??? (clearly define the weight matrix?)'
   * Avoid using "ith" in "the ith element of WZ is the proposition of cluster I's
* For spillover (S(Z)): write "gamma" outside of S(Z), where gamma is spillover magnitude (use as a scalar outside the neighbor/treatment product)
* Move error term "e ~ N(0,I)" into the prose (drop from the equation)
* "X is cluster incidence" --> "X is cluster-level incidence"
* Don't use parentheses unless absolutely necessary. Should be able to keep them in a main sentences
   * "e is independent standard normal noise, where tau is in units of standard deviation and rho is the spatial dependence of the outcomes (no parentheses)
* "the coefficient Beta sets how strongly incidence shifts the baseline outcome"
   * Is "baseline outcome" the right term? --> intercept or outcome shift
* Don't specify "We set B=1 in grid study" until the "Simulation study" subsection
   * Specify as an indicator??? Instead of 0/1
   * Specify "B" based on the outcome they desire?
   * In simulation, give Beta a pre-specified value
   * "Investigators can choose B to suit their own outcome" depending on whether they believe incidence influences the outcome? (chosen based on whether incidence is believed to be associated with the outcome)
* Estimate checkerboard with "tau - gamma" (no contrast)
* "We estimated tau,…incidence (separate thoughts?)
* Add citation to t-distribution confidence intervals (also cite the normal Z-score coverage intervals)
   * Not from a spatial paper?
   * What is "SE_hat sqrt(N/(N-p))"??? Cite this!
* "checkerboard is identified" --> "identifiable" (but not with rook)
   * Effects can be identified by checkerboard but are imprecise? (use identifiable instead of identified)
   * Drop "(cannot be estimated)" (implied by words before it)
* Drop "a design with only one assignment has no spread…" (obvious)
* Outcome Y, if depending on "X" --> can we estimate "B"?
   * Not intended to be estimable, instead is specified, can we estimate it?
   * Would tell us if there is association of the incidence & outcome when we have real data??? (if incidence & outcome are known, estimating Beta after a trial is run tells us if there is an association of incidence~outcome)
      * B=1 gives magnitude of association currently? (need to think about this)
         * FUTURE WORK

Designs

* No feedback currently

Simulation

* SUDDEN data is not used in simulation, only in application
* Simulation was regularized with 10x10 grid
* Where does "rho_X" come from?
   * Write down the full model for "X"
* 5 different incidence generation methods? Are they necessary?
   * Why is the domain from 0 to 1? Doesn't have to be?
   * Parameter of interest is beta more commonly
   * Should we change to "BZ" and then "TauX" (more common convention)
      * Strange because "tau" is the parameter of interest here!

Application

* "The observed rates were used as fixed inputs" (drop "rather than simulated")
* Not called "Application"? - Call it a simulation study still?
* Is outcome "Y" simulated here?
   * Inputs are not simulated
* Section name "Case study: Simulation with sudden unexpected death in North Carolina"
   * Or "Simulation 1" and "Simulation 2"
   * Decision: call it a "case study" instead of "application"
   * Can still stay in the methods
* Make blocks of counties/clusters more clear (K-means) --> check on this
   * Build compact clusters to build the "2x2 blocks"
   * Computed pairwise distances between clusters? --> grouped clusters

* Results
   * N/A for now
* Discussion
   * N/A for now

### Zoom AI summary (as received; contains transcription errors noted above)

**Quick recap.** The meeting focused on Feng-Chang reviewing Andrew's manuscript and providing feedback on language, notation, and structure. Feng-Chang suggested several editorial changes to improve clarity and consistency, such as standardizing the use of "checkerboard" and avoiding undefined terms like "trackable." They discussed the model's notation, with Feng-Chang recommending that the spillover effect be written as gamma times SZ rather than including gamma inside the SZ function, and suggested using beta for the treatment effect and tau for the incidence effect to align with common statistical conventions. Feng-Chang also questioned the use of "baseline outcome" and suggested more precise terminology, as the study is not longitudinal. They clarified that beta is a pre-specified parameter in the simulation rather than an estimated one. A significant portion of the discussion centered on the structure of the paper, where Feng-Chang advised renaming the "application" section to a "case study simulation" to better reflect that it uses real spatial data but simulated outcomes. They also discussed the use of k-means clustering on pairwise distances between clusters to adapt the checkerboard design to the irregular geography of the sudden death data. Andrew explained that the simulation results showed different incidence generation methods did not meaningfully change which design performed best.

**Next steps (Andrew).**

* Revise the introduction to cite more references in the first paragraph to strengthen the argument.
* Ensure consistent capitalization and usage of the term "checkerboard" throughout the manuscript.
* Revise the sentence "If I have a strict checkerboard with edge neighbors, a cluster's exposure to treated neighbors is fully determined by its own treatment" to be clearer and more precise.
* Be consistent in terminology regarding "block stratified sampling" and "simple checkerboard," either defining them as the same or using one term throughout.
* Clarify the reference to "that work" in the introduction by explicitly referring to the previous project or citation number 9.
* Consider deleting the paragraph explaining the meaning of complete randomization of a fixed number of clusters and revise the first sentence of the following paragraph to cover the motivation and intention of the current paper.
* Clarify the "additional strategies" section to explicitly state that the comparison is not limited to only two designs.
* Move the discussion of Rook vs. Queen continuity from the model definition section to a later subsection where it is more relevant.
* Rephrase the use of "i-th element" to avoid confusion since the indexing system is not clearly defined.
* Change the notation in the model from "gamma SZ" to "gamma XZ" (or similar) to keep gamma as a scalar coefficient outside the neighbor/treatment product.
* Combine the model equation into a single line and move the epsilon definition into the sentence or the approach section.
* Remove the word "so" from academic writing and use more formal transitions.
* Replace the term "baseline outcome" with "intercept shift" or "effect on the outcome" to avoid confusion with longitudinal study terminology.
* Move the discussion of beta values (e.g., beta = 1) from the model section to the simulation study/data settings section.
* Add a more specific explanation of why beta should be chosen/pre-specified rather than estimated in the simulation context.
* Break up long sentences with too many commas into separate, clearer thoughts.
* Find and add a citation for the t-distribution correction used for confidence intervals.
* Define "MSE" upon its first mention in the text.
* Revise the sentence about "a design with only one possible assignment has no spread" or remove it as it is too obvious.
* Rename the "application" section to "case study simulation" or similar to distinguish it from the actual pilot study application.
* Write down the spatial model for X (incidences) explicitly in the simulation section.
* Consider switching notation from "tau Z" to "beta Z" for the intervention effect and "tau X" for the nuisance parameter to align with common statistical conventions.
* Clarify the description of the k-means clustering approach, specifying that it is applied to pairwise distances between clusters to create compact groups.

**Next steps (Feng-Chang).** Send the document with comments and suggestions to Andrew.

**Summary.**

*Paper draft review.* Feng-Chang and Andrew reviewed a paper draft, focusing on language suggestions and clarity improvements. They discussed issues with terminology consistency, particularly around "trickle ball" design and the need to define terms clearly. Feng-Chang suggested reorganizing some content, moving certain explanations to later sections, and recommended avoiding complex sentences with multiple commands. They also discussed the proper placement of beta coefficient discussion and agreed to make the writing more formal and less conversational.

*Manuscript revision and statistical methods.* Andrew and Feng-Chang discussed manuscript revisions, with Andrew noting the 3,500 word limit and requesting feedback on areas that could be moved to supplementary material. They discussed statistical methods, particularly the use of t-distribution versus normal distribution for confidence intervals, with Feng-Chang suggesting a more general notation approach using quantiles. Andrew agreed to look into proper citations for the statistical methods used and to consider reporting both Z and t-distribution results in supplementary material.

*Methodology revision.* Feng-Chang and Andrew discussed clarifications and revisions to a paper's methodology section, focusing on treatment assignment designs and spillover effects. They agreed to remove a redundant sentence about designs with only one possible assignment having no spillover. The conversation outlined the paper's structure, which includes sections on the model, different treatment designs, a simulation study, and an application using Sunderland data, though the specific application method remains unclear.

*Section structure.* Andrew and Feng-Chang discussed the structure and terminology of a research section involving simulation and application of design methods. They clarified that both the simulation section using a 10x10 grid and the section applying designs to real sudden death data in North Carolina should be referred to as simulation studies, with the latter specifically using a spatial cluster structure rather than generated data. Feng-Chang suggested calling the second section "simulation study based on spatial clusters in sudden death" to avoid confusion with actual applications to pilot studies.

*Terminology.* Andrew and Feng-Chang discussed terminology for two different simulation studies, agreeing to use "case studies" instead of "simulation studies" to better distinguish them and avoid confusion with existing data analysis. They reviewed model specifications, with Feng-Chang suggesting changes to make beta Z the parameter of interest rather than tau X, and questioned the necessity of limiting incidence values to 0-1. The discussion also covered the k-means clustering method used to group clusters based on geographic distance, with Andrew agreeing to clarify this approach in the documentation. Feng-Chang committed to sending feedback on the document, which Andrew plans to incorporate into the proposal due that weekend.
