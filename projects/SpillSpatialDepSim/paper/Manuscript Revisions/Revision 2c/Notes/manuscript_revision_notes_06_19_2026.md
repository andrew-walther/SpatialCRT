# Manuscript Revision Notes — 6/19/2026 Zoom with Feng-Chang

---

## Pre-Meeting Notes

- Review the markup & notes included in "Revisions_V2b_markup"

---

## Key Issues to Address

### W Matrix Notation (Weight Matrix)
- W = sum_{ni}
- Y = cluster-level mean

### Cluster-Level Outcome (Critical)
- **Big key:** specify that y_i is a cluster-level outcome (not individual-level — we model cluster mean because individual outcomes in a cluster are identical)

### ICC (Intra-Cluster Correlation)
- Model for simulation doesn't consider ICC for estimation — add a comment in the paper
- Include in future work / limitations

---

## Equation Edits

- **Bigger picture:** can we drop the "general" section with equations 1 & 2 and move on to focusing on the specific implementation that starts with equation 3?

### EQN 1 (can be dropped — general form)
- N(i) isn't referenced anywhere (need to specify/define)
- "where j in N(i)" in bottom of EQN summation (instead of j !=i)

### EQN 2 (can be dropped — just a reiteration in matrix form)
- Move under EQN 3 with Psi*Z term?

### EQN 3 (the "best" model for this application)
- Drop k & k' (don't mention specific individuals — we specify at cluster-level)
- Drop section with equations 1 & 2, but move a vector/matrix form of EQN 2 with representation of EQN 3 (with Z*Psi term)

### EQN 5
- Drop k' (individual index)

### EQN 6
- Drop k's for individuals (model cluster level)

---

## Response Letter

- State intent to model cluster-level mean so that ICC (intra-cluster correlation) can't be explored
- Flag as current limitation

---

## Appendix

- Appendix A is good (cluster level)
- Appendic B is good (cluster level)
- Drop equation numbers from appendix sections if not referenced in main text
  - Add `*` to LaTeX equation calls

---

## Model Notation Convention

- **Ordering of equations:** write with alpha + linear predictor + spatial dependence (prefer spatial dependence at the end) — not a big deal though

---

## Nits / General Formatting

- Add comma "," after formulas
- Don't indent "where" after an equation
- Indent for a new subsection (after section heading)
- Check for indent/no indent to ensure consistency
- Drop em-dashes "—" (I think we already addressed this)

---

## Zoom Meeting Summary (6/19/26)

### Quick Recap
Feng-Chang and Andrew discussed notation improvements for the model in the paper. They focused on simplifying the mathematical notation by consolidating redundant equations and clarifying the use of cluster-level outcomes (Y_I) rather than individual-level outcomes (Y_IJ), which would help address reviewer concerns about intra-cluster correlation (ICC). Feng-Chang suggested reorganizing the equations so that Model 3 becomes Model 1, and adding spatial dependence terms in a more consistent manner. They also discussed formatting conventions for mathematical expressions, including proper indentation and use of commas after formula breaks. Andrew agreed to revise the notation based on their discussion and update both the paper and the response letter to reviewers.

---

## Next Steps

### Andrew
- Revise the model notation in the manuscript by making Equation 3 the new Equation 1, remove redundant equations (especially Equation 2), and update the notation to consistently use Y_I (cluster-level summary) instead of Y_IJ (individual-level) throughout
- Ensure all references to K and K' (e.g., K prime) are removed from the equations and related text
- Update the manuscript to explicitly state that the current model does not account for intra-cluster correlation (ICC) and discuss this as a limitation, suggesting it as a potential direction for future research
- Standardize the order of terms in equations (e.g., alpha, beta, Xi) for consistency across scalar and vector forms, with the spatial dependence term at the end in scalar forms as preferred by Feng-Chang
- Check and correct LaTeX formatting for equation breaks, indents, and punctuation (e.g., adding commas after equations that are part of a sentence, ensuring consistent paragraph breaks)
- Remove or consolidate redundant equations (e.g., Equation 7) and consider dropping equation numbers for equations not referenced in the main text
- Review and update the response letter (Word document) to reflect changes in notation and ICC discussion, ensuring consistency with the revised manuscript
- Send the updated version of the manuscript to Feng-Chang for further review

### Feng-Chang
- Review the updated manuscript and response letter to ensure all reviewer comments are fully addressed, especially regarding notation and ICC discussion

---

## Discussion Notes

### Cluster-Level Data Modeling Approaches
Feng-Chang recommended using Y_I to summarize individual outcomes into clusters, as this approach better handles spatial dependency and spillover effects without introducing intra-cluster correlation (ICC) complications. If ICC is discussed elsewhere in the paper, Andrew should clearly state that the current model does not consider ICC and could be a topic for future research.

### Equation Simplification
Feng-Chang suggested merging equations 3 and 1 and removing redundancy in equations 1 and 2. They also discussed modifying the formulation to include N(i) and adjusting the summation notation. Key focus was on ensuring the equations accurately reflect the model descriptions and are logically consistent.

### Mathematical Notation Revision
They agreed to move specific form content to the top of the section and remove redundant general form content. Feng-Chang suggested adding a vector form representation and recommended removing K and K prime from certain expressions.

### Model Updates and Limitations
Andrew confirmed he had updated the model to include comments about potential outcomes and specifications, while also addressing the interpretation of spillover effects. Equation 7 was identified as redundant with Equation 4.

### Appendix Formula Consistency
They discussed consistency in the ordering of alpha, beta, and spatial dependence terms in the formulas. Feng-Chang suggested putting the alpha term first and the spatial dependence term at the end. The vector form is to be kept as-is; individual scalar form adjusted for consistency.

### Paper Revisions and Formatting
Feng-Chang advised that unnecessary references could be dropped during revisions and provided guidance on proper formatting of mathematical formulas, including the use of commas and periods after formulas. Andrew will send an updated version of the paper over the weekend for further review.
