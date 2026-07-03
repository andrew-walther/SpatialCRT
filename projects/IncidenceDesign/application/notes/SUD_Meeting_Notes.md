# SUD Intervention Study — Meeting Notes

## 2026-04-27

### Notes
- Moving into the intervention phase.
- Exploring a county-level pilot in **Lenoir County**.
- The intervention will **not** be census-tract based — it operates at the county level, since counties are irregularly shaped and the intervention is applied county-wide.
- Likely to work with county-level incidence data.
  - All counties appear to have a Community College (CC) where the intervention could be applied.
  - Present coverage as a **network**: one county hosts the primary school, surrounding counties have satellite campuses (e.g., Durham Tech is also present in Orange County).
  - **Spillover** occurs because schools talk to each other.
  - We know which CC serves which county.

### Questions from Dr. Simpson
- Other dataset — needs to be rebuilt.
- Census tract data.

### Potential Study Intervention
- Data shows SUD is a **social problem**: driven by poor access to care and social determinants (access/money).
- Working with Lenoir County Community College ([directions](https://maps.app.goo.gl/bJe9Ma3VMNw42AkEA)).
- Working with workforce development, teaching disaster response.
  - Outreach beyond just transporting people to the hospital.
- Pilot could later spread to other community colleges.
- No intervention/evaluation design finalized yet.
- Goal: empower CC students/staff to get people to the hospital and share treatments.
- Potential spillover to other counties.
- Community colleges have satellite sites in other counties:
  - Kingston, NC
  - Lenoir County, NC
  - Jones County
- Hope is that other schools adopt the initiative.
- Developed an algorithm for estimating incidence rate via death certificates.
  - Not necessarily aiming for precise incidence — reasonably accurate estimates are sufficient (see Gan paper reference).
- Plan: pilot study first, then share findings with all other counties in the state.
- **Measuring effect**
  - Won't see SUD rates drop immediately.
  - Can observe process changes: medication adherence, timely clinic visits when needed.
- Lenoir County floods regularly — SUD rates spike after flooding events.
- The SUD program is structured so students earn credit for participation/work (get promoted) — worth sharing with the CC system.
- Crude SUD incidence rates vary substantially across census tracts statewide.
  - Varies by median income level.
  - Varies by urban vs. rural setting.
- CCs actively recruit from underserved/high-risk areas — targeting groups where intervention could do the most good.
- A census-tract-level intervention isn't feasible for the teaching intervention — it must be applied at the county level.
- **Open question:** How do we get SUD incidence data for all counties in NC?
  - Needed for the incidence-driven component of the design.
  - Consider applying this now — potential role for an epi PhD student.
  - We already know which CC (and its satellites) serve each county.

### Ideas
- Treat CCs as a spatial block/unit.
- Start with Lenoir CC — pilot a "teach" lecture there.
- NC has 58 community colleges spanning 100 counties.
  - [NC Community Colleges — list of colleges](https://www.nccommunitycolleges.edu/students/what-we-offer/colleges/)
  - [Wikipedia: NC Community College System](https://en.wikipedia.org/wiki/North_Carolina_Community_College_System#Colleges) — shows counties served by each school.
- Possible current EPI student collaborator: Ashkan Habib?

### Broader Notes
- Made real methodological progress on controlling for spillover/contamination.
- Want to make the method broadly usable by others post-publication.
- Consider submitting to a clinical trials journal — cluster designs are getting significant coverage currently.

---

## 2026-05-01

### Updates
- Worked on the applications section of the report.
- Report currently uses **synthetic** incidence data (designed to be plug-and-play with real data later).

### Meeting Notes
- NC death certificate data comes as a text file.
- Need to classify death certificates as SUD or not (at county or census-tract level).
- Know SUD counts by county (cases and population).
  - Aggregate counties into clusters (groups of counties covered by a single CC).
  - This yields total cases and population for a SUD rate per cluster.
- **Task breakdown**
  - Clean death certificate data → determine SUD case status.
  - Identify county of residence.
  - Group by: SUD cases & county → clusters.
- **Metrics of interest**
  - Sudden deaths / population (per 100k)
  - Sudden deaths / all deaths
- Consider reviewing historical high-incidence regions — potentially merge multiple CC coverage regions into a single cluster if the combined area is high-incidence (e.g., CC-1 & CC-2 → one larger cluster for intervention).

### Next Up
- Obtain real incidence data (confirm SUD formula).
- Apply real incidence data to the application study.
- Continue manuscript writing.
- Read Ashkan's MS paper.
  - Figure 2 is key (SUD/population rate).
  - Offer co-authorship to Ashkan for contributing SUD data.
- Dr. Lin to send a follow-up to Ashkan requesting tables for Figures 2 & 3 (for Andrew).
  - Schedule a meeting after 5/11 (Andrew to find a common time).
  - Possible slot: during the regular 2pm meeting on Friday, 5/8.

---

## 2026-05-08

### Ashkan's Work
- Collected death certificates from 2018–2021 (~100k deaths).
- Used an algorithm to classify deaths as **Sudden (SUDDEN)** or not.
  - This is an approximation of sudden death outcomes.
  - A physician panel reviewed cases to classify them — essentially a decision tree for inclusion/exclusion.
  - Reference paper: [Springer, J Gen Intern Med 2018](https://link.springer.com/article/10.1007/s11606-018-4771-5)
- Used death certificate zip codes to pool cases into counties.
  - Restricted to ages 18–64.
  - Population-at-risk estimated using SEER data (NCI cancer registry) for matching years.
    - Provides estimated county-level population by year.
- Calculated incidence rate from the above.

### Dataset Setup
Andrew to add:
- Community college coverage cluster tag
- Year
- County
- Population at risk (in county)
- Sudden deaths (primary outcome)
- ~~Deaths~~ (removed)
- Cardiovascular disease-related deaths (alternate outcome)
- Other covariates:
  - Income
  - Food deserts
  - Physicians (access)
- GIS info: average distance to nearest hospital (miles)

### Implementation Considerations
- Consider testing on 4 different years to demonstrate robustness, or make year-specific recommendations.
- Ashkan found SUD rate increased year-over-year across all counties.
- Alternative: average population & deaths across the 4 available years.

### IRB / Data Access
- Need access to the **"J Drive"** by joining Dr. Simpson's IRB.
  - Need to determine the process for this.

### Collaboration Plan
- Work together with Ashkan & Dr. Simpson.
- Acknowledge the SUDDEN IRB (go through the application process).
- Pass along the aggregated dataset.
- Send a note to Dr. Simpson to kick off the process.
  - Andrew to send the first draft/attempt.
