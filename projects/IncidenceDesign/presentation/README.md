# presentation/

Talk deck for the IncidenceDesign project (dissertation Project 2): the 20-minute talk
at the SUDDEN project meeting on 2026-10-29, which also seeds prelim Part II
(prelim 2026-12-07). This deck also serves as an initial draft for building the slides in
`bios-dissertation/slides/prelim/part-2/`.

## Files

| File | What it is |
|---|---|
| `sudden_talk_2026-10-29.qmd` | Deck source (Quarto revealjs, 1920×1080). Speaker notes hold the read-aloud script (one sentence per bullet) with per-slide timing; backup slides have short "use when asked" notes. |
| `sudden_talk_2026-10-29.html` | Rendered deck, self-contained (fonts, figures and MathJax embedded). Open in a browser; **S** opens the speaker view (notes, next slide, timer), **M** the slide menu, **F** full screen. |
| `OUTLINE.md` | Slide-by-slide plan: sections, timing, figure, slide text and prelim carry-over tag for every slide, plus open questions. |
| `outline.html` | Visual version of `OUTLINE.md` with figures inline (built by `tools/build_outline_html.py`). |
| `theme/unc-gillings.css`, `theme/prelim-overrides.css`, `theme/footer-part.html` | Copies of the prelim template theme from `bios-dissertation/slides/prelim/test/`, unchanged. `footer-part.html` writes each slide's `data-footer-part` (the section name) into the footer. |
| `theme/deck.css` | Everything specific to this deck (cards, glossary, tables, schematic). |
| `figures/` | PNGs used by the deck and outline. Manuscript figures are rasterized copies of the current PDFs; `nc_design_examples.png` is deck-only. Older PNGs from Draft 1 of the outline remain for reference. |
| `tools/refresh_figures.sh` | Re-copies the current manuscript figures as PNGs and redraws the NC example map. |
| `tools/deck_allocation_figure.R` | Slide-sized allocation-risk figure (CTJ Figure 3 data, larger text, SRS in dark grey). |
| `tools/outline_from_deck.py` | Regenerates the slide entries of `OUTLINE.md` from the deck source. |
| `tools/nc_design_examples.R` | Draws one example allocation per design (six designs + SRS) on the 58-college map, 2018, title-free and named without numbers. |
| `tools/build_outline_html.py` | Builds `outline.html` from `OUTLINE.md`. |

Reviews of the deck: `docs/reviews/2026-10-28-sudden-talk-*`.

## Build

From the project root:

```bash
bash presentation/tools/refresh_figures.sh
```

```bash
cd presentation && /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render sudden_talk_2026-10-29.qmd
```

```bash
python3 presentation/tools/outline_from_deck.py && python3 presentation/tools/build_outline_html.py
```

Use RStudio's bundled Quarto; `/usr/local/bin/quarto` fails on this machine. The
rendered HTML is about 8 MB because everything is embedded, so it can be presented
offline.

## Conventions

- Section names (Background, Methods, Simulation Results, Application, Discussion,
  Backup) appear in the footer; there are no section-divider slides.
- Slide titles are claims stated as sentences; figures carry no in-figure titles.
- Design names exactly as in `AGENTS.md`; no DIM; numbers trace to the V5 manuscript,
  Chapter 3 or `results/`.

## Adapting this deck for the prelim (Part II)

This deck is a standalone talk. The prelim deck (`bios-dissertation/slides/prelim/`, skeleton
in `test/slide-deck-test/prelim-test.qmd`) has one title, one agenda, Parts I–III, then one
closing slide and three backup sections. Copy this deck into
`bios-dissertation/slides/prelim/part-2/` only after Andrew approves it, then make these
changes (author, 2026-10-28):

1. **Cut the standalone frame.**
   - Drop the title slide and the Outline slide. The prelim has its own title and agenda.
   - Drop the closing thank-you slide. Its acknowledgements move to the single prelim
     "Thank you" slide after Part III, grouped with Parts I and III: Dr. Lin, Ashkan
     Habib, Dr. Simpson and the SUDDEN team, plus the IRB protocol 14-2036 note.
   - "Manuscript in preparation for *Clinical Trials*" goes on the last Part II
     Discussion slide.
2. **Open with the Part II divider.** Use the skeleton's
   `.section-divider .divider-carolina` slide with
   `data-footer-part="Part II · SpatialCRT: Incidence"` and `data-menu-title="Part II: …"`.
3. **Footer labels.** Prefix every section label with the part, so the footer tracks both:
   - `data-footer-part="Methods"` becomes `"Part II · Methods"`, and the same for every
     section;
   - give the menu titles a `II:` prefix, as the skeleton does.
4. **Keep, rewrite or drop each slide** by its carry-over tag in `OUTLINE.md`:
   - **prelim:** keep as is;
   - **technical:** rewrite for biostatisticians, with more model and estimation detail;
   - **talk-only:** drop. This covers the SUDDEN-specific framing and the questions for
     the group.
5. **Build on Part I, don't repeat it.** Part I comes first and introduces cluster trials,
   spillover, spatial dependence and the spatial lag model. Keep the Background content
   most relevant to Project 2 (heterogeneous incidence, SUD, the pilot) and frame it as
   the next step from Part I. Keep Part II's own model slide in full. Its notation (τ, γ on
   the share of treated queen neighbors, β for incidence, at the cluster level)
   deliberately differs from Part I's, and each part defines its own.
6. **Backup.** Move the backup slides into the second of the three backup sections, after
   Part III, under the skeleton's `Backup · Part II` divider. Prefix their titles with
   "Backup ·" and give them the `.backup-slide` class.
7. **Theme and figures.**
   - Keep the theme consistent across all three parts.
   - Merge `theme/deck.css` into the prelim theme, checking its class names against
     Part I's `deck.css` for clashes and namespacing them if needed.
   - Copy `figures/` into `part-2/figures/` and update the paths.
8. **Timing.** Aim for about 20 minutes per part, about 60 in total. Parts I and II may
   run a little longer while Part III is still brief. Re-time the speaker notes after the
   cuts and rewrites above.
