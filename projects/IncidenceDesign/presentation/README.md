# presentation/

Talk deck for the IncidenceDesign project (dissertation Project 2): the 20-minute talk
at the SUDDEN project meeting on 2026-10-29, which also seeds prelim Part II
(prelim 2026-12-07). The deck stays here until Andrew approves copying it to
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
