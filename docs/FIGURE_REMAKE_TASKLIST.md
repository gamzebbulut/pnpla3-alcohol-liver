# Figure remake list — the approved contrast notation

Written 8 September 2026. Companion to `PLAN_notation_and_repo.md` §1, which is
now unblocked: the glyph is settled.

**Approved, and implemented:**

| internal key (unchanged) | what the reader now sees |
|---|---|
| Alcohol on wild type | Alcohol alone |
| Alcohol on variant   | Alcohol (PNPLA3<sup>I148M</sup> bg) |
| Variant alone        | PNPLA3<sup>I148M</sup> alone |
| Variant on alcohol   | PNPLA3<sup>I148M</sup> (alcohol bg) |

I148M is a true superscript, drawn with plotmath. Every panel below inherits it
from `nc_style.R`; no panel script hardcodes the wording any more.

---

## 1 · What changed in the code

**`nc_style.R` — three edits, and everything else follows from them.**

1. **A display layer.** `CONTRAST_DISPLAY` (flat), `CONTRAST_PLOTMATH`
   (superscript, one line), `CONTRAST_PLOTMATH_2L` (superscript, two lines), and
   three helpers: `nc_label()`, `nc_label_expr()`, `nc_labeller()`. The internal
   keys and the `contrast` / `cmp` columns of the nineteen published
   `SourceData_*.csv` files **do not move** — the CSVs already sent out stay
   valid. A third notation change is one edit to this file.
2. **The two reds in `CONTRAST_COL` were exchanged** (the 1 Sep decision). One
   rule now holds across the whole figure set: **pale is the factor on its own,
   saturated is that factor on the other's background.** Only one current panel
   draws all four colours, so only one panel changes appearance from this —
   Supp lipid c.
3. **`legend.text` gained `hjust = 0` in `theme_nc()`.** Parsed labels differ in
   width and ggplot sizes every legend cell to the widest, so without it the
   short label drifts off its own swatch. Two scripts that override
   `legend.text` locally got the same treatment.

**A second round, 8 September: group labels.** A *contrast* is a comparison; a
*group* is one arm of the study, and the two had drifted apart — Fig 2d said
"Variant / Variant + alcohol", the Fig 4d boxplots said "No alcohol / Alcohol",
and Fig 2f said "WD / WD+var / WD+alc / WD+var+alc". Three vocabularies for the
same mice. `nc_style.R` now carries a group layer too: `ALC_TICKS`, `alc_axis()`
and `ARM_PLOTMATH` / `nc_arm_expr()`.

Two decisions inside it, both forced by measurement rather than taste:

- **The two-group boxplots name the background once, in the axis title.** Ticks
  are `-` and `+`; the title reads "Alcohol (all mice PNPLA3^I148M)". Printing
  the full name on both ticks collides at 8 pt — six facets across 180 mm is
  15 mm a label and the name needs about 15 mm on its own. Shortening "alcohol"
  to "Alc." does not help: the width is set by the gene name, on the other line.
  The ticks are the ASCII hyphen and plus, never U+2212 — a true minus is the
  same class of character as the Unicode `>=` that came out as an empty box in
  Illustrator.
- **Fig 2f drops "WD" from its arm names.** Every arm is on Western diet, so
  inside the panel it carries no information; spelled out in full the labels
  truncate mid-superscript even at the 6 pt floor. The arms read Control /
  PNPLA3^I148M / Alcohol / PNPLA3^I148M + alcohol, and the diet is stated in
  the legend.

Spacing uses `*" ("`, not `~`. Measured at real 6 pt: the quoted space is
0.589 mm, one normal word space; `~` is 1.100 mm, nearly twice as wide and the
cause of the gap in the 31 Aug preview. **`notation_preview/render.R` still has
the `~` form — it is the wording reference, not the code reference.**

---

## 2 · Panels to re-render — run these fourteen scripts

Each writes into its figure's `nc_panels/`. Run order does not matter; none of
them feeds another.

| # | script | panels it rewrites | what moved |
|---|---|---|---|
| 1 | `Fig 1/fig1_volcano_bars_nc.R` | **Fig 1c** `Fig1c_classcounts_both_contrast.svg` | key labels |
| 2 | `Fig 1/fig1d_lsea_classes_nc.R` | **Fig 1d** `Fig1d_lsea_classes.svg` | y axis |
| 3 | `Fig 2/fig2_abc_heatmaps_nc.R` | **Fig 2a, 2b, 2c** | y axis, all three |
| 4 | `Fig 2/fig2_d_lipotoxic_boxplots_nc.R` | **Fig 2d** `Fig2d_lipotoxic.svg` | ticks + axis title |
| 5 | `Fig 2/fig2_e_f_pc_composition_nc.R` | **Fig 2e** `Fig2e_PC_unsaturation.svg` **and Fig 2f** `Fig2f_PC_saturation_index.svg` | x axes **and both panel widths** |
| 6 | `Fig 3/fig3_grid_nc.R` | **Fig 3a, 3b, 3c** | facet strips, all three |
| 7 | `Fig 4/fig4_cyp_slim_panel_nc.R` | **Fig 4e** `Fig4e_cyp_strip_periportal.svg` | the two arm rows |
| 8 | `Fig 5/fig5_heatmap_nc.R` | **Fig 5a** `Fig5_integrated_pathways.svg` | x axis, two-line |
| 9 | `Fig 5/fig5_aminoacid_panel_nc.R` | **Fig 5b** `Fig5b_aminoacid_groups.svg` | caption, **6 -> 8 pt labels, 72 -> 84 mm** |
| 10 | `supp_lipid_story/supp_a_lsea_nc.R` | **Supp lipid a** | y axis |
| 11 | `supp_lipid_story/supp_b_slopes_nc.R` | **Supp lipid b** | y axis |
| 12 | `supp_lipid_story/supp_c_routing_nc.R` | **Supp lipid c** | key labels **and the two reds** |
| 13 | `supp_lipid_story/supp_d_enzymes_nc.R` | **Supp lipid d** | y axis |
| 14 | `Fig 4/fig4_single_population.R` | **the Fig 4d boxplots** (one file per pathway) | ticks + axis title |

Fig 2e and Fig 5a use the **two-line** form: four labels share the panel width
there and the one-line form needs about 32 mm each at 8 pt. Everywhere else the
label sits on a y axis or a facet strip and one line fits.

**Eighteen panel files plus the Fig 4d boxplot set, from fourteen scripts.**
Seven were run end to end here against the real inputs, and their SVG and PNG checked by eye:
Fig 1d, Fig 2a-c, Fig 3a-c, Fig 5b, and Supp lipid a, b and c. The rest read
data that lives outside the connected folder, so they were syntax-checked and,
where the change was a layout one, geometry-tested at true size with stand-in
data — Fig 2d, Fig 2e and Fig 2f were all sized that way. Everything parses
clean.

---

## 3 · Panels that do NOT change

Worth knowing so you do not re-place them for nothing.

| panel | why |
|---|---|
| **Fig 1a, 1b** volcanoes | no contrast name is drawn in R. If the .ai carries one as typeset text, that text needs the new wording — check the file, not the script |
| **Fig 4b, 4c** | colour is pathway category; no contrast label |
| **Fig 4d** KEGG maps | the maps themselves carry no contrast label; check the run's `INDEX_pathviews.csv` for the old wording. The BOXPLOTS beneath them DO change — row 14 above |

---

## 4 · Then Illustrator — the slow part

Re-place the sixteen SVGs above in `Fig1_*.ai`, `Fig2_*.ai`, `Fig3_*.ai`,
`Fig4_*.ai`, `Fig5_*.ai` and the supplementary lipid-story file. Panel letters
are placed in Illustrator (`NC_TAGS <- FALSE`), so they do not come back with
the panel and have to survive the swap.

**Two things to watch while placing.**

- **The y-label gutter grew on the 88 mm panels.** Fig 1d, Supp lipid b and
  Fig 2b/c spend a little more width on the labels than they did — measured at
  about 6 mm on Fig 1d in a wide test font, less in Arial. The plot area is
  correspondingly narrower. It looks fine; it is just not pixel-identical to
  what is in the .ai now, so do not eyeball-match to the old placement.
- **The Fig 2 e/f row changes shape: 108 + 72, not 120 + 60.** Panel f took
  the 12 mm because its four arm names need about 18 mm each; panel e's
  four two-line contrast labels still clear at 27 mm each. The row still totals
  180 mm, so nothing else in Figure 2 moves.
- **Fig 5b grows 72 -> 84 mm** for the 8 pt metabolite names. Figure 5 then
  totals 219 mm against the 225 mm maximum, leaving 6 mm for the gutter and the
  panel letters. If that is too tight, the lever is panel a's row height in
  `fig5_heatmap_nc.R` (`H <- 26 + 5.2 * n`) — 5.0 mm per row gives back 4 mm.
- **Fig 1c and Fig 1d are a matched pair** — the same height by design, with
  the constants in the two scripts referring to each other. Re-place both.

---

## 5 · Not figures, but the same change — for another day

- **The eight approved legends** in `Figures_v1_editGB1NatComm_legends_Aug31.pptx`
  name the comparisons in text.
- **Supp Fig 1** is the study-design figure whose whole purpose was defining
  C1-C4. With the numbering gone it names the four comparisons in plain
  language — and it should now use the approved wording, superscript included.
  That is a rewrite, not an edit (see `HANDOFF.md` §6).
- **Fig 4a**, the pipeline diagram, is a PowerPoint export. Check it for the
  old wording.
- **The manuscript body**, per `TODO_manuscript_text_updates.md` §2.
- **The repo.** `pnpla3-alcohol-liver/analysis/` holds byte-identical copies of
  all thirteen edited files. They have been updated in step, so the two trees
  have not diverged, but nothing is committed — the commit and the push are
  yours.

---

## 6 · One thing that is still open

The **zip for the co-authors cannot be built until the renders and the
Illustrator placement are done** — what exists today is the code that produces
the new panels, not the assembled figures. When the .ai files are updated, the
zip is the exported PDFs plus the `SourceData_*.csv` set, which is unchanged.
