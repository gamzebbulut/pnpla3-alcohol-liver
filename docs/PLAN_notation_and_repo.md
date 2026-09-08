# Plan — the notation change and the reproducibility repo

Written 31 August 2026, after reading the code rather than guessing at it.
Companion to `TODO_manuscript_text_updates.md`, which stays the list of text
fixes; this file covers only the two new pieces of work.

**Decisions already made:** notation moves to the `PNPLA3-I148M` form, not
`-mut`, keeping the Banini convention the eight legends now use. Repo is private
now, public at acceptance. **SETTLED 8 Sep 2026:** superscript `I148M`, and Task 1 is done — see
`figure making/FIGURE_REMAKE_TASKLIST.md`. (Was open: the exact glyph — superscript `mut`,
superscript `m`, or spelled-out `I148M`. spelled-out `I148M`. Nothing in Task 1 started until that was settled; see 1.6.)

---

# Task 1 · The contrast notation

## 1.1 What the reconnaissance found

The four labels are **not just display strings — they are the internal keys**
the scripts index by.

| where | how the string is used |
|---|---|
| `nc_style.R` (11 hits) | `CONTRAST_LABEL`, `CONTRAST_SHORT`, `CONTRAST_COL`, `CONTRAST_ORDER` |
| `supp_d_enzymes_nc.R` | names of the `CONTRASTS` vector — each maps to an edgeR contrast formula |
| `fig2_e_f_pc_composition_nc.R` | names of the `CON` list — each maps to a pair of group codes |
| `fig3_grid_nc.R` | the *values* of a `MAP` from Q-codes to labels |
| `fig1_volcano_bars_nc.R` | `PAIR_ORDER`, `ORDER_BY`, panel spec `name` |
| 9 `SourceData_*.csv` | a `cmp` or `contrast` **column value** in the published data |

So a blind find-and-replace would rename the internal keys and the published
Source Data columns at the same time. It would work — everything is consistent —
but it is fragile, and it means the Source Data files a reviewer downloads change
their contents for a cosmetic reason.

## 1.2 The recommendation: separate the key from the label

Add one display layer to `nc_style.R` and leave the keys alone:

```r
## Internal keys never change. Only what the reader sees changes.
CONTRAST_DISPLAY <- c(
  "Alcohol on wild type" = "alcohol (WT bg)",
  "Alcohol on variant"   = "alcohol (PNPLA3-I148M bg)",
  "Variant alone"        = "PNPLA3-I148M (no-alcohol bg)",
  "Variant on alcohol"   = "PNPLA3-I148M (alcohol bg)")

nc_label <- function(x) unname(CONTRAST_DISPLAY[as.character(x)])
```

Every panel script then calls `nc_label()` at the one place the label is drawn —
a facet labeller, a manual scale's `labels =`, or an axis. Roughly one to three
lines per script.

**Why this rather than a rename.** This is the *second* notation change (C1–C4
went to plain language in August). A third is likely. With this layer the third
change is one edit to one file and a re-render; without it, it is 24 scripts
again. It also leaves the Source Data columns stable, so the CSVs already sent
out stay valid.

## 1.3 The superscript problem — read this before choosing the glyph

`PNPLA3` with a superscript is not a string change. In the legends it is a Word
run property; in a figure it has to be **plotmath**, e.g.
`expression(PNPLA3^I148M~"(alcohol bg)")`.

Two facts from the codebase:

- `label_parsed` is **already used** in `fig4_boxplots_nc.R` and
  `fig5_aminoacid_boxplots_nc.R`, so the technique is established and needs no
  new package.
- `ggtext` was **deliberately removed** from `fig5_heatmap_nc.R` ("ggtext
  dropped, and with it the `install.packages()` call that ran on load"). Do not
  reintroduce it — plotmath is the house route.

Plotmath is fussy: every label in a scale must become an expression, spaces need
`~` or quoting, and legend keys and facet strips take different code paths.
Budget for it rather than being surprised by it.

## 1.4 Blast radius

| what | count | notes |
|---|---|---|
| `nc_style.R` | 1 | the display layer + plotmath variants |
| scripts with hardcoded labels | 10 | listed in 1.1; the ones needing edits |
| scripts sourcing `nc_style.R` | 24 live | the rest inherit the change free |
| panels to re-render and re-place | ~20 | Illustrator placement is the slow part |
| `SourceData_*.csv` | 9 | **unchanged** under the 1.2 approach |
| approved legends | 8 | text mentions of the comparisons |
| review deck | 8 slides | re-fill and re-stamp |
| manuscript body | many | folds into `TODO` §2 |

## 1.5 Order of work

1. Settle the glyph (1.6). **Blocking.**
2. Add `CONTRAST_DISPLAY` + `nc_label()` to `nc_style.R`, with a plotmath variant
   for the superscript form.
3. Convert one script end to end — `supp_c_routing_nc.R`, small and faceted — and
   look at the rendered SVG before touching anything else.
4. Roll through the remaining nine, one commit each.
5. Re-render all panels, diff old vs new SVG text to confirm only labels moved.
6. Re-place in Illustrator. Slowest step, no way around it.
7. Legends, deck, then manuscript.

## 1.6 The wording — proposed, out for confirmation (31 Aug)

Settled on this set after rendering it at real size. **Sent to the group for
sign-off; nothing changes until they reply.**

| internal key (unchanged) | proposed label |
|---|---|
| Alcohol on wild type | Alcohol alone |
| Alcohol on variant   | Alcohol (PNPLA3^I148M bg) |
| Variant alone        | PNPLA3^I148M alone |
| Variant on alcohol   | PNPLA3^I148M (alcohol bg) |

"Alone" is exact, not shorthand — alcohol on the wild-type background *is*
alcohol alone. The set reads as two pairs: one insult on two backgrounds, then
one variant on two backgrounds.

**The superscript is load-bearing, not cosmetic.** Written flat as
`PNPLA3-I148M`, the labels clip off the right edge of an 88 mm panel. Rendered
evidence in `figure making/notation_preview/COMPARE_panelc_labels.png`, form D.

Previews, and the script that makes them, are in
`figure making/notation_preview/`. They were rendered from
`SourceData_Supp_lipid_c_routing.csv`, so the geometry is the shipped figure's.

## 1.7 Two plotmath traps, already solved in the preview script

1. plotmath **swallows a leading space inside a quoted string** —
   `PNPLA3^"I148M"*" (bg)"` renders with no gap. Use `~` for the space.
2. Parsed labels **break legend alignment** — ggplot sizes legend cells to the
   widest label, so a short label drifts away from its own key. Fix with
   `legend.text = element_text(hjust = 0)`.

Both cost a render to find. Carry them into every script in the sweep rather
than rediscovering them ten times.

---

# Task 2 · The reproducibility repo

## 2.1 The architecture that makes this work

There are three tiers of input, and only the middle one has to be in the repo.

| tier | what | where now | in repo? |
|---|---|---|---|
| 0 · raw | `Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx`; `Data/pseudobulk_snRNA-seq_PNPLA3_alc/*.RDS` (17 cell types) | **outside** `writing skeleton` | no — deposit or cite |
| 1 · derived | `de_lipids_Q2_Q5.csv`, `lsea_*.csv`, `singlepop_*_genes.csv`, `shared_Cyps_*.csv`, the `GB_metabolomics_20260628` integration outputs | mixed | **yes** — this is the key move |
| 2 · figures | 24 `*_nc.R` scripts + `nc_style.R` + 19 `SourceData_*.csv` (392 KB total) | `writing skeleton/figure making` | yes |

**The point:** committing tier 1 means every figure reproduces from the repo
alone, with no raw data and no snRNA-seq objects. Whether Isabel releases the raw
files stops being a blocker for reproducibility — it only affects whether someone
can rebuild tier 1 themselves. That is the normal arrangement for this kind of
paper.

## 2.2 Two things already true that make this cheap

- **Paths are already portable.** 34 of 35 scripts read
  `Sys.getenv("LIPID_ROOT", <default>)`. One script hardcodes a path and needs
  the same treatment. Nothing else needs editing to run elsewhere.
- **`MANIFEST_all_figures.md` is already the README.** Panel → script → file →
  size, plus the upstream table. It needs converting, not writing.

## 2.3 Proposed layout

```
pnpla3-alcohol-liver/
  README.md              <- from MANIFEST_all_figures.md
  DATA.md                <- what tier 0 is, where it lives, how to request it
  R/
    nc_style.R
    fig1/ fig2/ fig3/ fig4/ fig5/ supp/
  data/
    derived/             <- tier 1, the CSVs the scripts read
    source_data/         <- tier 2, the 19 published SourceData_*.csv
  figures/               <- rendered SVG + PDF, committed so the repo is
                            checkable without running anything
  run_all.R              <- sources every script in dependency order
  sessionInfo.txt        <- exact package versions from the run that made the figures
  .gitignore             <- *.RDS, *.xlsx, .DS_Store, ~$* lock files
```

## 2.4 Order of work

1. `git init` locally, `.gitignore` first so no raw data is ever staged.
2. Copy tier 2 in, unchanged. Commit. This alone is a working repo.
3. Trace each script's inputs and copy tier 1 in beside them; fix the one
   hardcoded path.
4. Write `run_all.R` in dependency order — the upstream table in the MANIFEST
   already gives the order.
5. **Run it in a clean R session from a fresh clone.** This is the only real test
   and it will find missing inputs; expect one round of fixes.
6. `sessionInfo.txt`, README, DATA.md.
7. Create the private GitHub repo, push, invite Katia and Isabel.

## 2.5 One access problem to solve first

Only `writing skeleton` is connected to this session. Tier 0 and the
`GB_metabolomics_20260628` integration outputs sit **beside** it under the
project root, so I cannot see or copy them from here. To help with steps 3–5,
connect the project folder itself
(`~/Desktop/Lipidomics project with Dr Smirnova`) rather than the sub-folder.
Steps 1, 2, 4 and 6 can proceed without that.

---

# Sequencing for today

**Repo first, notation second.** The repo is unblocked and the notation is not —
it waits on the glyph, and starting it twice costs a day. The repo also has a
natural stopping point after step 2.4.2: a working, pushable repo that reproduces
every figure from committed CSVs, even if nothing else gets done.

A realistic half-day: repo through step 5, plus the `nc_style.R` display layer
and one converted script so the notation approach is proven before it is
repeated ten times.

**What would make today go wrong:** starting the Illustrator re-placement before
the glyph is fixed, or committing tier 1 before checking that nothing under
`GB_metabolomics_20260628` is a symlink into a location the repo cannot carry.
