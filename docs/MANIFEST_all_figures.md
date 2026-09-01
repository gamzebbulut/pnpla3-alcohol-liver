# Figure code map — every panel, the script that makes it, and what it reads

Written 31 August 2026. One row per panel. If a script is not in the panel
tables below or the "Upstream" table, it does not make anything in the current
manuscript.

`nc_style.R` (this folder) holds every shared color, type size, contrast name
and the SVG writer. Every panel script sources it. Paths come from environment
variables with defaults, so nothing needs editing to run.

**Current panels live only in `<Fig N>/nc_panels/`.** That directory also still
holds older renders under panel letters that have since moved — see
"Stale renders" at the bottom before you place anything in Illustrator.

---

## Figure 1 — the metabolome, two alcohol contrasts

Output: `Fig 1/nc_panels/`

| panel | what it shows | script | file | mm |
|---|---|---|---|---|
| a | Volcano, alcohol on wild type | `fig1_volcano_bars_nc.R` | `Fig1a_volcano_alcohol_on_wild_type.svg` | 88 x 76 |
| b | Volcano, alcohol on variant | same script | `Fig1b_volcano_alcohol_on_variant.svg` | 88 x 76 |
| c | Significant metabolites by class, both contrasts | same script | `Fig1c_classcounts_both_contrast.svg` | 88 x 50 |
| d | LSEA class heatmap, both contrasts | `fig1d_lsea_classes_nc.R` | `Fig1d_lsea_classes.svg` | 88 x 50 |

Reads: `Fig 1/FC_Q3_metabolites.csv`, `FC_Q4_metabolites.csv` (a-c);
`lsea_ALLQ_comparison.csv` (d, written by `Fig 3/lipid_LSEA_Q3_vs_Q4.R`).

c and d are deliberately the same height so the bottom row is even; the
constants reference each other in comments.

## Figure 2 — lipid classes and PC saturation

Output: `Fig 2/nc_panels/`

| panel | what it shows | script | file | mm |
|---|---|---|---|---|
| a | Glycerophospholipid heatmap, 4 contrasts | `fig2_abc_heatmaps_nc.R` | `Fig2a_glycerophospholipids.svg` | 180 x 44 |
| b | Sphingo- + glycerolipids | same script | `Fig2b_sphingo_glycerolipids.svg` | 88 x 44 |
| c | Amino acids + other classes | same script | `Fig2c_aminoacids_other.svg` | 88 x 44 |
| d | Six lipotoxic readouts, boxplots | `fig2_d_lipotoxic_boxplots_nc.R` | `Fig2d_lipotoxic.svg` | 180 x 42 |
| e | PC fold-change by unsaturation bin | `fig2_e_f_pc_composition_nc.R` | `Fig2e_PC_unsaturation.svg` | 120 x 62 |
| f | PC saturation index, one point per mouse | same script | `Fig2f_PC_saturation_index.svg` | 60 x 62 |

Reads: `Fig 1/FC_Q2..Q5_metabolites.csv` (a-c);
`Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx` (d-f, or a CSV cache via `OWL_CACHE`).

The three a-c panels come from one `GROUPS` list in one script — the split into
180 / 88 / 88 is a layout decision in that list, not three scripts.

## Figure 3 — chain length and unsaturation

Output: `Fig 3/nc_panels/`

| panel | what it shows | script | file | mm |
|---|---|---|---|---|
| a | Fold-change by total double bonds | `fig3_grid_nc.R` | `Fig3a_unsaturation.svg` | 180 x 52 |
| b | Mean fold-change by chain length | same script | `Fig3b_chainlength.svg` | 180 x 52 |
| c | Per-species triglycerides | same script | `Fig3c_TG_species.svg` | 180 x 58 |

Reads: `Fig 3/outs/de_lipids_Q2_Q5.csv`, `Fig 3/outs/lsea_results_Q2_Q5.csv`
(both written by `Fig 3/lipid_LSEA_Q3_vs_Q4.R`).

## Figure 4 — multi-omics integration

Output: `Fig 4/nc_panels/` and `Fig 4/pathview_outs/`

| panel | what it shows | script | file | mm |
|---|---|---|---|---|
| a | Pipeline diagram | none — PowerPoint export | `Fig4a_pipeline.pdf` | placed in Illustrator |
| b | Gene arm vs metabolite arm Z, joint hits | `fig4_concordance_scatter_nc.R` | `Fig4b_concordance_scatter.svg` + `Fig4_category_key.svg` | 74 x 86, key 180 x 8 |
| c | Gene arm ranking bars | `fig4b_arm_bars_nc.R` | `Fig4c_gene_arm_slim.svg` | 54 x 86 |
| d | KEGG pathway maps, periportal | `fig4_single_population.R` | `pathview_outs/20260830_periportal/native_png/`, `.../vector_pdf/` | KEGG raster / Graphviz |
| e | Cyp pathway membership + periportal DE | `fig4_cyp_slim_panel_nc.R` | `Fig4e_cyp_strip_periportal.svg` | 170 x 45 |
| — | Boxplots beneath the maps | `fig4_single_population.R` | `pathview_outs/20260830_periportal/boxplots/` | 180 x 80 |
| supp | Metabolite arm ranking | `fig4b_arm_bars_nc.R` | `SuppFig_metabolite_arm.svg` | 88 x 86 |

Reads: `GB_metabolomics_20260628/Q3_StoufferPMAX_pooled/*.csv` (b, c, d);
`Fig 4/pathway_categories.csv` (b, c); `Fig 4/C1_joint_hits.csv` (b);
`Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv` and `fig4_outs/cyp_specificity_qlf_long.csv` (e);
`Data/pseudobulk_snRNA-seq_PNPLA3_alc/*.RDS` (d and the boxplots).

Panel d and the boxplots come from ONE script and ONE test — periportal
hepatocytes, edgeR quasi-likelihood F — so nothing in the figure can disagree
with anything else in it. `pathview_outs/README.md` explains the run folders,
the color scheme and the PPAR-signaling extra; each run folder carries an
`INDEX_pathviews.csv` naming its provenance.

**`Fig 4/MANIFEST.md` is out of date.** It still lists
`fig4_animal_level_limma.R` and `fig4_cyp_panel_nc.R` as trusted, from before
Fig 4 moved to the single-population route.

## Figure 5 — integration heatmap and amino acids

Output: `Fig 5/nc_panels/`

| panel | what it shows | script | file | mm |
|---|---|---|---|---|
| a | Cross-contrast pathway heatmap, 21 pathways | `fig5_heatmap_nc.R` | `Fig5_integrated_pathways.svg` | 180 x 135 |
| b | 37 amino-acid + sulfur metabolites by pathway | `fig5_aminoacid_panel_nc.R` | `Fig5b_aminoacid_groups.svg` | 180 x 72 |

Reads: `GB_metabolomics_20260628/Q*_StoufferPMAX_pooled/*.csv` (a);
`Fig 1/FC_Q3_metabolites.csv` (b).

`fig5_aminoacid_boxplots_nc.R` is the earlier per-metabolite boxplot version of
panel b (180 x 56). Kept, not placed.

## Supplementary figure — the lipid story

Combines the former Supp Fig 2 and Supp Fig 3. Output:
`supp_lipid_story/nc_panels/`. Panel letters go on in Illustrator
(`NC_TAGS <- FALSE`), not in R.

| panel | what it shows | script | file | mm |
|---|---|---|---|---|
| a | LSEA across all four contrasts | `supp_a_lsea_nc.R` | `Supp_lipid_a_lsea.svg` | 180 x 46 |
| b | Slope of fold-change on double bonds | `supp_b_slopes_nc.R` | `Supp_lipid_b_slopes.svg` | 88 x 58 |
| c | Lipid routing, grouped bars | `supp_c_routing_nc.R` | `Supp_lipid_c_routing.svg` | 88 x 58 |
| d | Enzyme heatmap, periportal edgeR QLF | `supp_d_enzymes_nc.R` | `Supp_lipid_d_enzymes.svg` | 180 x 44 |

Layout: a (46) + [b | c] (58) + d (44) = 156 mm. See
`supp_lipid_story/LAYOUT_Supp_lipid_story.md`.

---

## Upstream scripts — these write CSVs, not panels

| script | writes | feeds |
|---|---|---|
| `Fig 3/lipid_LSEA_Q3_vs_Q4.R` | `de_lipids_Q2_Q5.csv`, `lsea_results_Q2_Q5.csv`, `lsea_ALLQ_comparison.csv` | Fig 1d, Fig 3a-c, Supp a/b/c |
| `Fig 4/shared_cyps_table.R` | `shared_Cyps_AA_LA_PPAR_C1.csv` | Fig 4e |
| `Fig 4/fig4_cyp_specificity_qlf.R` | `fig4_outs/cyp_specificity_qlf_long.csv` | Fig 4e |
| `Fig 4/fig4_single_population.R` | `fig4_outs/singlepop_*_genes.csv`, `pathview_outs/<run>/INDEX_pathviews.csv` | Fig 4d, Fig 4e, the boxplots |

Method checks — run to answer a reviewer, not to build a panel:
`Fig 4/fig4_steve_edgeR_qlf.R` (per-population QLF, the sanctioned route),
`Fig 4/fig4_animal_level_limma.R` (animal-level limma, the earlier route),
`Fig 4/fig4_ppara_programme_and_qc.R` (PPAR-alpha programme + design QC),
`Fig 4/fig4_explore_two_arms.R` (why the two alcohol arms differ).

---

## Stale renders sitting in `nc_panels/` — check before placing

These are older files under panel letters that have since moved. The risk is
real: in Figure 1 the same letter exists twice with the contrasts swapped.

| file | why it is stale |
|---|---|
| `Fig 1/nc_panels/Fig1a_volcano_alcohol_on_variant.svg` | before the 31 Aug left/right swap — a is now wild type |
| `Fig 1/nc_panels/Fig1b_volcano_alcohol_on_wild_type.svg` | same swap — b is now variant |
| `Fig 1/nc_panels/Fig1b_classcounts_alcohol_on_variant.svg` | 29 Aug four-panel layout, before the merge |
| `Fig 1/nc_panels/Fig1c_volcano_alcohol_on_wild_type.svg` | 29 Aug layout |
| `Fig 1/nc_panels/Fig1d_classcounts_alcohol_on_wild_type.svg` | 29 Aug layout |
| `Fig 4/nc_panels/Fig4d_cyp_panel.svg` | earlier Cyp panel |
| `Fig 4/nc_panels/Fig4e_cyp_panel.svg` | earlier Cyp panel, relettered |
| `Fig 4/nc_panels/Fig4e_cyp_boxplots.svg` | six-hepatocyte boxplots, replaced by the single-population set |
| `Fig 4/nc_panels/Fig4e_cyp_population_grid.svg` | population grid, dropped when Cyps were parked |
| `Fig 4/nc_panels/Explore_a_*.svg`, `Explore_b_*.svg` | exploratory, never a panel |
| `supp_lipid_story/nc_panels/Supp_lipid_b_reconciliation.svg` | merged into `Supp_lipid_b_slopes.svg` |

Their scripts are still around and still run; nothing here is broken. They are
just not the current panels.
