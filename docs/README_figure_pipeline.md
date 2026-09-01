# Figure pipeline — what makes each panel, and what it reads

**Cleaned 2026-06-26.** The multi-omics engine (`Fig 4/Stouffer_rmd_pathview_6hep_clean.Rmd`)
is now a **data producer only**: it computes the CT-weighted Stouffer + Pmax integration and
**writes CSVs**. Every figure is drawn by a **standalone script** that reads those CSVs (or reads
the metabolomics / snRNA-seq source data directly). No figure depends on objects living in the
Rmd session anymore.

Contrast naming: **Q3 = C1, Q4 = C2, Q2 = C3, Q5 = C4.**

---

## Run order

1. **Knit the Stouffer Rmd once per contrast** (set `QUESTION` in the `config` chunk to
   `Q2 / Q3 / Q4 / Q5`). Each knit writes a folder `<Q>_StoufferPMAX_pooled/` of CSVs under
   `OUT_BASE` (see config). Figures C1–C4 need the matching `Q` folder.
2. **Run the figure scripts** (any order). Each is `source()`-able in R; set the path/`PMAX_CSV`
   at the top to point at the contrast you want.

---

## CSVs the Stouffer Rmd writes (per contrast folder)

| File | Made in | Used by |
|---|---|---|
| `genes_liver_ct_stoufferDec12.csv` | §5 gene arm | Supp Table 1 |
| `metabolite_limma_ZDec12.csv` | §6 metabolite arm | Supp Table 2 |
| `modalities_pmax_gene_plus_metDec12.csv` | §7 **integration (the key file)** | Fig 4B, Fig 4C/D, Fig 5, Supp Table 3 |
| `concordance_labels_FDRlt0.10.csv` | §8 | the 13-pathway record (Supp Table 3) |
| `gene_members_concordant_FDR0.10.csv` / `metabolite_members_concordant_FDR0.10.csv` | §9 | member reference |
| `met_mat_contrast.csv` | §6 (decoupling export) | **Fig 4C/D** (metabolite pathview colors) |
| `sample_groups_met.csv` | §6 (decoupling export) | **Fig 4C/D** (REF/ALT sample pick) |
| `kegg_pathway_gene_symbols.csv` | §4.2 (decoupling export) | **Fig 4C/D** (boxplot member genes) |
| `ct_weights_pooled.csv`, `ct_limma_all_pathways.csv`, `genes_ct_Z_wide.csv`, `ct_group_mean_props_*.csv` | §2–§5 | QC / supplementary |

---

## Figure → code → inputs → outputs

### Figure 2 — metabolite heatmaps + indices (C1–C4)  *(independent of Stouffer)*
| Panel | Script | Reads | Writes |
|---|---|---|---|
| Heatmaps | `Fig 2/fig2_metabolite_heatmaps_C1C4.R` | metabolomics per-metabolite FC CSVs | heatmap PNG/PDF |
| Slim boxplots | `Fig 2/fig2_C1_slim_boxplots.R` | `todos.xlsx` (Data per metabolite) | boxplot PNG/PDF |
| 3-group panels | `Fig 2/Fig2_metabolite_panels_3groups.R` | metabolomics CSV | panels |

### Figure 3 — chain length / saturation / LSEA  *(independent of Stouffer)*
| Panel | Script | Reads | Writes |
|---|---|---|---|
| LSEA engine | `Fig 3/lipid_LSEA_Q3_vs_Q4.R` | `todos.xlsx` | `de_lipids_Q2_Q5.csv`, LSEA tables |
| 2×4 bars | `Fig 3/fig2_2x4_bars.R` | sources LSEA engine | bars |
| Banini arrows | `Fig 3/fig2_banini_arrows.R` | sources LSEA engine | arrow plot |
| LSEA heatmap | `Fig 3/fig2_LSEA_slim_heatmap.R` | LSEA CSV | heatmap |

### Figure 4 — multi-omics integration
| Panel | Script | Reads | Writes |
|---|---|---|---|
| A schematic | (PowerPoint, manual) | — | — |
| B pathway bars | `Fig 4/fig4_pathway_bars.R` | `modalities_pmax…csv` | `Fig4_pathway_bars_C1.png/pdf` |
| C pathview maps | `Fig 4/fig4_consolidated_pathview_boxplots.R` | `met_mat_contrast.csv`, `sample_groups_met.csv`, `kegg_pathway_gene_symbols.csv`, `modalities_pmax…csv`, hepatocyte RDS | `mmu*_<Q>_6hep*.png` |
| D member boxplots | *(same script)* | *(same)* | `boxplots_mmu*_<Q>.png`, `allHep_limma_genes_<Q>.csv` |

`fig4_consolidated_pathview_boxplots.R` is **standalone**: it reads the four CSVs above from
`DATA_DIR = dirname(PMAX_CSV)`, recomputes the 6-hepatocyte limma itself, and errors clearly if
any input CSV is missing (re-knit the Rmd for that contrast to regenerate them).

### Figure 5 — cross-contrast Pmax heatmap
| Panel | Script | Reads | Writes |
|---|---|---|---|
| Heatmap | `Fig 5/heatmap.Rmd` | per-Q `modalities_pmax…csv` (all four contrasts) | cross-contrast heatmap |

### Supplemental Figure 3 — enzyme / chain-length / reconciliation  *(standalone)*
| Panel | Script | Reads | Writes |
|---|---|---|---|
| Enzyme heatmap (4 contrasts) | `suppfig3/fig2_enzyme_heatmap.R` | hepatocyte RDS | `fig2_enzyme_heatmap_4contrasts.png/pdf`, `enzyme_limma_allHep_4contrasts.csv` |
| TG vs PL reconciliation | `suppfig3/fig2_reconciliation_TG_vs_PL.R` | `de_lipids_Q2_Q5.csv` (from Fig 3 LSEA) | `fig2_reconciliation_TG_vs_PL.png/pdf` |
| (gene list) | `suppfig3/chainlength_genes.R` | — | sourced helper only |

### Supplemental Tables 1–3
`supp Table 1/SuppTable1_gene_arm_stouffer_top15_C1.docx` — built from
`genes_liver_ct_stoufferDec12.csv` (T1), `metabolite_limma_ZDec12.csv` (T2),
`modalities_pmax_gene_plus_metDec12.csv` (T3).

---

## Notes
- **Decoupling (Fig 4C/D):** previously the consolidated script required `met_mat`,
  `sample_groups_met`, and `kegg` to exist in the Rmd session. It now reads them from CSV, so the
  Rmd no longer `source()`s it — run it yourself after knitting.
- **Pathway selection vs display (decision #13):** the 13 pathways are *selected* from the
  whole-liver, all-cell-type integration; Fig 4C/D *displays* them with the 6-hepatocyte model to
  avoid the zonation confound. Disclose in Methods.
- **Pmax ≈ P_met:** for all 13 joint pathways `Pmax` equals the metabolite p-value — the gene arm
  is the significance gate, the metabolite arm sets the rank.
