# Supplementary tables

| Table | Script | What it lists |
|---|---|---|
| Supplementary Table 1 | `supp_table1_alcohol_both_backgrounds.R` | the 9 metabolites altered by alcohol on both the PNPLA3 I148M and the wild-type background |
| Supplementary Table 2 | `supp_table2_alcohol_variant_only.R` | the 98 metabolites altered by alcohol on the PNPLA3 I148M background only |
| Supplementary Table 3 | `supp_table3_joint_pathways.R` | the 14 KEGG pathways perturbed in both the gene and metabolite arms (Fig. 4a rule) |

Tables 1 and 2 read `analysis/Fig 1/FC_Q3_metabolites.csv` and `FC_Q4_metabolites.csv`
(Welch P < 0.05 and Storey q < 0.10 in each comparison). Table 3 reads the integration output
`GB_metabolomics_20260628/Q3_StoufferPMAX_pooled/modalities_pmax_gene_plus_metDec12.csv`
verbatim. Each script stops if its counts drift from the manuscript (107 / 32 / 9 / 98; 14 pathways).

```r
source("setup.R")
for (f in list.files(file.path(Sys.getenv("FIGROOT"), "supp_tables"), "^supp_table.*\\.R$", full.names = TRUE))
  source(f, chdir = TRUE, local = new.env())
```

Output (`.csv`, `.xlsx`, `.docx`) goes to `analysis/supp_tables/out/`, or to `SUPPTAB_OUT` if set.
Needs `openxlsx`, `flextable` and `officer`.
