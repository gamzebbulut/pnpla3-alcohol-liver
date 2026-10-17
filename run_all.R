## Rebuild every panel in the paper, in dependency order.
##
## Upstream scripts write the CSVs that panel scripts read, so they come first.
## Anything requiring raw data that is not distributed is skipped with a note
## rather than failed, so a fresh clone still produces the full figure set.
## Find setup.R beside this file, so run_all.R works from any working directory.
## %||% is not available in older R, hence the explicit test.
local({
  here <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { here <- dirname(normalizePath(of)); break }
  }
  if (is.null(here)) here <- normalizePath(".")
  source(file.path(here, "setup.R"))
})
FIG <- Sys.getenv("FIGROOT")

run <- function(rel) {
  f <- file.path(FIG, rel)
  if (!file.exists(f)) { message("skip (absent): ", rel); return(invisible(NULL)) }
  message("\n=== ", rel)
  tryCatch(source(f, chdir = TRUE, local = new.env()),
           error = function(e) {
             m <- conditionMessage(e)
             message("   FAILED: ", m)
             ## The two failures a clone can expect, named so nobody has to
             ## reverse-engineer them from a file path.
             if (grepl("EtOH.xlsx", m, fixed = TRUE))
               message("   -> needs the raw metabolomics workbook, which is not ",
                       "distributed. See DATA.md; set OWL_CACHE to run this panel ",
                       "from the cached sheets instead.")
             if (grepl("no package called", m))
               message("   -> missing R package. install.packages() or ",
                       "BiocManager::install() it, then re-run.")
           })
}

## -- panels ----------------------------------------------------------------
run("Fig 1/fig1_volcano_bars_nc.R")
run("Fig 1/fig1d_lsea_classes_nc.R")
run("Fig 2/fig2_abc_heatmaps_nc.R")
run("Fig 2/fig2_d_lipotoxic_boxplots_nc.R")
run("Fig 2/fig2_e_f_pc_composition_nc.R")
run("Fig 3/fig3_grid_nc.R")
run("Fig 4/fig4b_arm_bars_nc.R")
run("Fig 4/fig4_concordance_scatter_nc.R")
run("Fig 4/fig4_cyp_slim_panel_nc.R")
run("Fig 5/fig5_heatmap_nc.R")
run("Fig 5/fig5_aminoacid_panel_nc.R")
run("supp_lipid_story/supp_a_lsea_nc.R")
run("supp_lipid_story/supp_b_slopes_nc.R")
run("supp_lipid_story/supp_c_routing_nc.R")
run("supp_lipid_story/supp_d_enzymes_nc.R")

## -- supplementary tables (need openxlsx, flextable, officer) ----------------
run("supp_tables/supp_table1_alcohol_both_backgrounds.R")
run("supp_tables/supp_table2_alcohol_variant_only.R")
run("supp_tables/supp_table3_joint_pathways.R")

message("\ndone. panels are under analysis/<Fig N>/nc_panels/, tables under analysis/supp_tables/out/")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
