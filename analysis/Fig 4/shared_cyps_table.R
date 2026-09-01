# =============================================================================
# Shared Cyp table — Arachidonic acid (mmu00590), Linoleic acid (mmu00591),
# PPAR signaling (mmu03320), contrast C1.
# Which Cyp genes are shared vs unique across the three pathways, with the
# 6-hepatocyte alcohol logFC/FDR. (PPAR is gene-only, so it is NOT in the
# member CSV of the 13 — we read full KEGG membership instead.)
#
# INPUTS (both written by the pipeline):
#   - kegg_pathway_gene_symbols.csv  (Stouffer Rmd §4.2; full KEGG sets,
#                                      independent of min_size)
#   - allHep_limma_genes_<Q>.csv     (fig4_consolidated; per-gene logFC/FDR)
# HOW TO RUN: set the two paths below, then source() this file.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})

BASE     <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova")
# Full KEGG membership. After re-knitting the cleaned Stouffer Rmd this also
# appears in the publication Q3 folder; until then use the exploratory copy
# (membership is identical — written before any min_size filtering).
KEGG_CSV <- file.path(BASE, "writing skeleton/figure making/exploratory_gene_minsize",
                      "Q3_minGene10_StoufferPMAX/kegg_pathway_gene_symbols.csv")
DE_CSV   <- file.path(BASE, "writing skeleton/figure making/Fig 4/fig4_outs/allHep_limma_genes_Q3.csv")
OUT_DIR  <- file.path(BASE, "writing skeleton/figure making/Fig 4")

paths   <- c(AA = "mmu00590", LA = "mmu00591", PPAR = "mmu03320")
PREFIX  <- "^Cyp"            # gene family to tabulate (change to "^Ugt", etc.)

## ---- membership matrix: Cyp x {AA, LA, PPAR} --------------------------------
mem <- read_csv(KEGG_CSV, show_col_types = FALSE) %>%
  filter(KEGG_id %in% paths, str_detect(symbol, PREFIX)) %>%
  mutate(pathway = names(paths)[match(KEGG_id, paths)]) %>%
  distinct(symbol, pathway) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = pathway, values_from = present, values_fill = FALSE)
for (p in names(paths)) if (!p %in% names(mem)) mem[[p]] <- FALSE   # ensure all 3 cols

## ---- join per-gene alcohol effect (6-hep limma) ----------------------------
de <- read_csv(DE_CSV, show_col_types = FALSE) %>%
  transmute(symbol = SYMBOL, logFC, FDR = adj.P.Val)

cyp <- mem %>%
  left_join(de, by = "symbol") %>%
  rowwise() %>%
  mutate(sets      = paste(c("AA","LA","PPAR")[c(AA, LA, PPAR)], collapse = "+"),
         n_pathway = sum(c(AA, LA, PPAR)),
         direction = case_when(is.na(logFC) ~ "not detected",
                               logFC < 0   ~ "down",
                               TRUE        ~ "up")) %>%
  ungroup() %>%
  arrange(is.na(logFC), logFC) %>%                 # most-down first; undetected last
  select(symbol, AA, LA, PPAR, sets, n_pathway, logFC, FDR, direction)

## ---- save + print -----------------------------------------------------------
write_csv(cyp, file.path(OUT_DIR, "shared_Cyps_AA_LA_PPAR_C1.csv"))

cat("\nCyps per pathway:  AA =", sum(cyp$AA), " LA =", sum(cyp$LA),
    " PPAR =", sum(cyp$PPAR), " | union =", nrow(cyp), "\n")
cat("Shared AA+LA+PPAR :", paste(cyp$symbol[cyp$n_pathway == 3], collapse = ", "), "\n")
cat("AA+PPAR (Cyp4a)   :", paste(cyp$symbol[cyp$AA & cyp$PPAR & !cyp$LA], collapse = ", "), "\n")
cat("AA+LA  (Cyp2c/2j) :", paste(cyp$symbol[cyp$AA & cyp$LA & !cyp$PPAR], collapse = ", "), "\n")
cat("PPAR-only (bile)  :", paste(cyp$symbol[cyp$PPAR & !cyp$AA & !cyp$LA], collapse = ", "), "\n\n")
print(cyp, n = Inf)
cat("\nSaved -> ", file.path(OUT_DIR, "shared_Cyps_AA_LA_PPAR_C1.csv"), "\n")
