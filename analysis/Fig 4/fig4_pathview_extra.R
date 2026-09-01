# =============================================================================
# EXTRA pathview maps — gene-only pathways NOT in the integrated 13
# -----------------------------------------------------------------------------
# Renders a KEGG pathview for any pathway you list (e.g. PPAR signaling), colored
# by the SAME 6-hepatocyte limma logFC used in Fig 4C. Reads the gene table the
# consolidated Fig 4 script already wrote (allHep_limma_genes_<Q>.csv) — it does
# NOT recompute limma and does NOT touch the publication Fig 4 script.
#
# Use case: PPAR signaling (mmu03320) is the strongest gene-arm pathway but is
# gene-only (no metabolite corroboration), so it is absent from the 13 integrated
# maps. This renders it for a side-by-side "gene-specific highlights" panel.
#
# HOW TO RUN: set QUESTION + PATHWAYS below, then source() this file.
# =============================================================================
suppressWarnings(suppressMessages({
  library(readr); library(dplyr)
  library(pathview); library(org.Mm.eg.db); library(AnnotationDbi)
}))

## ----------------------------- CONFIG ---------------------------------------
BASE      <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova")
QUESTION  <- "Q3"                                   # = C1
REF       <- "MUT_WSD_NAL"; ALT <- "MUT_WSD_PAL"

# Gene table written by fig4_consolidated_pathview_boxplots.R (SYMBOL, logFC, ...)
GENE_CSV  <- file.path(BASE, "writing skeleton/figure making/Fig 4/fig4_outs",
                       paste0("allHep_limma_genes_", QUESTION, ".csv"))

# OPTIONAL metabolite layer. If these exist (written by the cleaned Stouffer Rmd),
# the cpd nodes get colored too; otherwise the map is rendered gene-only.
STOUFFER_DIR <- file.path(BASE, "GB_metabolomics_20260628/Q3_StoufferPMAX_pooled")
MET_CSV      <- file.path(STOUFFER_DIR, "met_mat_contrast.csv")
SGM_CSV      <- file.path(STOUFFER_DIR, "sample_groups_met.csv")

# Pathways to render: named vector  KEGG_id = title  (add more if you want)
PATHWAYS  <- c(mmu03320 = "PPAR signaling pathway")

OUT_DIR   <- file.path(BASE, "writing skeleton/figure making/Fig 4/fig4_outs/extra_pathviews")
COLOR_SCHEME <- "kegg_default"   # match Fig 4C: genes green/red, metabolites blue/yellow
KEGG_NATIVE  <- TRUE
CEX          <- 0.30
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ----------------------------- gene.data ------------------------------------
stopifnot("Gene CSV not found — run fig4_consolidated first" = file.exists(GENE_CSV))
de <- readr::read_csv(GENE_CSV, show_col_types = FALSE)
map <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(de$SYMBOL),
                             keytype = "SYMBOL", columns = "ENTREZID")
map <- map[!is.na(map$ENTREZID), ]; map <- map[!duplicated(map$SYMBOL), ]
de_e <- merge(de[, c("SYMBOL","logFC")], map, by = "SYMBOL")
gene.data <- setNames(de_e$logFC, de_e$ENTREZID)
gene.data <- gene.data[!is.na(names(gene.data))]
cat(sprintf("gene.data: %d ENTREZ-mapped genes\n", length(gene.data)))

## ----------------------------- cpd.data (optional) --------------------------
cpd.data <- NULL
if (file.exists(MET_CSV) && file.exists(SGM_CSV)) {
  .mm <- readr::read_csv(MET_CSV, show_col_types = FALSE)
  met_mat <- as.matrix(.mm[, -1, drop = FALSE]); rownames(met_mat) <- .mm$cpd_id
  sgm <- readr::read_csv(SGM_CSV, show_col_types = FALSE)
  pal <- sgm$sample[sgm$group == ALT]; nal <- sgm$sample[sgm$group == REF]
  cpd_fc <- log2(rowMeans(met_mat[, pal, drop=FALSE], na.rm=TRUE) + 1e-6) -
            log2(rowMeans(met_mat[, nal, drop=FALSE], na.rm=TRUE) + 1e-6)
  cpd.data <- setNames(as.numeric(cpd_fc), rownames(met_mat))
  cat(sprintf("cpd.data: %d metabolites\n", length(cpd.data)))
} else {
  cat("cpd.data: none found -> rendering gene-only (fine for PPAR signaling)\n")
}

## ----------------------------- colors ---------------------------------------
if (COLOR_SCHEME == "manuscript") {
  pv_low <- list(gene="#92C5DE", cpd="#92C5DE"); pv_mid <- list(gene="#FFFFFF", cpd="#FFFFFF")
  pv_high<- list(gene="#F4A582", cpd="#F4A582"); pv_tag <- paste0(QUESTION, "_extra")
} else {  # kegg_default
  pv_low <- list(gene="green", cpd="blue"); pv_mid <- list(gene="gray", cpd="gray")
  pv_high<- list(gene="red",   cpd="yellow"); pv_tag <- paste0(QUESTION, "_extra_keggcol")
}

## ----------------------------- render ---------------------------------------
old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(OUT_DIR)
for (kid in names(PATHWAYS)) {
  pid <- sub("^mmu", "", kid)
  cat(sprintf("Rendering %s — %s ...\n", kid, PATHWAYS[[kid]]))
  tryCatch(
    pathview(gene.data = gene.data, cpd.data = cpd.data,
             pathway.id = pid, species = "mmu",
             out.suffix = pv_tag, kegg.native = KEGG_NATIVE,
             same.layer = !KEGG_NATIVE, cex = CEX,
             limit = list(gene = 3, cpd = 3),
             low = pv_low, mid = pv_mid, high = pv_high),
    error = function(e) message("  pathview mmu", pid, " FAILED: ", conditionMessage(e)))
}
setwd(old)
cat("\nDONE -> ", OUT_DIR, "\n",
    "  files: mmu", paste(sub("^mmu","",names(PATHWAYS)), collapse=", "),
    ".", pv_tag, ".png\n", sep = "")
