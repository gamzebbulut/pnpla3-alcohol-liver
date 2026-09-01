# =============================================================================
# FIGURE 4 — CONSOLIDATED: 6-hepatocyte limma  ->  pathview + matching boxplots
# -----------------------------------------------------------------------------
# ONE model drives everything (resolves Fig4-Problems #1-#8, #15, #16):
#   * 6 hepatocyte zones, RUV log2cpm_rb, ~0 + group + zone,
#     duplicateCorrelation(mouse), eBayes(trend, robust)   [Steve's adapted limma]
#   * pathview gene colors  = that model's logFC  (ENTREZ-mapped, NA-guarded)
#   * boxplot significance  = that model's BH FDR (NOT Wilcoxon)
#   * pathway list          = FDR_Pmax < 0.10 from the integration (no hardcode)
#   * blue = down (#1565C0), red = up (#C62828) for genes AND metabolites
#   * ONE explicit OUT_DIR; pathview wd handled with on.exit
#
# RUN STANDALONE: source() this file directly. It reads everything it needs from
# the per-contrast folder the Stouffer Rmd writes (met_mat_contrast.csv,
# sample_groups_met.csv, kegg_pathway_gene_symbols.csv, modalities_pmax...csv =
# DATA_DIR = dirname(PMAX_CSV)) and computes the 6-hepatocyte limma itself.
# No in-session objects required. Set PMAX_CSV/OUT_DIR for the contrast you want.
#
# NOTE (#13): pathways are SELECTED from the whole-liver integration but
# DISPLAYED with the 6-hepatocyte model — disclose this in the methods.
# =============================================================================
suppressWarnings(suppressMessages({
  library(limma); library(dplyr); library(tibble); library(ggplot2)
  library(readr); library(patchwork)
  library(pathview); library(org.Mm.eg.db); library(AnnotationDbi)
}))

## ----------------------------- CONFIG ---------------------------------------
BASE      <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova")
CT_DIR    <- file.path(BASE, "Data/pseudobulk_snRNA-seq_PNPLA3_alc")
EXPR_SLOT <- "log2cpm_rb"                     # RUV-corrected log2 CPM
QUESTION  <- if (exists("QUESTION")) QUESTION else "Q3"   # = C1
REF       <- if (exists("REF")) REF else "MUT_WSD_NAL"    # MUT, no alcohol
ALT       <- if (exists("ALT")) ALT else "MUT_WSD_PAL"    # MUT + alcohol
# Authoritative integrated pathway table (KEGG_id, KEGG_title, FDR_Pmax, concordance):
PMAX_CSV  <- file.path(BASE, "GB_metabolomics_20260628/Q3_StoufferPMAX_pooled/modalities_pmax_gene_plus_metDec12.csv")  # canonical run (matches Fig 4/5)
# ── EXPLORATORY TOGGLE (Isabel's question): use ALL snRNA-seq cell types? ─────
USE_ALL_CELLS <- TRUE        # TRUE = all cell types (stellate/Kupffer/etc.); FALSE = 6 hepatocyte zones
CELL_TAG      <- if (USE_ALL_CELLS) "allcells" else "6hep"
PATHWAYS_KEEP <- NULL   # NULL = all significant pathways (full supp set); or c("mmu00590",...) to restrict
OUT_DIR   <- file.path(BASE, "writing skeleton/figure making/Fig 4",
                       if (USE_ALL_CELLS) "fig4_outs_ALLCELLS" else "fig4_outs")

HEP_ZONES <- c("Periportal hepatocyte","Mid-zone hepatocyte","Pericentral hepatocyte",
               "ISG expressing hepatocyte","Migratory_IR hepatocyte","Remodeling hepatocyte")
FDR_PATH  <- 0.10        # pathway selection (FDR_Pmax)
FDR_GENE  <- 0.10        # boxplot star: limma adj.P.Val (matches figure FDR<0.10)
# Pathview fills must be LIGHT so the black KEGG gene labels stay readable
# (light ends of a red/blue diverging scale; boxplots keep the saturated Fig-1 colors)
## Colours come from nc_style.R so the supplement cannot drift from Figures 1-3
## or from the six-hepatocyte main figure. Sourced defensively.
.nc <- file.path(path.expand("~/Desktop/Lipidomics project with Dr Smirnova"),
                 "writing skeleton/figure making/nc_style.R")
if (file.exists(.nc)) source(.nc) else { NC_RED <- "#C62828"; NC_BLUE <- "#1565C0"
                                         NC_RED_TINT <- "#E59E9E"; NC_BLUE_TINT <- "#95BAE2"
                                         NC_GREY <- "#BBBBBB"; NC_INK <- "#1A1A1A" }

## PATHVIEW fills: TINTS of the Figure 3 hues (see the main-figure script).
COL_DOWN  <- NC_BLUE_TINT; COL_MID <- "#FFFFFF"; COL_UP <- NC_RED_TINT
## BOXPLOT fills: grey = no alcohol, blue = alcohol, matching Figure 2d and the
## main figure. Not red/blue — that pair means decreased/increased in the map.
BOX_DOWN  <- NC_GREY;  BOX_UP  <- NC_BLUE
# Pathview color scheme:
#   "kegg_default" = pathview defaults (genes green/red, metabolites blue/yellow) — the OLD look co-authors saw
#   "manuscript"   = light blue=down / red=up for BOTH layers (matches Fig 1-3, readable labels)
COLOR_SCHEME <- "manuscript"
# TRUE  = native KEGG-map bitmap PNG (the look you approved; labels NOT editable).
# FALSE = vector PDF (Graphviz layout) with editable text + adjustable CEX,
#         openable/resizable in Illustrator (Save As SVG there). Needs Rgraphviz + KEGGgraph.
KEGG_NATIVE <- TRUE
CEX         <- 0.30   # label size — only affects the FALSE (vector PDF) render
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- inputs: read from the per-contrast folder written by the Stouffer Rmd ---
## STANDALONE: no in-session objects needed. The Stouffer Rmd writes these three
## CSVs into the same folder as PMAX_CSV (DATA_DIR). Just point DATA_DIR there.
DATA_DIR <- dirname(PMAX_CSV)

.need_csv <- c("met_mat_contrast.csv", "sample_groups_met.csv",
               "kegg_pathway_gene_symbols.csv", basename(PMAX_CSV))
.missing  <- .need_csv[!file.exists(file.path(DATA_DIR, .need_csv))]
if (length(.missing))
  stop("Missing input CSV(s) in DATA_DIR:\n  ", DATA_DIR,
       "\n  -> ", paste(.missing, collapse = ", "),
       "\n  Re-knit the Stouffer Rmd for this contrast to (re)write them.")

# met_mat: compounds x samples (rownames = KEGG compound IDs)
.mm <- readr::read_csv(file.path(DATA_DIR, "met_mat_contrast.csv"), show_col_types = FALSE)
met_mat <- as.matrix(.mm[, -1, drop = FALSE]); rownames(met_mat) <- .mm$cpd_id
storage.mode(met_mat) <- "double"

# sample_groups_met: tibble(sample, group)
sample_groups_met <- readr::read_csv(file.path(DATA_DIR, "sample_groups_met.csv"),
                                     show_col_types = FALSE)

# kegg$sets_symbol: rebuilt from long CSV (KEGG_id, symbol) -> named list of gene vectors
.ks  <- readr::read_csv(file.path(DATA_DIR, "kegg_pathway_gene_symbols.csv"),
                        show_col_types = FALSE)
kegg <- list(sets_symbol = split(.ks$symbol, .ks$KEGG_id))

## ===== 1. all-hepatocyte limma =============================================
.ct_name <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
hep <- if (USE_ALL_CELLS) rds else rds[.ct_name(rds) %in% HEP_ZONES]
stopifnot("No cell-type RDS found in CT_DIR" = length(hep) > 0)
cat(sprintf("Cell types used (%s): %d  [%s]\n", CELL_TAG, length(hep),
            paste(.ct_name(hep), collapse = ", ")))

mats <- list(); metas <- list()
for (fp in hep) {
  z <- .ct_name(fp); o <- readRDS(fp)
  xm <- as.matrix(o$dge_list_obj[[EXPR_SLOT]]); storage.mode(xm) <- "double"
  colnames(xm) <- gsub("_.*", "", colnames(xm))
  sm  <- o$dge_list_obj$samples
  s2g <- setNames(as.character(sm$cmm_var), gsub("_.*", "", as.character(sm$sample)))
  xm  <- xm[, colnames(xm) %in% names(s2g), drop = FALSE]
  colnames(xm) <- paste(z, colnames(xm), sep = "::")
  mats[[z]]  <- xm
  metas[[z]] <- data.frame(col = colnames(xm), zone = z,
                           mouse = sub(".*::", "", colnames(xm)),
                           group = unname(s2g[sub(".*::", "", colnames(xm))]),
                           stringsAsFactors = FALSE)
}
common <- Reduce(intersect, lapply(mats, rownames))
mat  <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
meta <- do.call(rbind, metas); rownames(meta) <- meta$col; meta <- meta[colnames(mat), ]
meta <- meta[!is.na(meta$group), ]; mat <- mat[, rownames(meta)]
meta$group <- factor(meta$group)
meta$zone  <- if (USE_ALL_CELLS) factor(meta$zone) else factor(meta$zone, levels = HEP_ZONES)
cat(sprintf("Stacked: %d genes x %d pseudo-samples (%d mice x %d zones)\n",
            nrow(mat), ncol(mat), length(unique(meta$mouse)), nlevels(meta$zone)))

design <- model.matrix(~ 0 + group + zone, data = meta)
colnames(design) <- make.names(sub("^group", "", colnames(design)))
corfit <- duplicateCorrelation(mat, design, block = meta$mouse)
cat(sprintf("duplicateCorrelation (mouse): %.3f\n", corfit$consensus.correlation))
fit <- lmFit(mat, design, block = meta$mouse, correlation = corfit$consensus.correlation)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
ctr <- makeContrasts(contrasts = paste0(make.names(ALT), " - ", make.names(REF)), levels = design)
de  <- topTable(eBayes(contrasts.fit(fit, ctr), trend = TRUE, robust = TRUE), number = Inf) |>
       tibble::rownames_to_column("SYMBOL")
write_csv(de, file.path(OUT_DIR, paste0("allHep_limma_genes_", QUESTION, ".csv")))
cat(sprintf("limma de: %d genes | Pla2g4a logFC=%s (FDR=%s)\n", nrow(de),
            { v <- de$logFC[de$SYMBOL=="Pla2g4a"]; if (length(v)) sprintf("%+.2f", v[1]) else "NA" },
            { v <- de$adj.P.Val[de$SYMBOL=="Pla2g4a"]; if (length(v)) signif(v[1],2) else "NA" }))

## ===== 2. gene.data for pathview (ENTREZ-named logFC) =======================
map <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(de$SYMBOL),
                             keytype = "SYMBOL", columns = "ENTREZID")
map <- map[!is.na(map$ENTREZID), ]; map <- map[!duplicated(map$SYMBOL), ]
de_e <- merge(de[, c("SYMBOL","logFC","adj.P.Val")], map, by = "SYMBOL")
gene.data <- setNames(de_e$logFC, de_e$ENTREZID); gene.data <- gene.data[!is.na(names(gene.data))]
cat(sprintf("gene.data: %d ENTREZ-mapped genes\n", length(gene.data)))
stopifnot("gene.data empty — ENTREZ mapping failed" = length(gene.data) > 0)

## ===== 3. cpd.data for pathview (metabolite log2FC, ALT - REF) ==============
# sample_groups_met is a tibble (cols: sample, group), not a named vector
pal <- sample_groups_met$sample[sample_groups_met$group == ALT]
nal <- sample_groups_met$sample[sample_groups_met$group == REF]
cpd_fc   <- log2(rowMeans(met_mat[, pal, drop=FALSE], na.rm=TRUE) + 1e-6) -
            log2(rowMeans(met_mat[, nal, drop=FALSE], na.rm=TRUE) + 1e-6)
cpd.data <- setNames(as.numeric(cpd_fc), rownames(met_mat))

## ===== 4. significant pathways (FDR_Pmax < 0.10, both arms) — no hardcode ===
sig <- readr::read_csv(PMAX_CSV, show_col_types = FALSE) |>
       dplyr::filter(FDR_Pmax < FDR_PATH, concordance != "One-missing") |>
       dplyr::arrange(FDR_Pmax)
if (!is.null(PATHWAYS_KEEP)) sig <- dplyr::filter(sig, KEGG_id %in% PATHWAYS_KEEP)
cat(sprintf("Pathways selected (%s): %d  [%s]\n", CELL_TAG, nrow(sig),
            paste(sig$KEGG_id, collapse = ", ")))
stopifnot("No significant pathways — check PMAX_CSV path" = nrow(sig) > 0)

## ===== 5. pathview maps (blue/red, fixed limit) ============================
# pick color scheme + a filename tag so the two versions don't overwrite each other
if (COLOR_SCHEME == "manuscript") {
  pv_low <- list(gene=COL_DOWN, cpd=COL_DOWN); pv_mid <- list(gene=COL_MID, cpd=COL_MID)
  pv_high<- list(gene=COL_UP,   cpd=COL_UP);   pv_tag <- CELL_TAG
} else {  # kegg_default — OLD colors (genes green/red, metabolites blue/yellow)
  pv_low <- list(gene="green", cpd="blue"); pv_mid <- list(gene="gray", cpd="gray")
  pv_high<- list(gene="red",   cpd="yellow"); pv_tag <- paste0(CELL_TAG, "_keggcol")
}

old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(OUT_DIR)
for (i in seq_len(nrow(sig))) {
  pid <- sub("^mmu", "", sig$KEGG_id[i])
  tryCatch(
    pathview(gene.data = gene.data, cpd.data = cpd.data, pathway.id = pid, species = "mmu",
             out.suffix = paste0(QUESTION, "_", pv_tag),
             kegg.native = KEGG_NATIVE,
             same.layer  = !KEGG_NATIVE,              # native PNG needs FALSE; Graphviz PDF uses TRUE
             cex = CEX,
             limit = list(gene = 3, cpd = 3),
             low = pv_low, mid = pv_mid, high = pv_high),
    error = function(e) message("  pathview mmu", pid, " FAILED: ", conditionMessage(e)))
}
setwd(old)

## ===== 6. per-mouse 6-hep expression (mean across zones) ====================
mouse_grp <- setNames(meta$group[!duplicated(meta$mouse)],
                      meta$mouse[!duplicated(meta$mouse)])
permouse <- sapply(names(mouse_grp), function(mo)
  rowMeans(mat[, rownames(meta)[meta$mouse == mo], drop = FALSE]))   # gene x mouse
grp_of <- as.character(mouse_grp[colnames(permouse)])

## ===== 7. boxplots per pathway (members present; show significant ones) =====
qlook <- setNames(de$adj.P.Val, de$SYMBOL)
gene_box <- function(sym) {
  v <- permouse[sym, ]; g <- grp_of
  keep <- g %in% c(REF, ALT); v <- v[keep]; g <- factor(g[keep], levels = c(REF, ALT))
  q <- qlook[sym]; star <- if (!is.na(q) && q < FDR_GENE) "*" else ""
  ggplot(data.frame(g, v), aes(g, v, fill = g)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.5) +
    geom_jitter(width = 0.12, size = 1.6, shape = 21, colour = "black", stroke = 0.3) +
    annotate("text", x = 1.5, y = max(v), label = star, size = 6,
             colour = BOX_UP, fontface = "bold", vjust = 0) +
    scale_fill_manual(values = setNames(c(BOX_DOWN, BOX_UP), c(REF, ALT))) +
    scale_x_discrete(labels = c("-EtOH", "+EtOH")) +
    labs(title = sym, x = NULL, y = paste0("log2CPM (", CELL_TAG, ")")) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold.italic", size = 10, hjust = 0.5))
}
for (i in seq_len(nrow(sig))) {
  syms <- kegg$sets_symbol[[ sig$KEGG_id[i] ]]
  if (is.null(syms)) next
  present <- intersect(syms, rownames(permouse))
  sigm <- present[vapply(present, function(s){q<-qlook[s]; !is.na(q) && q < FDR_GENE}, logical(1))]
  show <- if (length(sigm)) sigm else present     # significant-only if any, else all present
  if (!length(show)) next
  panels <- lapply(show, gene_box)
  nc <- min(8, length(panels))   # up to 8 per row (so the 6 AA CYPs sit on one row)
  fig <- patchwork::wrap_plots(panels, ncol = nc) +
    patchwork::plot_annotation(title = sprintf("%s — %s (C1; * = limma FDR<%.2f)",
                                               sig$KEGG_id[i], sig$KEGG_title[i], FDR_GENE))
  ggsave(file.path(OUT_DIR, sprintf("boxplots_%s_%s.png", sig$KEGG_id[i], QUESTION)),
         fig, width = nc*2.3, height = ceiling(length(panels)/nc)*2.5 + 0.5,
         dpi = 200, limitsize = FALSE)
}

cat("\nDONE -> ", OUT_DIR, "\n",
    "  pathview: mmu*.", QUESTION, "_6hep.png (verify Pla2g4a = red in mmu00590)\n",
    "  boxplots: boxplots_mmu*_", QUESTION, ".png (stars = limma FDR)\n", sep = "")
