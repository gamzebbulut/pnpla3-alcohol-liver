# =============================================================================
#  Is the Cyp hit specific to the variant background?  —  all cell types
#
#  THE QUESTION, AND THE ONLY HONEST WAY TO ASK IT
#  "The Cyps go down when alcohol is given on the variant background, and they
#  do not on the wild-type background" is NOT a demonstration of specificity.
#  Two separate tests, one significant and one not, differ in significance;
#  they have not been shown to differ from each other. That is the same fallacy
#  as the white cells in Figure 5a and the "C4 amplifies the alcohol response"
#  sentence in the Supp Fig 1 legend.
#
#  The design settles it directly. The cell-means model carries all four
#  groups, so the DIFFERENCE of the two alcohol effects is estimable in one
#  contrast, per cell type, with the mouse as the unit throughout:
#
#      alcohol on variant    :  MUT_WSD_PAL - MUT_WSD_NAL
#      alcohol on wild type  :  CTL_WSD_PAL - CTL_WSD_NAL
#      the extra hit         : (MUT_WSD_PAL - MUT_WSD_NAL) - (CTL_WSD_PAL - CTL_WSD_NAL)
#
#  The third is the claim. Quote the first two side by side only as description.
#
#  EXPECT IT TO BE UNDERPOWERED. A difference of two differences has roughly
#  half the effective power of either one, on 4-5 animals per group. If nothing
#  clears FDR, that is a real answer -- "we cannot demonstrate specificity" --
#  and it is a far better position than asserting specificity from a pair of
#  separate tests and having a reviewer take it apart. The descriptive claim
#  that survives either way is the paired one: across the 24 Cyps, is the
#  alcohol effect more negative on the variant background than on wild type?
#  That is a statement about the gene set, of the same standing as the 18-of-24
#  sign test, and it is printed at the end.
#
#  NOTE ON THE INTERACTION ARM. Katia's memo cut the formal interaction analysis
#  from the main manuscript. This contrast is the gene-arm version of the same
#  idea, so if any of it is quoted, that decision has to be revisited
#  deliberately rather than by accident.
#
#  ALL CELL TYPES. Cyps are largely hepatocyte genes, so in stellate, Kupffer
#  and endothelial populations most will fail detection and drop out. That is
#  informative in its own right -- but "absent because undetected" and "present
#  and unchanged" are different statements, and the n_detected column keeps them
#  apart.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
USE_ALL_CELLS <- TRUE      # TRUE = every cell type; FALSE = the 6 hepatocyte populations
DROP_RUV      <- FALSE     # keep Steve's W_* terms (see fig4_steve_edgeR_qlf.R)
FDR_CUT       <- 0.05
G_MUT_NAL <- "MUT_WSD_NAL"; G_MUT_PAL <- "MUT_WSD_PAL"
G_CTL_NAL <- "CTL_WSD_NAL"; G_CTL_PAL <- "CTL_WSD_PAL"
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
CT_DIR  <- Sys.getenv("FIG4_CT",   file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4",
                      if (USE_ALL_CELLS) "fig4_outs_ALLCELLS" else "fig4_outs"))
CYP     <- Sys.getenv("FIG4_CYP",  file.path(FIGROOT, "Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv"))

suppressPackageStartupMessages({ library(edgeR); library(limma) })
dir.create(OUTS, showWarnings = FALSE, recursive = TRUE)

HEP <- c("Periportal hepatocyte", "Mid-zone hepatocyte", "Pericentral hepatocyte",
         "ISG expressing hepatocyte", "Migratory_IR hepatocyte", "Remodeling hepatocyte")
.ct_name <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
if (!USE_ALL_CELLS) rds <- rds[.ct_name(rds) %in% HEP]
stopifnot("No RDS found in CT_DIR" = length(rds) > 0)

CONTRASTS <- list(
  variant   = sprintf("%s - %s", G_MUT_PAL, G_MUT_NAL),
  wild_type = sprintf("%s - %s", G_CTL_PAL, G_CTL_NAL),
  extra_hit = sprintf("(%s - %s) - (%s - %s)", G_MUT_PAL, G_MUT_NAL, G_CTL_PAL, G_CTL_NAL),
  ## The baseline gap, added 30 Aug. If alcohol moves far fewer genes on the
  ## variant background, one explanation is occlusion: the variant liver on
  ## western diet is already displaced, so alcohol has less room to act. That
  ## is testable -- it predicts a LARGE baseline contrast. If this comes back
  ## small, occlusion is not the explanation and something else is going on.
  baseline  = sprintf("%s - %s", G_MUT_NAL, G_CTL_NAL))

cat(sprintf("\n  edgeR glmQLFTest, %d cell type%s, 3 contrasts each\n",
            length(rds), if (length(rds) == 1) "" else "s"))
for (nm in names(CONTRASTS)) cat(sprintf("    %-10s %s\n", nm, CONTRASTS[[nm]]))
cat("\n")

cyp <- utils::read.csv(CYP, stringsAsFactors = FALSE)
CYPS <- unique(cyp$symbol)
SETS <- stats::setNames(cyp$sets, cyp$symbol)

fit_of <- function(o) {
  fit <- o$fitted_model
  colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
  if (DROP_RUV) {
    d2  <- fit$design[, !grepl("^W_.+", colnames(fit$design)), drop = FALSE]
    fit <- edgeR::glmQLFit(o$dge_list_obj, d2, robust = TRUE)
    colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
  }
  fit
}

LONG <- list(); SUMM <- list()
for (fp in rds) {
  ct <- .ct_name(fp); o <- readRDS(fp)
  fit <- fit_of(o)
  sm  <- o$dge_list_obj$samples
  n <- vapply(c(G_MUT_NAL, G_MUT_PAL, G_CTL_NAL, G_CTL_PAL),
              function(g) sum(sm$cmm_var == g, na.rm = TRUE), numeric(1))
  have <- colnames(fit$design)
  need <- c(G_MUT_NAL, G_MUT_PAL, G_CTL_NAL, G_CTL_PAL)
  if (!all(need %in% have)) {
    cat(sprintf("  %-30s missing group(s): %s -- skipped\n", ct,
                paste(setdiff(need, have), collapse = ", ")))
    next
  }
  de <- c()
  for (nm in names(CONTRASTS)) {
    ctr <- limma::makeContrasts(contrasts = CONTRASTS[[nm]], levels = fit$design)
    tt  <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = ctr), n = Inf)$table
    de[nm] <- sum(tt$FDR < FDR_CUT, na.rm = TRUE)
    i <- match(CYPS, rownames(tt))
    LONG[[length(LONG) + 1]] <- data.frame(
      cell_type = ct, contrast = nm, symbol = CYPS, sets = unname(SETS[CYPS]),
      logFC = tt$logFC[i], PValue = tt$PValue[i], FDR = tt$FDR[i],
      detected = !is.na(i), stringsAsFactors = FALSE)
  }
  SUMM[[length(SUMM) + 1]] <- data.frame(
    cell_type = ct, n_MUT_NAL = n[[1]], n_MUT_PAL = n[[2]],
    n_CTL_NAL = n[[3]], n_CTL_PAL = n[[4]],
    DE_variant = de[["variant"]], DE_wild_type = de[["wild_type"]],
    DE_extra_hit = de[["extra_hit"]], DE_baseline = de[["baseline"]],
    Cyps_detected = sum(CYPS %in% rownames(o$dge_list_obj$counts)),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-30s %d/%d vs %d/%d mice   DE: variant %5d | wild type %5d | extra hit %5d | baseline %5d\n",
              ct, n[[1]], n[[2]], n[[3]], n[[4]],
              de[["variant"]], de[["wild_type"]], de[["extra_hit"]], de[["baseline"]]))
}
stopifnot("No cell type carried all four groups" = length(LONG) > 0)

L <- do.call(rbind, LONG); S <- do.call(rbind, SUMM)
utils::write.csv(L, file.path(OUTS, "cyp_specificity_qlf_long.csv"), row.names = FALSE)
utils::write.csv(S, file.path(OUTS, "cyp_specificity_qlf_celltype_summary.csv"), row.names = FALSE)

## ---- the Cyp picture, averaged over the cell types that detected them ------
agg <- function(nm, col, fun = mean) {
  x <- L[L$contrast == nm & L$detected, ]
  stats::setNames(vapply(CYPS, function(s) {
    v <- x[[col]][x$symbol == s]; v <- v[!is.na(v)]
    if (!length(v)) NA_real_ else fun(v) }, numeric(1)), CYPS)
}
n_det <- stats::setNames(vapply(CYPS, function(s)
  sum(L$detected & L$symbol == s & L$contrast == "variant"), numeric(1)), CYPS)

tab <- data.frame(
  symbol = CYPS, sets = unname(SETS[CYPS]), n_celltypes = unname(n_det),
  logFC_variant   = unname(agg("variant",   "logFC")),
  logFC_wild_type = unname(agg("wild_type", "logFC")),
  logFC_extra_hit = unname(agg("extra_hit", "logFC")),
  FDR_extra_min   = unname(agg("extra_hit", "FDR", function(v) min(v, na.rm = TRUE))),
  stringsAsFactors = FALSE)
tab <- tab[tab$n_celltypes > 0, ]
tab <- tab[order(tab$logFC_variant), ]
utils::write.csv(tab, file.path(OUTS, "cyp_specificity_qlf_summary.csv"), row.names = FALSE)

cat(sprintf("\n  %-11s %-9s %5s %10s %10s %10s %11s\n",
            "gene", "set", "cells", "variant", "wild type", "extra hit", "min FDR ext"))
cat("  --------------------------------------------------------------------------\n")
for (i in seq_len(nrow(tab)))
  cat(sprintf("  %-11s %-9s %5d %+10.2f %+10.2f %+10.2f %11.3g\n",
              tab$symbol[i], tab$sets[i], tab$n_celltypes[i], tab$logFC_variant[i],
              tab$logFC_wild_type[i], tab$logFC_extra_hit[i], tab$FDR_extra_min[i]))

## ---- the paired statement, which is what the figure can claim --------------
d  <- tab$logFC_variant - tab$logFC_wild_type
ok <- !is.na(d)
w  <- suppressWarnings(stats::wilcox.test(tab$logFC_variant[ok], tab$logFC_wild_type[ok],
                                          paired = TRUE))
b  <- stats::binom.test(sum(d[ok] < 0), sum(ok))
cat(sprintf("\n  Across %d Cyps, alcohol effect on variant minus on wild type:\n", sum(ok)))
cat(sprintf("    mean difference          : %+.2f\n", mean(d[ok])))
cat(sprintf("    more negative on variant : %d of %d genes (sign test P = %s)\n",
            sum(d[ok] < 0), sum(ok), format.pval(b$p.value, digits = 2)))
cat(sprintf("    Wilcoxon signed rank     : P = %s\n", format.pval(w$p.value, digits = 2)))
## A minimum FDR taken across 16 cell types is a selection over 16 tests, not a
## test. One gene clearing 0.05 that way is what chance produces. Read this
## number as "nothing unless it is several".
cat(sprintf("    genes with extra-hit FDR < %.2f in any cell type : %d  (min over %d cell types -- not a test)\n",
            FDR_CUT, sum(tab$FDR_extra_min < FDR_CUT, na.rm = TRUE),
            length(unique(L$cell_type))))
cat(sprintf("\n  Baseline gap (variant vs wild type, no alcohol), summed over cell types: %d genes\n",
            sum(S$DE_baseline)))
cat("    A large baseline gap alongside a small variant-arm response supports\n")
cat("    occlusion -- already displaced, less room for alcohol to act. A small\n")
cat("    one rules that explanation out.\n")

cat("\n  HOW TO SAY IT\n")
cat("  If the extra-hit contrast clears FDR for some genes, specificity is\n")
cat("  demonstrated for those genes and can be stated outright.\n")
cat("  If it does not, the honest sentence is descriptive and paired -- \"across\n")
cat("  the 24 Cyps the alcohol effect was more negative on the variant\n")
cat("  background\" with the Wilcoxon P -- and NOT \"the Cyps were affected only\n")
cat("  on the variant background\", which the data would not support.\n")
cat("  Genes are not independent either, so that Wilcoxon is of the same\n")
cat("  standing as the 18-of-24 sign test: descriptive support, not proof.\n")
cat(sprintf("\n  written -> %s\n", file.path(OUTS, "cyp_specificity_qlf_summary.csv")))
