# =============================================================================
#  Is it the whole PPARa programme, or only the Cyp4a genes?  + two QC checks
#
#  THE HYPOTHESIS THIS TESTS (Gamze, 30 Aug)
#  The variant retains lipid. Fatty acids are PPARa ligands, so on the variant
#  background PPARa stays ligand-activated and its target genes -- the Cyp4a
#  omega-hydroxylases among them -- are less suppressed by alcohol. On the wild
#  type there is less retained lipid to keep PPARa on, so the decrease is
#  bigger. That is a real mechanism and it makes a sharp prediction:
#
#      the blunting should show up across the PPARa TARGET PROGRAMME,
#      not only in the four Cyp4a genes.
#
#  If Acox1, Ehhadh, Hmgcs2, Cpt1a, Acot1/2, Fgf21 and the rest move the same
#  way, the hypothesis has a programme behind it. If only Cyp4a does, the
#  hypothesis is not supported and the Cyp4a result is something narrower.
#  Either answer is worth having before the Discussion commits to a mechanism.
#
#  FDR UNIVERSE -- READ THIS BEFORE COMPARING TO ISABEL'S NUMBERS
#  Isabel's May 2026 CYP report computed FDR over the CYP genes only (748 tests
#  for Q3), which is why it reports 119 hits at FDR < 0.05 where a
#  transcriptome-wide analysis reports a handful. Neither is wrong; they answer
#  different questions and MUST NOT be mixed in one paper -- that is the trap
#  that turned a q of 0.034 into 0.19 in the amino-acid panel. This script
#  reports BOTH, in labeled columns, so the choice is explicit.
#
#  SECTION 2 is unrelated QC: why does alcohol move ~40x more genes on the wild
#  type than on the variant, in every cell type? Before that goes in the paper
#  it is worth knowing whether the variant mice are simply more variable, and
#  whether the RUV factors track genotype.
# =============================================================================

USE_ALL_CELLS <- FALSE     # FALSE = the six hepatocyte populations
FDR_CUT       <- 0.05
G_MUT_NAL <- "MUT_WSD_NAL"; G_MUT_PAL <- "MUT_WSD_PAL"
G_CTL_NAL <- "CTL_WSD_NAL"; G_CTL_PAL <- "CTL_WSD_PAL"

## PPARa targets, canonical mouse set, split so the question can be answered.
## Cyp4a is deliberately separate: it is the family the result rests on, and a
## programme-level claim must not be circular.
PPARA_CORE <- c("Acox1","Ehhadh","Acadm","Acadvl","Cpt1a","Cpt2","Slc25a20",
                "Hmgcs2","Fgf21","Acot1","Acot2","Plin5","Pdk4","Angptl4",
                "Cd36","Lpl","Scp2","Crat","Aldh3a2","Hadha","Hadhb","Etfdh")
PPARA_CYP4A <- c("Cyp4a10","Cyp4a14","Cyp4a31","Cyp4a32","Cyp4a12a","Cyp4a12b")
SETS <- list(`PPARa targets (non-Cyp)` = PPARA_CORE, `Cyp4a family` = PPARA_CYP4A)

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
CT_DIR  <- Sys.getenv("FIG4_CT",   file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4",
                      if (USE_ALL_CELLS) "fig4_outs_ALLCELLS" else "fig4_outs"))

suppressPackageStartupMessages({ library(edgeR); library(limma) })
dir.create(OUTS, showWarnings = FALSE, recursive = TRUE)

HEP <- c("Periportal hepatocyte","Mid-zone hepatocyte","Pericentral hepatocyte",
         "ISG expressing hepatocyte","Migratory_IR hepatocyte","Remodeling hepatocyte")
.ct <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
if (!USE_ALL_CELLS) rds <- rds[.ct(rds) %in% HEP]
stopifnot("No RDS found" = length(rds) > 0)

CTRS <- c(variant   = sprintf("%s - %s", G_MUT_PAL, G_MUT_NAL),
          wild_type = sprintf("%s - %s", G_CTL_PAL, G_CTL_NAL),
          extra_hit = sprintf("(%s - %s) - (%s - %s)",
                              G_MUT_PAL, G_MUT_NAL, G_CTL_PAL, G_CTL_NAL))
GENES <- unlist(SETS, use.names = FALSE)
setof <- stats::setNames(rep(names(SETS), lengths(SETS)), GENES)

L <- list(); QC <- list()
for (fp in rds) {
  ct <- .ct(fp); o <- readRDS(fp)
  fit <- o$fitted_model
  colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
  if (!all(c(G_MUT_NAL, G_MUT_PAL, G_CTL_NAL, G_CTL_PAL) %in% colnames(fit$design))) next

  for (nm in names(CTRS)) {
    ctr <- limma::makeContrasts(contrasts = CTRS[[nm]], levels = fit$design)
    tt  <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = ctr), n = Inf)$table
    i   <- match(GENES, rownames(tt))
    ## FDR within the gene set as well as transcriptome-wide, both labeled.
    q_set <- rep(NA_real_, length(GENES))
    ok <- !is.na(i); q_set[ok] <- stats::p.adjust(tt$PValue[i[ok]], "BH")
    L[[length(L) + 1]] <- data.frame(
      cell_type = ct, contrast = nm, symbol = GENES, set = unname(setof[GENES]),
      logFC = tt$logFC[i], PValue = tt$PValue[i],
      FDR_transcriptome = tt$FDR[i], FDR_within_set = q_set,
      stringsAsFactors = FALSE)
  }

  ## ---- QC: how variable is each group, and do the RUV factors track genotype?
  X  <- o$dge_list_obj$log2cpm_rb
  colnames(X) <- gsub("_cluster.+", "", colnames(X))
  sm <- o$dge_list_obj$samples
  g  <- stats::setNames(as.character(sm$cmm_var), as.character(sm$sample))
  g  <- g[colnames(X)]
  spread <- vapply(c(G_MUT_NAL, G_MUT_PAL, G_CTL_NAL, G_CTL_PAL), function(k) {
    cols <- names(g)[!is.na(g) & g == k]
    if (length(cols) < 2) return(NA_real_)
    ## mean pairwise correlation distance between mice in that group: a direct,
    ## model-free measure of within-group heterogeneity
    mean(1 - stats::cor(X[, cols, drop = FALSE])[lower.tri(diag(length(cols)))])
  }, numeric(1))
  ## FIXED 30 Aug: this read the genotype via g[rownames(sm)], but the samples
  ## table is row-numbered 1..n, not named by sample, so every lookup was NA and
  ## the column printed -Inf. Take the genotype from cmm_var directly, and
  ## restrict to the four western-diet groups so the chow arm cannot drive it.
  W  <- fit$design[, grepl("^W_", colnames(fit$design)), drop = FALSE]
  cv <- as.character(sm$cmm_var)
  keep <- cv %in% c(G_MUT_NAL, G_MUT_PAL, G_CTL_NAL, G_CTL_PAL)
  geno <- as.numeric(grepl("^MUT", cv[keep]))
  wcor <- NA_real_
  if (ncol(W) && sum(keep) > 2 && stats::sd(geno) > 0) {
    cc <- suppressWarnings(stats::cor(W[keep, , drop = FALSE], geno))
    if (any(is.finite(cc))) wcor <- max(abs(cc[is.finite(cc)]))
  }
  QC[[length(QC) + 1]] <- data.frame(
    cell_type = ct, spread_MUT_NAL = spread[[1]], spread_MUT_PAL = spread[[2]],
    spread_CTL_NAL = spread[[3]], spread_CTL_PAL = spread[[4]],
    max_abs_cor_W_genotype = wcor, stringsAsFactors = FALSE)
}
stopifnot("No cell type carried all four groups" = length(L) > 0)

D <- do.call(rbind, L); Q <- do.call(rbind, QC)
utils::write.csv(D, file.path(OUTS, "ppara_programme_long.csv"), row.names = FALSE)
utils::write.csv(Q, file.path(OUTS, "ppara_qc_group_spread.csv"), row.names = FALSE)

## ---- 1. the programme question ---------------------------------------------
cat("\n  ===== Is the blunting a PPARa programme, or only Cyp4a? =====\n")
cat(sprintf("  %-24s %5s %10s %10s %10s %8s\n",
            "gene set", "genes", "variant", "wild type", "difference", "n<0 var"))
cat("  --------------------------------------------------------------------------\n")
for (s in names(SETS)) {
  gs <- SETS[[s]]
  mv <- vapply(gs, function(x) mean(D$logFC[D$symbol == x & D$contrast == "variant"],   na.rm = TRUE), numeric(1))
  mw <- vapply(gs, function(x) mean(D$logFC[D$symbol == x & D$contrast == "wild_type"], na.rm = TRUE), numeric(1))
  ok <- !is.na(mv) & !is.na(mw)
  cat(sprintf("  %-24s %5d %+10.2f %+10.2f %+10.2f %5d/%d\n", s, sum(ok),
              mean(mv[ok]), mean(mw[ok]), mean(mv[ok] - mw[ok]),
              sum(mv[ok] < mw[ok]), sum(ok)))
  if (sum(ok) >= 5) {
    w <- suppressWarnings(stats::wilcox.test(mv[ok], mw[ok], paired = TRUE))
    cat(sprintf("      paired Wilcoxon across genes: P = %s\n",
                format.pval(w$p.value, digits = 2)))
  }
}
cat("\n  READ IT LIKE THIS: the hypothesis predicts BOTH rows lean the same way\n")
cat("  (difference positive = less decrease on the variant). If only the Cyp4a\n")
cat("  row does, the PPARa-ligand mechanism is not what is driving it.\n")

## ---- 2. the QC question ----------------------------------------------------
cat("\n  ===== Why does the wild-type arm move more genes? =====\n")
cat(sprintf("  %-26s %9s %9s %9s %9s %10s\n", "cell type", "MUT-NAL", "MUT-PAL",
            "CTL-NAL", "CTL-PAL", "cor(W,geno)"))
for (i in seq_len(nrow(Q)))
  cat(sprintf("  %-26s %9.3f %9.3f %9.3f %9.3f %10.2f\n", Q$cell_type[i],
              Q$spread_MUT_NAL[i], Q$spread_MUT_PAL[i], Q$spread_CTL_NAL[i],
              Q$spread_CTL_PAL[i], Q$max_abs_cor_W_genotype[i]))
cat("\n  Columns are mean within-group correlation distance between mice: bigger\n")
cat("  means the animals in that group resemble each other less. edgeR pools one\n")
cat("  dispersion per gene across all groups, so if the two MUT groups are much\n")
cat("  more variable than the two CTL groups, the variant contrast is tested\n")
cat("  against a variance the CTL groups largely set, and vice versa -- which\n")
cat("  would be a technical explanation for the 40-fold gap rather than biology.\n")
cat("  cor(W, genotype) near 1 would mean the RUV factors absorb genotype, which\n")
cat("  would blunt the variant contrast for a second, separate reason.\n")
