# =============================================================================
#  Animal-level differential expression, and a direct comparison against the
#  stacked population-level model.
#
#  WHY THIS EXISTS
#  The consolidated script fits limma on pseudo-samples stacked as
#  mouse x hepatocyte population (~0 + group + zone), blocking on mouse through
#  duplicateCorrelation. That is the correct limma idiom for repeated measures
#  and it is NOT naive pseudoreplication -- mouse is in the model. But it is
#  under-corrected here, for three reasons that compound:
#
#   1. duplicateCorrelation estimates ONE consensus correlation for every gene
#      in the matrix. That consensus came out at 0.210, dragged down by the tens
#      of thousands of genes with no mouse-level structure. For a gene like
#      Cyp4a10, where every population in a mouse moves together, the true intra-mouse
#      correlation is far higher. Underestimating it makes limma treat the six
#      populations as closer to independent than they are, so the standard error is
#      too small and the P value too small.
#   2. Each pseudo-sample is aggregated over thousands of nuclei, so technical
#      noise is negligible and nearly all residual variance IS the between-mouse
#      variance that step 1 is under-weighting.
#   3. Residual d.f. is set by mice x populations, not by mice.
#
#  The experiment has a handful of animals per group. No test on those animals
#  can support an FDR of 1e-13; that exponent is the six populations being counted.
#  This script computes the statistic at the level the experiment actually
#  randomised -- one value per mouse, the mean across its six hepatocyte populations,
#  which is exactly what the panel e boxplots plot -- and puts the two side by
#  side so the size of the gap is visible rather than argued about.
#
#  WHAT TO REPORT
#  The animal-level number. A boxplot showing one dot per mouse annotated with a
#  P value computed on mice x populations is an internal inconsistency a reviewer will
#  find, and Sanyal has already raised replication units once on this manuscript.
#  The population-level model stays as the pathway-level engine and as a sensitivity
#  analysis; it is not the number that goes next to the dots.
# =============================================================================

ROOT     <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT  <- Sys.getenv("FIGROOT",  file.path(ROOT, "writing skeleton/figure making"))
CT_DIR   <- Sys.getenv("FIG4_CT",  file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUTS     <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
CYP      <- Sys.getenv("FIG4_CYP",  file.path(FIGROOT, "Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv"))

suppressPackageStartupMessages({ library(limma); library(nlme) })

QUESTION  <- "Q3"
BOX_REF   <- "MUT_WSD_NAL"
BOX_ALT   <- "MUT_WSD_PAL"
EXPR_SLOT <- "log2cpm_rb"

## TERMINOLOGY. The R variables below are still called zone / HEP_ZONES because
## that is how the RDS filenames and the upstream Stouffer Rmd name them. In
## anything a reader sees they are hepatocyte POPULATIONS -- annotated snRNA-seq
## clusters, not microdissected zones. Three of the six carry zonal names
## (periportal, mid-zone, pericentral) because that is how they annotate; the
## other three (interferon-stimulated, migratory/IR, remodelling) are states.
## No zonation analysis was done and no spatial information was used, so
## "six zones" in a label or a legend overclaims.
HEP_ZONES <- c("Periportal hepatocyte","Mid-zone hepatocyte","Pericentral hepatocyte",
               "ISG expressing hepatocyte","Migratory_IR hepatocyte","Remodeling hepatocyte")

## ---- 1. stack the zones, exactly as the consolidated script does ------------
.ct_name <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
hep <- rds[.ct_name(rds) %in% HEP_ZONES]
stopifnot("No hepatocyte-population RDS found in CT_DIR" = length(hep) > 0)

mats <- list(); metas <- list()
for (fp in hep) {
  z <- .ct_name(fp); o <- readRDS(fp)
  xm <- as.matrix(o$dge_list_obj[[EXPR_SLOT]]); storage.mode(xm) <- "double"
  colnames(xm) <- gsub("_.*", "", colnames(xm))
  sm  <- o$dge_list_obj$samples
  s2g <- stats::setNames(as.character(sm$cmm_var), gsub("_.*", "", as.character(sm$sample)))
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

## ---- 2. collapse to ONE value per mouse ------------------------------------
mouse_grp <- stats::setNames(meta$group[!duplicated(meta$mouse)], meta$mouse[!duplicated(meta$mouse)])
permouse  <- sapply(names(mouse_grp), function(mo)
  rowMeans(mat[, rownames(meta)[meta$mouse == mo], drop = FALSE]))
grp_of    <- as.character(mouse_grp[colnames(permouse)])

keep <- grp_of %in% c(BOX_REF, BOX_ALT)
X    <- permouse[, keep, drop = FALSE]
g    <- factor(grp_of[keep], levels = c(BOX_REF, BOX_ALT))
n1   <- sum(g == BOX_REF); n2 <- sum(g == BOX_ALT)

cat(sprintf("\n  pseudo-samples in the stacked model : %d  (%d mice x %d populations)\n",
            ncol(mat), length(unique(meta$mouse)), length(unique(meta$zone))))
cat(sprintf("  animals in this contrast            : %d  (%d %s vs %d %s)\n",
            n1 + n2, n1, BOX_REF, n2, BOX_ALT))

## The exact permutation floor. With these group sizes there are choose(n, n2)
## distinct labellings, so no rank-based test on these animals can return a P
## below 1/that, however large the effect. It is the ceiling on what the design
## can evidence, and it is worth printing next to any 1e-13.
nperm <- choose(n1 + n2, n2)
cat(sprintf("  exact permutation floor             : 1/%d = %.2g one-sided, %.2g two-sided\n",
            nperm, 1 / nperm, min(1, 2 / nperm)))

## ---- 3. animal-level limma --------------------------------------------------
## Still moderated -- variance is borrowed across genes, which is legitimate and
## is what limma is for. What changes is the unit: d.f. now comes from animals.
des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
fit <- lmFit(X, des)
ctr <- makeContrasts(contrasts = paste(BOX_ALT, "-", BOX_REF), levels = des)
fit <- eBayes(contrasts.fit(fit, ctr), trend = TRUE, robust = TRUE)
an  <- topTable(fit, number = Inf)
an  <- data.frame(SYMBOL = rownames(an), logFC_animal = an$logFC,
                  P_animal = an$P.Value, FDR_animal = an$adj.P.Val,
                  stringsAsFactors = FALSE)

## An unmoderated Welch t as a no-shrinkage floor, so the effect of borrowing
## variance across genes is also visible.
welch <- apply(X, 1, function(v) {
  tt <- try(stats::t.test(v[g == BOX_ALT], v[g == BOX_REF]), silent = TRUE)
  if (inherits(tt, "try-error")) NA_real_ else tt$p.value })
an$P_welch <- unname(welch[an$SYMBOL])

utils::write.csv(an, file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION)),
                 row.names = FALSE)

## ---- 4. side by side, for the Cyps -----------------------------------------
st <- utils::read.csv(file.path(OUTS, sprintf("allHep_limma_genes_%s.csv", QUESTION)),
                      stringsAsFactors = FALSE)
st <- data.frame(SYMBOL = st$SYMBOL, logFC_zone = st$logFC, FDR_zone = st$adj.P.Val,
                 stringsAsFactors = FALSE)

cyp <- utils::read.csv(CYP, stringsAsFactors = FALSE)
cmp <- merge(merge(cyp[, c("symbol","sets")], st, by.x = "symbol", by.y = "SYMBOL"),
             an, by.x = "symbol", by.y = "SYMBOL")
cmp <- cmp[order(cmp$logFC_zone), ]

## Intra-mouse correlation for THESE genes, against the single consensus limma
## used for all genes. If these run much higher than the consensus, that is the
## under-correction, measured rather than asserted.
icc <- vapply(cmp$symbol, function(s) {
  if (!s %in% rownames(mat)) return(NA_real_)
  v  <- mat[s, ]; mo <- meta$mouse; gg <- meta$group
  v  <- v - stats::ave(v, gg)                      # remove the group effect
  mm <- stats::ave(v, mo)                          # mouse means
  vb <- stats::var(mm[!duplicated(mo)]); vw <- stats::var(v - mm)
  if (is.na(vb) || is.na(vw) || (vb + vw) == 0) NA_real_ else vb / (vb + vw)
}, numeric(1))
cmp$ICC_mouse <- unname(icc)

## A proper random-intercept mixed model on the full stacked data, per gene.
## This is the arbiter: it uses every population, estimates the mouse variance from
## THAT gene rather than from a global consensus, and is what duplicateCorrelation
## is approximating. If P_lme tracks the animal-level column and not the zone
## column, the collapse to per-mouse means is sound and the consensus is the
## problem. (In a balanced design the two are near-identical; they separate a
## little if some mice are missing populations.)
mm <- meta[meta$group %in% c(BOX_REF, BOX_ALT), ]
mm$group <- factor(mm$group, levels = c(BOX_REF, BOX_ALT))
p_lme <- vapply(cmp$symbol, function(s) {
  if (!s %in% rownames(mat)) return(NA_real_)
  dd <- data.frame(v = mat[s, rownames(mm)], group = mm$group,
                   zone = factor(mm$zone), mouse = factor(mm$mouse))
  out <- try(nlme::lme(v ~ group + zone, random = ~ 1 | mouse, data = dd, method = "REML"),
             silent = TRUE)
  if (inherits(out, "try-error")) NA_real_
  else summary(out)$tTable[paste0("group", BOX_ALT), "p-value"]
}, numeric(1))
cmp$P_lme <- unname(p_lme)

cat("\n  gene         set        logFC   FDR (population model)   FDR (animal)    P mixed    P Welch   ICC\n")
cat(  "  -------------------------------------------------------------------------------------------\n")
for (i in seq_len(nrow(cmp))) if (!is.na(cmp$logFC_zone[i]))
  cat(sprintf("  %-11s %-9s %+6.2f   %14.3g   %12.3g   %9.3g   %8.3g   %4.2f\n",
              cmp$symbol[i], cmp$sets[i], cmp$logFC_zone[i], cmp$FDR_zone[i],
              cmp$FDR_animal[i], cmp$P_lme[i], cmp$P_welch[i], cmp$ICC_mouse[i]))

sz <- sum(cmp$FDR_zone < 0.05, na.rm = TRUE); sa <- sum(cmp$FDR_animal < 0.05, na.rm = TRUE)
cat(sprintf("\n  Cyps at FDR < 0.05 : %d by the population model, %d at the animal level\n", sz, sa))
cat(sprintf("  median ICC for these Cyps : %.2f   (limma used one consensus of 0.210 for every gene)\n",
            stats::median(cmp$ICC_mouse, na.rm = TRUE)))
cat("\n  READ IT LIKE THIS: if P mixed sits with FDR (animal) and far from FDR (pop.\n")
cat("  model), the animal-level number is the honest one and the consensus\n")
cat("  correlation is what inflated the population model. If the ICC column runs well\n")
cat("  above 0.210, that is the mechanism, measured.\n")
cat(sprintf("\n  written -> %s\n", file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION))))
