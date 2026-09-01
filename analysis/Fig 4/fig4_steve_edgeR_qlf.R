# =============================================================================
#  Figure 4, gene arm — the sanctioned route, for comparison
#
#  WHY THIS SCRIPT EXISTS
#  Steve's README (READMESept26.html, "Differential Expression Analysis")
#  specifies edgeR glmQLFTest on the DGEGLM already fitted inside each RDS, one
#  cell type at a time. Nothing in this repository used it: the gene arm has
#  been limma throughout, and the consolidated version stacked the six
#  hepatocyte populations into one model with duplicateCorrelation blocking on
#  mouse. That stacking is ours, not his -- in his workflow each cell type's
#  DGEList already has ONE COLUMN PER MOUSE, so his alcohol contrast is 5 vs 4
#  animals and the replication-unit question never arises.
#
#  This script runs his route unmodified so the two can be compared. It reports
#  nothing on its own; it exists to answer "do the sanctioned per-population
#  tests agree with the animal-level limma we propose to report?"
#
#  WHAT IT DOES NOT DO, DELIBERATELY
#  It does not combine the six populations into one P value. Stouffer or Fisher
#  across six populations measured in the SAME nine animals treats correlated
#  tests as independent evidence -- the identical error, in a different place,
#  that made the stacked model return an FDR of 1e-13. The honest summary of six
#  per-population tests is how many of them a gene clears, which is descriptive,
#  and that is what the k_sig column is. If a combined P is wanted later it needs
#  a correlation-aware method, not a plain Stouffer.
#
#  Each population's FDR is transcriptome-wide (BH over every gene in that
#  population's DGEList), the same universe as the animal-level column. Do not
#  recompute q over the 24 Cyps alone -- see the Fig 5b note in HANDOFF.md §8.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
QUESTION  <- "Q3"                 # alcohol on the variant background
BOX_REF   <- "MUT_WSD_NAL"
BOX_ALT   <- "MUT_WSD_PAL"
DROP_RUV  <- FALSE   # FALSE = keep the W_* factors of unwanted variation, which
                     # is Steve's first example and the design his fitted_model
                     # already carries. TRUE re-fits without them; his README
                     # notes this usually yields fewer DEGs. Both are printed
                     # either way -- this switch only picks which one is primary.
FDR_CUT   <- 0.05
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
CT_DIR  <- Sys.getenv("FIG4_CT",   file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
CYP     <- Sys.getenv("FIG4_CYP",  file.path(FIGROOT, "Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv"))

suppressPackageStartupMessages({ library(edgeR); library(limma) })
dir.create(OUTS, showWarnings = FALSE, recursive = TRUE)

## TERMINOLOGY: these are hepatocyte POPULATIONS (annotated snRNA-seq clusters),
## not zones. Three carry zonal names because that is how they annotate; no
## zonation analysis was done and no spatial information was used.
HEP <- c("Periportal hepatocyte", "Mid-zone hepatocyte", "Pericentral hepatocyte",
         "ISG expressing hepatocyte", "Migratory_IR hepatocyte", "Remodeling hepatocyte")

.ct_name <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
rds <- rds[.ct_name(rds) %in% HEP]
stopifnot("No hepatocyte-population RDS found in CT_DIR" = length(rds) > 0)

cat(sprintf("\n  edgeR glmQLFTest, %s - %s, one test per hepatocyte population\n",
            BOX_ALT, BOX_REF))
cat(sprintf("  populations found : %d of %d\n\n", length(rds), length(HEP)))

## ---- one QLF test per population, exactly as the README does it ------------
qlf_one <- function(obj, drop_ruv) {
  fit <- obj$fitted_model
  colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
  if (drop_ruv) {
    d2 <- fit$design[, !grepl("^W_.+", colnames(fit$design)), drop = FALSE]
    fit <- edgeR::glmQLFit(obj$dge_list_obj, d2, robust = TRUE)
    colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
  }
  if (!all(c(BOX_REF, BOX_ALT) %in% colnames(fit$design))) return(NULL)
  ctr <- limma::makeContrasts(contrasts = paste(BOX_ALT, "-", BOX_REF),
                              levels = fit$design)
  tt  <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = ctr), n = Inf)$table
  data.frame(SYMBOL = rownames(tt), logFC = tt$logFC, PValue = tt$PValue,
             FDR = tt$FDR, stringsAsFactors = FALSE)
}

RES <- list(); RES_NORUV <- list(); NMICE <- list()
for (fp in rds) {
  z <- .ct_name(fp); o <- readRDS(fp)
  sm <- o$dge_list_obj$samples
  n_ref <- sum(sm$cmm_var == BOX_REF, na.rm = TRUE)
  n_alt <- sum(sm$cmm_var == BOX_ALT, na.rm = TRUE)
  NMICE[[z]] <- c(ref = n_ref, alt = n_alt)

  a <- qlf_one(o, FALSE); b <- qlf_one(o, TRUE)
  if (is.null(a)) { cat(sprintf("  %-26s contrast not estimable, skipped\n", z)); next }
  RES[[z]] <- a; RES_NORUV[[z]] <- b
  cat(sprintf("  %-26s %d vs %d mice   %5d DE at FDR<%.2f  (%5s without the W_ terms)\n",
              z, n_ref, n_alt, sum(a$FDR < FDR_CUT, na.rm = TRUE), FDR_CUT,
              if (is.null(b)) "NA" else sum(b$FDR < FDR_CUT, na.rm = TRUE)))
}
stopifnot("No population produced an estimable contrast" = length(RES) > 0)

PRIMARY <- if (DROP_RUV) RES_NORUV else RES
PRIMARY <- PRIMARY[!vapply(PRIMARY, is.null, logical(1))]

## ---- the Cyp panel, population by population -------------------------------
cyp <- utils::read.csv(CYP, stringsAsFactors = FALSE)
cyp <- cyp[!is.na(cyp$logFC), c("symbol", "sets", "logFC", "FDR")]
names(cyp) <- c("symbol", "sets", "logFC_stacked", "FDR_stacked")

grab <- function(sym, what) vapply(PRIMARY, function(d) {
  i <- match(sym, d$SYMBOL); if (is.na(i)) NA_real_ else d[[what]][i] }, numeric(1))

tab <- cyp
for (z in names(PRIMARY)) {
  tab[[paste0("logFC.", z)]] <- vapply(tab$symbol,
                                       function(s) unname(grab(s, "logFC")[z]), numeric(1))
  tab[[paste0("FDR.",   z)]] <- vapply(tab$symbol,
                                       function(s) unname(grab(s, "FDR")[z]),   numeric(1))
}
fdr_cols <- grep("^FDR\\.", names(tab), value = TRUE)
fc_cols  <- grep("^logFC\\.", names(tab), value = TRUE)
tab$k_sig    <- rowSums(tab[, fdr_cols, drop = FALSE] < FDR_CUT, na.rm = TRUE)
tab$k_tested <- rowSums(!is.na(tab[, fdr_cols, drop = FALSE]))
tab$logFC_median_pop <- apply(tab[, fc_cols, drop = FALSE], 1, stats::median, na.rm = TRUE)
tab$FDR_min_pop      <- suppressWarnings(apply(tab[, fdr_cols, drop = FALSE], 1, min, na.rm = TRUE))
tab$FDR_min_pop[!is.finite(tab$FDR_min_pop)] <- NA_real_

## ---- the animal-level limma we propose to report ---------------------------
af <- file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION))
if (file.exists(af)) {
  an <- utils::read.csv(af, stringsAsFactors = FALSE)
  i <- match(tab$symbol, an$SYMBOL)
  tab$logFC_animal <- an$logFC_animal[i]
  tab$FDR_animal   <- an$FDR_animal[i]
} else {
  tab$logFC_animal <- NA_real_; tab$FDR_animal <- NA_real_
  cat("\n  NOTE: run fig4_animal_level_limma.R first to fill the animal-level columns\n")
}

tab <- tab[order(tab$logFC_stacked), ]
utils::write.csv(tab, file.path(OUTS, sprintf("steve_edgeR_qlf_Cyps_%s.csv", QUESTION)),
                 row.names = FALSE)
for (z in names(PRIMARY))
  utils::write.csv(PRIMARY[[z]],
    file.path(OUTS, sprintf("steve_edgeR_qlf_%s_%s.csv", QUESTION, gsub("[^A-Za-z0-9]+", "_", z))),
    row.names = FALSE)

## ---- the comparison, which is the whole point ------------------------------
cat(sprintf("\n  %-11s %-9s %7s %12s %12s %10s\n",
            "gene", "set", "logFC", "FDR stacked", "FDR animal", "pops sig"))
cat("  ---------------------------------------------------------------------\n")
for (i in seq_len(nrow(tab)))
  cat(sprintf("  %-11s %-9s %+7.2f %12.3g %12.3g %6d / %d\n",
              tab$symbol[i], tab$sets[i], tab$logFC_stacked[i], tab$FDR_stacked[i],
              tab$FDR_animal[i], tab$k_sig[i], tab$k_tested[i]))

ns <- sum(tab$FDR_stacked < FDR_CUT, na.rm = TRUE)
na_ <- sum(tab$FDR_animal  < FDR_CUT, na.rm = TRUE)
nq <- sum(tab$k_sig > 0)
cat(sprintf("\n  Cyps at FDR < %.2f\n", FDR_CUT))
cat(sprintf("    stacked limma (populations as samples) : %d of %d\n", ns, nrow(tab)))
cat(sprintf("    animal-level limma (per-mouse means)   : %d of %d\n", na_, nrow(tab)))
cat(sprintf("    edgeR QLF, significant in >= 1 population : %d of %d\n", nq, nrow(tab)))
cat(sprintf("    edgeR QLF, significant in >= 3 populations: %d of %d\n",
            sum(tab$k_sig >= 3), nrow(tab)))

cat("\n  READ IT LIKE THIS: the sanctioned per-population route already uses the\n")
cat("  mouse as its unit, so if its k_sig column tracks the animal-level FDR and\n")
cat("  not the stacked FDR, the two independent routes agree and the Methods\n")
cat("  sentence becomes 'both approaches gave the same answer'. If a gene is\n")
cat("  significant in five of six populations but fails the animal-level test,\n")
cat("  that is the moderated-vs-unmoderated difference, not a contradiction:\n")
cat("  each population test borrows variance across genes within that population.\n")
cat(sprintf("\n  written -> %s\n", file.path(OUTS, sprintf("steve_edgeR_qlf_Cyps_%s.csv", QUESTION))))
