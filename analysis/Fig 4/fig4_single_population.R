# =============================================================================
#  One hepatocyte population, one model, one set of numbers
#     -> pathview map AND the boxplots beneath it, from the same test
#
#  WHY ONE POPULATION
#  Pooling six populations means the map's colours and the boxplot's dots come
#  from an average, and every reader question ("which populations? weighted
#  how?") lands on the legend. Running one population removes all of it: the
#  map is coloured by that population's log2 fold-change, the boxplot shows that
#  population's per-mouse values, and the number beside the dots is that
#  population's FDR. Nothing is pooled, nothing is stacked, and the replication
#  unit is the mouse throughout because Steve's per-cell-type DGEList already
#  has one column per mouse.
#
#  HOW THE POPULATION IS CHOSEN -- THIS MATTERS
#  By ABUNDANCE, not by effect size. Picking the population with the biggest
#  effect and then showing it as the result is selection on the outcome, and it
#  is the same class of error as reading blank cells as specificity. Abundance
#  is decided by the data before any contrast is fitted, so it cannot be
#  accused of that. The ranking is printed on every run and belongs in the
#  legend: "the most abundant hepatocyte population".
#  POP_OVERRIDE exists for when there is a stated biological reason to choose
#  differently -- put the reason in the legend if you use it.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
QUESTION     <- "Q3"
BOX_REF      <- "MUT_WSD_NAL"; BOX_ALT <- "MUT_WSD_PAL"
POP_OVERRIDE <- NULL      # e.g. "Periportal hepatocyte"; NULL = most abundant
RUN_TAG      <- "20260830_periportal"   # the folder these maps are filed under
FDR_GENE     <- 0.10
## Pathways to render on top of the joint list. Keyed by KEGG id.
EXTRA_PATHWAYS <- c(mmu03320 = "PPAR signaling pathway")
LIMIT        <- 3         # colour limit for GENES
## Compounds need their own, smaller limit. The measured metabolite effects run
## to about |1.6|, so on a +-3 scale arachidonate at -0.49 renders almost white
## and the map reads as though no metabolite were measured at all. pathview
## draws a separate legend for compounds, so a different limit is not a
## sleight of hand -- the scale is printed beside the map.
CPD_LIMIT    <- 1.5
## DEFAULT IS "kegg_default" (changed 30 Aug), for one specific reason: in the
## manuscript scheme the MIDPOINT is white, so a metabolite that was measured and
## did not change is indistinguishable from a metabolite that was never measured.
## KEGG's own scheme puts grey at the midpoint, so every measured compound is
## visibly filled and only the unmeasured ones stay white. On a map where two of
## roughly thirty compounds are on our platform, that difference is the whole
## readability of the metabolite layer.
##
## The cost is real and should be stated in the legend: KEGG's gene ramp is
## green-to-red, which is the pairing colour-vision reviewers object to, and it
## contradicts the house rule that red means increased and blue decreased. These
## are supplementary reference maps, which is where breaking that rule costs
## least. If you want the house palette back, set "manuscript" -- but then also
## change the compound midpoint away from white, or the metabolites vanish again.
##
## "manuscript"  = blue/red for both, matching Figures 1-3, colourblind-safe.
## "kegg_default"= KEGG's own look: genes green->red, compounds blue->yellow.
##   That is what the CURRENT Figure 4 and Supp Fig 5 legends already describe,
##   so choosing it makes those legends correct with no rewrite -- but it breaks
##   the house rule that red means increased and blue decreased everywhere a
##   direction is shown, and red/green is the one pairing colour-vision
##   reviewers reliably object to.
COLOR_SCHEME <- "kegg_default"
CEX_NATIVE   <- 0.25      # KEGG-style raster
CEX_VECTOR   <- 0.55      # Graphviz PDF: 0.25 is why the labels looked lost
N_BOX        <- 8         # genes boxplotted per pathway
WRITE_MAPS   <- TRUE
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
CT_DIR  <- Sys.getenv("FIG4_CT",   file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
MET_CSV <- Sys.getenv("FIG4_MET",  file.path(FIGROOT, "Fig 1",
                      sprintf("FC_%s_metabolites.csv", QUESTION)))
PMAX_CSV <- Sys.getenv("FIG4_PMAX", file.path(ROOT, "GB_metabolomics_20260628",
                      "Q3_StoufferPMAX_pooled", "modalities_pmax_gene_plus_metDec12.csv"))

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(ggplot2)
  library(pathview); library(org.Mm.eg.db); library(AnnotationDbi) })
source(file.path(FIGROOT, "nc_style.R"))

HEP <- c("Periportal hepatocyte","Mid-zone hepatocyte","Pericentral hepatocyte",
         "ISG expressing hepatocyte","Migratory_IR hepatocyte","Remodeling hepatocyte")
.ct <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
rds <- rds[.ct(rds) %in% HEP]
stopifnot("No hepatocyte-population RDS found" = length(rds) > 0)

## ---- 1. abundance ranking, before any contrast is fitted -------------------
## Total pseudobulk library size is the available proxy for how much of the
## tissue a population represents; nuclei counts are not stored in the DGEList.
ab <- vapply(rds, function(fp) sum(readRDS(fp)$dge_list_obj$samples$lib.size), numeric(1))
names(ab) <- .ct(rds); ab <- sort(ab, decreasing = TRUE)
cat("\n  hepatocyte populations by total pseudobulk library size\n")
for (i in seq_along(ab))
  cat(sprintf("    %d. %-28s %12.0f  (%4.1f%%)\n", i, names(ab)[i], ab[i],
              100 * ab[i] / sum(ab)))
POP <- if (!is.null(POP_OVERRIDE)) POP_OVERRIDE else names(ab)[1]
if (!is.null(POP_OVERRIDE))
  cat(sprintf("\n  POP_OVERRIDE in use: %s -- give the reason in the legend\n", POP))
cat(sprintf("\n  population used: %s\n\n", POP))
## The run tag names the output folder, so a tag that disagrees with the
## population silently mislabels every map filed under it.
if (!grepl(tolower(gsub("[^A-Za-z]", "", strsplit(POP, " ")[[1]][1])), tolower(RUN_TAG)))
  warning(sprintf("RUN_TAG '%s' does not name the population '%s' -- maps will be filed under a tag that does not describe them", RUN_TAG, POP), call. = FALSE)
stopifnot("POP not among the hepatocyte populations" = POP %in% names(ab))

o <- readRDS(rds[match(POP, .ct(rds))])

## ---- 2. edgeR QLF for that population, Steve's route -----------------------
fit <- o$fitted_model
colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
stopifnot("contrast not estimable in this population" =
            all(c(BOX_REF, BOX_ALT) %in% colnames(fit$design)))
ctr <- limma::makeContrasts(contrasts = paste(BOX_ALT, "-", BOX_REF), levels = fit$design)
de  <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = ctr), n = Inf)$table
de  <- data.frame(SYMBOL = rownames(de), logFC = de$logFC, PValue = de$PValue,
                  FDR = de$FDR, stringsAsFactors = FALSE)
POPTAG <- gsub("[^A-Za-z0-9]+", "-", POP)
utils::write.csv(de, file.path(OUTS, sprintf("singlepop_%s_%s_genes.csv", POPTAG, QUESTION)),
                 row.names = FALSE)
cat(sprintf("  %d genes tested, %d at FDR < 0.05, %d at FDR < %.2f\n",
            nrow(de), sum(de$FDR < 0.05, na.rm = TRUE),
            sum(de$FDR < FDR_GENE, na.rm = TRUE), FDR_GENE))

## ---- 3. per-mouse expression for the boxplots, same population -------------
X <- o$dge_list_obj$log2cpm_rb
colnames(X) <- gsub("_cluster.+", "", colnames(X))
sm <- o$dge_list_obj$samples
grp <- stats::setNames(as.character(sm$cmm_var), as.character(sm$sample))[colnames(X)]
keep <- !is.na(grp) & grp %in% c(BOX_REF, BOX_ALT)
X <- X[, keep, drop = FALSE]; grp <- factor(grp[keep], levels = c(BOX_REF, BOX_ALT))
cat(sprintf("  boxplot samples: %d %s vs %d %s (one dot per mouse)\n\n",
            sum(grp == BOX_REF), BOX_REF, sum(grp == BOX_ALT), BOX_ALT))

## ---- 4. pathways: the joint list, plus the gene-arm extras ------------------
## The joint integration is METABOLITE-LIMITED by construction -- a pathway can
## only qualify if it carries measured metabolites -- so the strongest pathway in
## the gene arm never enters it. PPAR signalling is Z = -6.67, the most decreased
## pathway in the whole gene arm, and it is absent from the 13. Adding it back as
## an explicit extra is not cherry-picking as long as the figure and the index
## say which pathways were joint hits and which were added; that is the `extra`
## column below.
## WHICH PATHWAYS GET A MAP -- and why it is 13, not 14 (noted 31 Aug)
## The manuscript defines the joint list by the two-threshold rule (gene
## FDR < 0.10 AND metabolite P < 0.05), which gives 14 and is what Figure 4
## reports. THIS filter is stricter: FDR_Pmax < 0.10, a BH correction on the
## larger of the two arms' P values. It keeps 13 of those 14.
##
## The one it drops is mmu04964, proximal tubule bicarbonate reclamation --
## discordant, and a kidney pathway. We decided not to render it, and the
## Supplementary Fig. 3 legend says "13 of the 14" accordingly.
##
## TO RENDER ALL 14 INSTEAD: replace the FDR_Pmax filter below with the joint
## list itself --
##     J   <- utils::read.csv(file.path(FIGROOT, "Fig 4/C1_joint_hits.csv"))
##     sig <- sig[sig$KEGG_id %in% J$KEGG_id, ]
## -- and change the Supp Fig 3 legend back to "the 14 joint pathways".
sig <- utils::read.csv(PMAX_CSV, stringsAsFactors = FALSE)
sig <- sig[!is.na(sig$FDR_Pmax) & sig$FDR_Pmax < 0.10 & sig$concordance != "One-missing", ]
sig <- sig[order(sig$FDR_Pmax), ]
sig$extra <- FALSE
if (length(EXTRA_PATHWAYS)) {
  add <- names(EXTRA_PATHWAYS)[!names(EXTRA_PATHWAYS) %in% sig$KEGG_id]
  if (length(add)) {
    ex <- data.frame(KEGG_id = add, KEGG_title = unname(EXTRA_PATHWAYS[add]),
                     FDR_Pmax = NA_real_, concordance = "gene arm only",
                     extra = TRUE, stringsAsFactors = FALSE)
    keep <- intersect(names(sig), names(ex))
    sig  <- rbind(sig[, keep], ex[, keep])
    cat(sprintf("  extra pathways added (not joint hits): %s\n",
                paste(sprintf("%s %s", add, unname(EXTRA_PATHWAYS[add])), collapse = "; ")))
  }
}
kegg_sets <- pathview::korg  # touch the package so a broken install fails early
slug <- function(x) { x <- gsub("\\s*-\\s*Mus musculus.*$", "", x)
  gsub("^-|-$", "", gsub("[^A-Za-z0-9]+", "-", substr(x, 1, 60))) }

## ---- 5. maps, coloured by THIS population ----------------------------------
## Everything lands under Fig 4/pathview_outs/<run tag>/ so the maps are filed by
## which cells coloured them, and the KEGG downloads are shared by every run
## rather than duplicated per folder. See pathview_outs/README.md.
PV_ROOT <- Sys.getenv("FIG4_PVROOT", file.path(FIGROOT, "Fig 4/pathview_outs"))
EXPORT <- file.path(PV_ROOT, RUN_TAG)
D_PNG <- file.path(EXPORT, "native_png"); D_PDF <- file.path(EXPORT, "vector_pdf")
D_CAC <- file.path(PV_ROOT, "kegg_cache"); D_BOX <- file.path(EXPORT, "boxplots")
for (d in c(D_PNG, D_PDF, D_CAC, D_BOX)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

mp <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(de$SYMBOL),
                            keytype = "SYMBOL", columns = "ENTREZID")
mp <- mp[!is.na(mp$ENTREZID) & !duplicated(mp$SYMBOL), ]
ge <- merge(de[, c("SYMBOL", "logFC")], mp, by = "SYMBOL")
gene.data <- stats::setNames(ge$logFC, ge$ENTREZID)

PV <- if (identical(COLOR_SCHEME, "kegg_default")) {
  list(low  = list(gene = "green", cpd = "blue"),
       mid  = list(gene = "gray",  cpd = "gray"),
       high = list(gene = "red",   cpd = "yellow"))
} else {
  list(low  = list(gene = NC_BLUE_TINT, cpd = NC_BLUE_TINT),
       mid  = list(gene = "#FFFFFF",    cpd = "#FFFFFF"),
       high = list(gene = NC_RED_TINT,  cpd = NC_RED_TINT))
}
cat(sprintf("  colour scheme: %s\n", COLOR_SCHEME))

cpd.data <- NULL
if (file.exists(MET_CSV)) {
  m <- utils::read.csv(MET_CSV, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("KEGG", "log2_FC") %in% names(m))) {
    k <- trimws(as.character(m$KEGG))
    ok <- !is.na(k) & nzchar(k) & k != "NA" & !is.na(m$log2_FC)
    cc <- tapply(m$log2_FC[ok], sub("[;, ].*$", "", k[ok]), mean)
    cpd.data <- stats::setNames(as.numeric(cc), names(cc))
    cat(sprintf("  cpd.data: %d KEGG compounds, |log2FC| up to %.2f\n",
                length(cpd.data), max(abs(cpd.data), na.rm = TRUE)))
  }
}
if (is.null(cpd.data))
  cat("  cpd.data: NONE -- maps will carry one legend, not two\n")

## The index is the map's provenance: which population coloured it, under which
## run tag and colour scheme, and -- the column a reader will ask for -- whether
## the pathway was a joint metabolite+gene hit or added by hand from the gene arm.
IDX <- list()
if (WRITE_MAPS) {
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(D_CAC)
  for (i in seq_len(nrow(sig))) {
    kid <- sig$KEGG_id[i]; pid <- sub("^mmu", "", kid)
    base <- sprintf("Fig4_%s_%s_%s_%s", kid, slug(sig$KEGG_title[i]), POPTAG, QUESTION)
    for (nat in c(TRUE, FALSE)) {
      sfx <- if (nat) "native" else "graphviz"
      args <- list(gene.data = gene.data, pathway.id = pid, species = "mmu",
                   out.suffix = sfx, kegg.native = nat, same.layer = !nat,
                   cex = if (nat) CEX_NATIVE else CEX_VECTOR,
                   limit = list(gene = LIMIT, cpd = CPD_LIMIT),
                   low = PV$low, mid = PV$mid, high = PV$high)
      if (!is.null(cpd.data)) args$cpd.data <- cpd.data
      ok <- tryCatch({ do.call(pathview, args); TRUE },
                     error = function(e) { message("  ", kid, " ", sfx, ": ",
                                                   conditionMessage(e)); FALSE })
      if (ok) {
        src <- sprintf("mmu%s.%s.%s", pid, sfx, if (nat) "png" else "pdf")
        dst <- file.path(if (nat) D_PNG else D_PDF,
                         paste0(base, if (nat) ".png" else ".pdf"))
        if (file.exists(src)) file.rename(src, dst)
      }
    }
    cat(sprintf("  map %-10s %-46s %s\n", kid, substr(slug(sig$KEGG_title[i]), 1, 46),
                if (isTRUE(sig$extra[i])) "EXTRA (gene arm only)" else ""))
    IDX[[length(IDX) + 1]] <- data.frame(
      KEGG_id = kid, KEGG_title = sig$KEGG_title[i], extra = isTRUE(sig$extra[i]),
      FDR_Pmax = sig$FDR_Pmax[i], population = POP, run = RUN_TAG,
      color_scheme = COLOR_SCHEME, stringsAsFactors = FALSE)
  }
  setwd(old)
  if (length(IDX)) {
    idx <- do.call(rbind, IDX)
    idx <- idx[order(idx$extra, idx$FDR_Pmax), ]
    utils::write.csv(idx, file.path(EXPORT, "INDEX_pathviews.csv"), row.names = FALSE)
    cat(sprintf("\n  index: %d maps -> %s\n", nrow(idx),
                file.path(basename(EXPORT), "INDEX_pathviews.csv")))
    cat(sprintf("         %d joint hits, %d added from the gene arm (extra = TRUE)\n",
                sum(!idx$extra), sum(idx$extra)))
  }
}

## ---- 6. boxplots, SAME population, SAME test -------------------------------
## The number under each gene is this population's FDR, printed exactly rather
## than as a star: the manuscript reports exact values (revision task J8), and
## with one population at n = 5 vs 4 a star would hide how close most of these
## are to the line.
gset <- tryCatch(pathview::pathway.gene.list, error = function(e) NULL)  # may be absent
for (i in seq_len(nrow(sig))) {
  kid <- sig$KEGG_id[i]
  syms <- ge$SYMBOL[ge$SYMBOL %in% rownames(X)]
  d_i  <- de[de$SYMBOL %in% syms, ]
  d_i  <- d_i[order(d_i$FDR), ]
  ## Without a pathway->gene map in this session, fall back to the strongest
  ## genes overall; supply KEGG_GENES externally to restrict to the pathway.
  pick <- head(d_i$SYMBOL, N_BOX)
  if (!length(pick)) next
  df <- do.call(rbind, lapply(pick, function(s) data.frame(
    gene = s, value = as.numeric(X[s, ]), group = grp, stringsAsFactors = FALSE)))
  df$gene <- factor(df$gene, levels = pick)
  lab <- stats::setNames(sprintf("%s\nFDR = %s", pick,
           fmt_p(de$FDR[match(pick, de$SYMBOL)])), pick)
  p <- ggplot(df, aes(group, value, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, colour = NC_INK) +
    geom_jitter(width = 0.12, size = 0.8, shape = 16, colour = NC_INK) +
    facet_wrap(~ gene, scales = "free_y", nrow = 2,
               labeller = labeller(gene = lab)) +
    scale_fill_manual(values = stats::setNames(c(NC_GREY, NC_BLUE),
                                               c(BOX_REF, BOX_ALT)), guide = "none") +
    scale_x_discrete(labels = c("No alcohol", "Alcohol")) +
    labs(x = NULL, y = expression("Log"[2]~"CPM")) +
    theme_nc() +
    theme(strip.text = element_text(size = NC_MIN, face = "italic"),
          axis.text.x = element_text(size = NC_MIN))
  save_nc(p, file.path(D_BOX, sprintf("Boxplots_%s_%s_%s", kid, POPTAG, QUESTION)),
          180, 80)
}

cat(sprintf("\n  written to %s\n", EXPORT))
cat("    native_png/  raster, KEGG's own resolution -- cannot be increased\n")
cat("    vector_pdf/  Graphviz layout, LIVE TEXT: restyle in Illustrator,\n")
cat("                 do not retype -- retyped labels cannot be re-checked\n")
cat("    boxplots/    same population, same test as the map above them\n")
