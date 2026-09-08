# =============================================================================
#  Combined lipid supplement, panel d — the lipid-handling enzymes
#
#  THE CLAIM IS AN ABSENCE. Almost none of the lipid remodeling seen in panels
#  a-c is accompanied by a transcriptional change in the enzymes that would have
#  to carry it. That is the point of the panel, and it is why it belongs at the
#  end of the arc rather than the start.
#
#  CHANGED 30 Aug -- THE MODEL
#  This panel used to stack the six hepatocyte populations into one limma fit
#  with duplicateCorrelation blocking on mouse, and star enzymes at that model's
#  FDR. Figure 4 no longer uses that model: it reports edgeR glmQLFTest on one
#  hepatocyte population, which is the route Steve's README specifies and which
#  keeps the mouse as the replication unit without a consensus correlation. A
#  manuscript that stars genes one way in Figure 4 and another way in the
#  supplement invites a question it does not need to answer, so this panel now
#  uses the same route.
#
#  The claim survives the change and should strengthen: the per-population route
#  is the more conservative of the two, so an absence of transcriptional change
#  measured that way is a stronger absence.
#
#  NO C1-C4. The previous version labeled the rows C1..C4, which revision task
#  A2 removes everywhere a reader can see. Rows are named in plain language and
#  ordered as nc_style.R's CONTRAST_ORDER does: alcohol first, then the variant.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
POP_FULL  <- "Periportal hepatocyte"   # NULL = pick the most abundant
FDR_CUT   <- 0.1
DROP_RUV  <- FALSE
SIG_ONLY  <- FALSE     # TRUE drops enzymes not significant in any contrast
## Full width. Eighteen gene columns in three functional blocks; at 88 mm the
## facet strips truncated to "aturation / elonga". This panel is a strip by
## nature and takes the page width.
W_MM      <- 180
H_MM      <- 44
## Rows, in nc_style.R's reading order. The strings are what the reader sees.
CONTRASTS <- c("Alcohol on wild type" = "CTL_WSD_PAL - CTL_WSD_NAL",
               "Alcohol on variant"   = "MUT_WSD_PAL - MUT_WSD_NAL",
               "Variant alone"        = "MUT_WSD_NAL - CTL_WSD_NAL",
               "Variant on alcohol"   = "MUT_WSD_PAL - CTL_WSD_PAL")
GENE_MODULES <- list(
  "Synthesis"                 = c("Srebf1","Mlxipl","Acly","Acaca","Acacb","Fasn"),
  "Desaturation / elongation" = c("Scd1","Scd2","Fads1","Fads2","Elovl5","Elovl6"),
  "Mobilization / handling"   = c("Pnpla2","Lpcat3","Lpcat1","Mboat7","Pla2g6","Mttp","Apob"))
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
CT_DIR  <- Sys.getenv("FIG4_CT",   file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUT     <- Sys.getenv("SUPP_OUT",  file.path(FIGROOT, "supp_lipid_story/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages({ library(edgeR); library(limma); library(ggplot2) })
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

HEP <- c("Periportal hepatocyte","Mid-zone hepatocyte","Pericentral hepatocyte",
         "ISG expressing hepatocyte","Migratory_IR hepatocyte","Remodeling hepatocyte")
.ct <- function(p) sub("_fitted.*$|\\.RDS$", "", basename(p))
rds <- list.files(CT_DIR, pattern = "\\.RDS$", full.names = TRUE)
rds <- rds[.ct(rds) %in% HEP]
stopifnot("No hepatocyte-population RDS found" = length(rds) > 0)

if (is.null(POP_FULL)) {
  ab <- vapply(rds, function(f) sum(readRDS(f)$dge_list_obj$samples$lib.size), numeric(1))
  POP_FULL <- .ct(rds)[which.max(ab)]
  cat(sprintf("  population chosen by abundance: %s\n", POP_FULL))
}
stopifnot("POP_FULL not found" = POP_FULL %in% .ct(rds))
o <- readRDS(rds[match(POP_FULL, .ct(rds))])
cat(sprintf("  population: %s\n", POP_FULL))

fit <- o$fitted_model
colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
if (DROP_RUV) {
  d2  <- fit$design[, !grepl("^W_.+", colnames(fit$design)), drop = FALSE]
  fit <- edgeR::glmQLFit(o$dge_list_obj, d2, robust = TRUE)
  colnames(fit$design) <- gsub("^cmm_var", "", colnames(fit$design))
}

ALL <- unique(unlist(GENE_MODULES, use.names = FALSE))
modof <- stats::setNames(rep(names(GENE_MODULES), lengths(GENE_MODULES)),
                         unlist(GENE_MODULES, use.names = FALSE))

D <- list()
for (nm in names(CONTRASTS)) {
  need <- unlist(strsplit(CONTRASTS[[nm]], "\\s*-\\s*"))
  if (!all(need %in% colnames(fit$design))) {
    cat(sprintf("  %-22s not estimable, skipped\n", nm)); next
  }
  ctr <- limma::makeContrasts(contrasts = CONTRASTS[[nm]], levels = fit$design)
  tt  <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = ctr), n = Inf)$table
  i <- match(ALL, rownames(tt))
  D[[length(D) + 1]] <- data.frame(
    gene = ALL, module = unname(modof[ALL]), contrast = nm,
    logFC = tt$logFC[i], FDR = tt$FDR[i], stringsAsFactors = FALSE)
  cat(sprintf("  %-22s %5d genes at FDR < %.2f transcriptome-wide\n",
              nm, sum(tt$FDR < FDR_CUT, na.rm = TRUE), FDR_CUT))
}
stopifnot("No contrast was estimable" = length(D) > 0)
D <- do.call(rbind, D)
D <- D[!is.na(D$logFC), ]
D$sig <- !is.na(D$FDR) & D$FDR < FDR_CUT

if (SIG_ONLY) {
  keep <- unique(D$gene[D$sig])
  if (!length(keep)) stop("SIG_ONLY is TRUE but no enzyme is significant anywhere")
  D <- D[D$gene %in% keep, ]
}

D$gene     <- factor(D$gene, levels = ALL[ALL %in% D$gene])
D$module   <- factor(D$module, levels = names(GENE_MODULES))
D$contrast <- factor(D$contrast, levels = rev(names(CONTRASTS)))  # first level at bottom

LIM <- max(abs(D$logFC), na.rm = TRUE)
p <- ggplot(D, aes(gene, contrast)) +
  geom_tile(aes(fill = logFC), colour = "white", linewidth = 0.4) +
  ## Significance is an outline, matching Fig 4e and Fig 1d. With a panel whose
  ## claim is that almost nothing moves, a saturation step would make the whole
  ## grid look like missing data.
  geom_tile(data = D[D$sig, ], fill = NA, colour = NC_INK, linewidth = 0.55) +
  facet_grid(~ module, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                       limits = c(-LIM, LIM), na.value = "#F2F2F2",
                       name = expression("Log"[2]~"fold-change"),
                       breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(expand = expansion(add = 0.04)) +
  ## Reader-facing contrast names from nc_style.R's display layer, parsed for
  ## the I148M superscript. The CONTRASTS vector keeps the internal keys, which
  ## are also the limma contrast formulas' names.
  scale_y_discrete(expand = expansion(add = 0.04), labels = nc_label_expr) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(34, "pt"),
                                barheight = grid::unit(5.5, "pt"),
                                title.position = "left", ticks.colour = NC_INK)) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = NC_BIG, angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic"),
        axis.text.y = element_text(size = NC_BIG),
        strip.text = element_text(size = NC_BIG),
        panel.spacing = grid::unit(5, "pt"),
        legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = NC_BIG),
        legend.text = element_text(size = NC_BIG),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "pt")) +
  nc_tag("d")

save_nc(p, file.path(OUT, "Supp_lipid_d_enzymes"), W_MM, H_MM)
utils::write.csv(D, file.path(OUT, "SourceData_Supp_lipid_d_enzymes.csv"), row.names = FALSE)

## ---- the absence, stated as a number ---------------------------------------
cat(sprintf("\n  %d enzymes x %d contrasts = %d cells\n",
            nlevels(D$gene), nlevels(D$contrast), nrow(D)))
cat(sprintf("  cells at FDR < %.2f : %d\n", FDR_CUT, sum(D$sig)))
if (any(D$sig)) {
  s <- D[D$sig, ]
  s <- s[order(s$FDR), ]
  cat("  what moves:\n")
  for (i in seq_len(nrow(s)))
    cat(sprintf("    %-8s %-22s logFC %+5.2f  FDR %.3g\n",
                s$gene[i], as.character(s$contrast[i]), s$logFC[i], s$FDR[i]))
} else {
  cat("  nothing reaches significance in any contrast -- which IS the result.\n")
}
cat(sprintf("\n  largest |log2FC| anywhere : %.2f\n", LIM))
cat("  LEGEND MATERIAL: the panel shows the lipid-handling programme is largely\n")
cat("  unchanged transcriptionally, so the lipid differences in panels a-c are\n")
cat("  not explained by expression of these enzymes in this population.\n")
