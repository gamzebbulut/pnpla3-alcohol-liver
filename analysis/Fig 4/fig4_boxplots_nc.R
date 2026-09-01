# =============================================================================
#  Figure 4 panel e -- per-mouse expression of the significant Cyps
#
#  The last step of the chain. Panel b puts arachidonic and linoleic acid
#  metabolism concordantly down; panel c puts PPAR signalling most decreased in
#  the gene arm; panel d shows those three sets are largely the same Cyps. This
#  panel is the underlying data for that claim -- one point per mouse, so a
#  reader can see the effect is not carried by one animal.
#
#  THE FDR IN THE STRIP IS THE ANIMAL-LEVEL ONE, from
#  fig4_animal_level_limma.R. The dots are mice, so the number beside them must
#  be computed on mice. The population-stacked model in the consolidated script fits
#  mouse x population pseudo-samples and returns FDRs down to 1e-13, which no
#  experiment with this many animals can support; annotating a per-mouse boxplot
#  with it is an inconsistency a reviewer will find. Set FDR_SOURCE to "zone" to
#  see the old numbers.
#
#  SELF-CONTAINED and cheap. It rebuilds the per-mouse 6-hepatocyte expression
#  matrix from the population RDS files, then reads the FDRs that
#  fig4_consolidated_pathview_boxplots.R already wrote to
#  fig4_outs/allHep_limma_genes_<Q>.csv. It does NOT refit limma and does NOT
#  touch pathview, so it runs in seconds and can be re-run freely while you
#  tune the panel. Run the consolidated script once first, for that CSV.
#
#  REPLACES the per-pathway boxplot grids in the consolidated script, which
#  render every significant member of every pathway at 2.3 inches per gene --
#  fine as a working readout, far too large for a 180 mm figure.
#
#  Colour follows Figure 2d: grey = no alcohol, blue = alcohol. NOT the red/blue
#  of panels b-d, which mean increased/decreased there. One figure must not use
#  one colour pair for two different things.
# =============================================================================

ROOT     <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT  <- Sys.getenv("FIGROOT",  file.path(ROOT, "writing skeleton/figure making"))
CT_DIR   <- Sys.getenv("FIG4_CT",  file.path(ROOT, "Data/pseudobulk_snRNA-seq_PNPLA3_alc"))
OUTS     <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
CYP      <- Sys.getenv("FIG4_CYP",  file.path(FIGROOT, "Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv"))
OUT      <- Sys.getenv("FIG4_OUT",  file.path(FIGROOT, "Fig 4/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

QUESTION  <- "Q3"                     # = contrast C1, alcohol on the variant
BOX_REF   <- "MUT_WSD_NAL"            # variant, no alcohol
BOX_ALT   <- "MUT_WSD_PAL"            # variant + alcohol
EXPR_SLOT <- "log2cpm_rb"             # RUV-corrected log2 CPM
FDR_CUT   <- 0.05
NCOL      <- 3
FDR_SOURCE <- "animal"   # "animal" (correct unit) or "pop" (the old inflated one)

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

## ---- 1. per-mouse expression, averaged across the six hepatocyte populations ------
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

mouse_grp <- stats::setNames(meta$group[!duplicated(meta$mouse)], meta$mouse[!duplicated(meta$mouse)])
permouse  <- sapply(names(mouse_grp), function(mo)
  rowMeans(mat[, rownames(meta)[meta$mouse == mo], drop = FALSE]))
grp_of    <- as.character(mouse_grp[colnames(permouse)])
cat(sprintf("  %d genes x %d mice (%d %s, %d %s)\n", nrow(permouse), ncol(permouse),
            sum(grp_of == BOX_REF), BOX_REF, sum(grp_of == BOX_ALT), BOX_ALT))

## ---- 2. which genes ---------------------------------------------------------
## The significant Cyps from panel d, in the same order, so the two panels read
## together. Override GENES here to show a different set.
if (FDR_SOURCE == "animal") {
  af <- file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION))
  if (!file.exists(af)) stop("Run fig4_animal_level_limma.R first, or set FDR_SOURCE <- \"zone\".")
  de <- utils::read.csv(af, stringsAsFactors = FALSE)
  qlook <- stats::setNames(de$FDR_animal, de$SYMBOL)
} else {
  de <- utils::read.csv(file.path(OUTS, sprintf("allHep_limma_genes_%s.csv", QUESTION)),
                        stringsAsFactors = FALSE)
  qlook <- stats::setNames(de$adj.P.Val, de$SYMBOL)
}

cyp  <- utils::read.csv(CYP, stringsAsFactors = FALSE)
cyp  <- cyp[!is.na(cyp$FDR) & cyp$FDR < FDR_CUT, ]
GENES <- cyp$symbol[order(cyp$logFC)]
GENES <- intersect(GENES, rownames(permouse))
cat(sprintf("  plotting %d genes: %s\n", length(GENES), paste(GENES, collapse = ", ")))

## ---- 3. long frame ----------------------------------------------------------
keep <- grp_of %in% c(BOX_REF, BOX_ALT)
long <- do.call(rbind, lapply(GENES, function(s) data.frame(
  gene = s, mouse = colnames(permouse)[keep], v = permouse[s, keep],
  grp = factor(grp_of[keep], levels = c(BOX_REF, BOX_ALT)), stringsAsFactors = FALSE)))
long$gene <- factor(long$gene, levels = GENES)

## Facet strips carry the gene name AND its FDR, so no asterisk column and no
## separate significance legend is needed. Plotmath, via fmt_p_math from
## nc_style.R -- fmt_p alone would print a literal "10^-10" in the strip.
lab <- stats::setNames(
  sprintf("italic('%s')~~'FDR'~%s", GENES,
          vapply(GENES, function(s) {
            q <- qlook[[s]]
            if (is.na(q)) "'NA'" else if (q >= 1e-4) sprintf("'%s'", fmt_p(q)) else fmt_p_math(q)
          }, character(1))),
  GENES)

p <- ggplot(long, aes(grp, v, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.25, colour = NC_INK) +
  geom_jitter(width = 0.12, height = 0, size = 0.7, shape = 21,
              colour = NC_INK, stroke = 0.2, fill = "white") +
  facet_wrap(~ gene, ncol = NCOL, scales = "free_y",
             labeller = ggplot2::as_labeller(lab, default = ggplot2::label_parsed)) +
  scale_fill_manual(values = stats::setNames(c(NC_GREY, NC_BLUE), c(BOX_REF, BOX_ALT)),
                    guide = "none") +
  scale_x_discrete(labels = c("No alcohol", "Alcohol")) +
  ## Both groups carry the variant. Left implicit, "no alcohol vs alcohol" reads
  ## as though the groups were pooled across genotypes; naming the background on
  ## the axis makes it unambiguous that this is contrast C1.
  labs(x = "PNPLA3 I148M background",
       y = expression("Log"[2]~"CPM, six hepatocyte populations")) +
  theme_nc(grid_y = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = NC_MIN, hjust = 0,
                                  margin = margin(b = 1)),
        panel.spacing = grid::unit(4, "pt"),
        axis.text.x = element_text(size = NC_MIN),
        axis.title.x = element_text(size = NC_MIN, colour = NC_INK2)) +
  nc_tag("e")

save_nc(p, file.path(OUT, "Fig4e_cyp_boxplots"), width_mm = 104, height_mm = 96)

utils::write.csv(long, file.path(OUT, "SourceData_Fig4e_boxplots.csv"), row.names = FALSE)
