# =============================================================================
#  Figure 5 panel b — amino-acid metabolites under alcohol on the variant
#
#  Replaces fig5_aminoacid_C1_slim_boxplots.R.
#
#  The heatmap in panel a says an amino-acid block moves together on the variant
#  background. This panel is the underlying data for six of those metabolites,
#  one point per mouse, chosen to span the pathways that block is made of:
#  glycine, threonine and sarcosine from glycine/serine/threonine metabolism,
#  cysteine from cysteine and methionine metabolism, phenylalanine from
#  Phe/Tyr/Trp biosynthesis, valine from the 2-oxocarboxylic acid map.
#  Sarcosine is the glycine N-methylation product, so it reports on one-carbon
#  flux directly.
#
#  CHANGED from the previous version:
#    - COLOUR. It used blue for no-alcohol and red for alcohol, i.e. the
#      increased/decreased palette spent on group identity. Now grey and blue,
#      the same pair Figure 2d and Figure 4e use for exactly this contrast.
#    - The q calculation had a real bug. It computed the metabolome-wide p
#      vector, then for each panel appended that metabolite's own p to the
#      vector and re-ran qvalue, taking the last element. But that p is ALREADY
#      in the vector, so every panel's q was estimated from 471 values with one
#      duplicated. q is now computed once over the 470 and looked up.
#    - Exact P values through fmt_p (task J8) instead of a "p<0.001" floor.
#    - patchwork dropped in favour of facet_wrap, so this panel is built the same
#      way as Figure 4e and cannot drift from it.
#    - 180 mm, 5-7 pt, vector SVG.
#    - Each strip shows q as well as P, and BOTH ARE READ FROM THE FIGURE 1
#      TABLE rather than recomputed here. This matters more than it sounds.
#
#  WHY q IS NOT RECOMPUTED ON THIS MATRIX
#  An earlier version estimated q over the metabolites in met_mat_contrast.csv,
#  the 119 KEGG-mapped compounds the integration uses. That produced q = 0.17 to
#  0.24 for these six and made it look as though none of them cleared the
#  manuscript's q < 0.10 bar. They all do. The artefact is pure denominator:
#  glycine is rank 15 of 470 on the full platform, and 13 of the 14 metabolites
#  ahead of it are lipids. Restricting to the KEGG-mapped subset deletes those
#  lipids, glycine becomes rank 1 of 119, and BH returns p x n = 0.19.
#  P is identical either way -- same data, same Welch test. Only the multiple-
#  testing universe moved. Figure 1's universe is the manuscript's universe, so
#  this panel uses it, and the two figures now agree by construction. The script
#  recomputes P anyway and stops if it disagrees with the table.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("FIG5_MET",   file.path(ROOT, "GB_metabolomics_20260628/Q3_StoufferPMAX_pooled"))
FC      <- Sys.getenv("FIG5_FC",    file.path(FIGROOT, "Fig 1/FC_Q3_metabolites.csv"))
OUT     <- Sys.getenv("FIG5_OUT",   file.path(FIGROOT, "Fig 5/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

REF <- "MUT_WSD_NAL"   # variant, no alcohol
TRT <- "MUT_WSD_PAL"   # variant + alcohol

## KEGG compound id -> display name. Edit here to change the panel set.
PANELS <- data.frame(
  cpd  = c("C00037", "C00188", "C00213", "C00097", "C00079", "C00183"),
  name = c("Glycine", "Threonine", "Sarcosine", "Cysteine", "Phenylalanine", "Valine"),
  stringsAsFactors = FALSE)

mat <- utils::read.csv(file.path(IN, "met_mat_contrast.csv"), check.names = FALSE,
                       stringsAsFactors = FALSE)
grp <- utils::read.csv(file.path(IN, "sample_groups_met.csv"), stringsAsFactors = FALSE)
g   <- stats::setNames(grp$group, grp$sample)

samp <- names(g)[g %in% c(REF, TRT)]
samp <- intersect(samp, colnames(mat))
gg   <- factor(unname(g[samp]), levels = c(REF, TRT))
cat(sprintf("  %d metabolites x %d mice (%d without alcohol, %d with)\n",
            nrow(mat), length(samp), sum(gg == REF), sum(gg == TRT)))

## ---- P and q from the Figure 1 table, with a consistency check --------------
X <- as.matrix(mat[, samp, drop = FALSE]); rownames(X) <- mat$cpd_id

fc <- utils::read.csv(FC, check.names = FALSE, stringsAsFactors = FALSE)
key <- trimws(as.character(fc$INDIVIDUAL.NOTATION))
P_TAB <- stats::setNames(fc$p, key)
Q_TAB <- stats::setNames(fc$qvalue, key)
cat(sprintf("  P and q read from %s (%d metabolites, the Figure 1 universe)\n",
            basename(FC), nrow(fc)))

## ---- long frame for the six panels -----------------------------------------
PANELS <- PANELS[PANELS$cpd %in% rownames(X), ]
miss <- setdiff(PANELS$name, names(P_TAB))
if (length(miss)) stop("Not found in the Figure 1 table by name: ", paste(miss, collapse = ", "))
PANELS$P <- unname(P_TAB[PANELS$name])
PANELS$q <- unname(Q_TAB[PANELS$name])

## Guard: the values plotted come from met_mat_contrast.csv, the P and q from
## the Figure 1 table. They should be the same underlying measurements. Verify
## rather than assume -- a silent mismatch here would put one experiment's P
## beside another experiment's boxes.
chk <- vapply(seq_len(nrow(PANELS)), function(i) {
  v <- X[PANELS$cpd[i], samp]
  stats::t.test(v[gg == TRT], v[gg == REF], var.equal = FALSE)$p.value
}, numeric(1))
dev <- abs(log10(chk) - log10(PANELS$P))
if (any(dev > 0.02)) stop("Recomputed P disagrees with the Figure 1 table for: ",
                          paste(PANELS$name[dev > 0.02], collapse = ", "))
cat(sprintf("  consistency check passed: recomputed P matches the table for all %d (max %.1f%% deviation)\n",
            nrow(PANELS), 100 * max(abs(chk / PANELS$P - 1))))
long <- do.call(rbind, lapply(seq_len(nrow(PANELS)), function(i) data.frame(
  name = PANELS$name[i], sample = samp, val = X[PANELS$cpd[i], samp],
  grp = gg, stringsAsFactors = FALSE)))

## Strip carries the metabolite, its exact P AND its q. Showing P alone would be
## misleading here: every one of these six is nominally significant while none
## clears q < 0.10. The inferential claim in this figure is the PATHWAY-level
## result in panel a; these six illustrate what that block is made of. Printing
## both numbers keeps that honest without a sentence of hedging in the legend.
lab <- stats::setNames(
  vapply(seq_len(nrow(PANELS)), function(i)
    sprintf("atop(%s, italic(P)~'='~'%s'*','~italic(q)~'='~'%s')",
            PANELS$name[i], fmt_p(PANELS$P[i]), fmt_p(PANELS$q[i])),
    character(1)), PANELS$name)
long$name <- factor(long$name, levels = PANELS$name)

p <- ggplot(long, aes(grp, val, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.25, colour = NC_INK) +
  geom_jitter(width = 0.12, height = 0, size = 0.8, shape = 21,
              colour = NC_INK, stroke = 0.2, fill = "white") +
  facet_wrap(~ name, nrow = 1, scales = "free_y",
             labeller = ggplot2::as_labeller(lab, default = ggplot2::label_parsed)) +
  scale_fill_manual(values = stats::setNames(c(NC_GREY, NC_BLUE), c(REF, TRT)), guide = "none") +
  scale_x_discrete(labels = c("No alcohol", "Alcohol")) +
  labs(x = NULL, y = "Relative abundance") +
  theme_nc(grid_y = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = NC_MIN, hjust = 0.5, lineheight = 0.9,
                                  margin = margin(b = 1)),
        panel.spacing = grid::unit(4, "pt"),
        axis.text.x = element_text(size = NC_MIN, angle = 30, hjust = 1)) +
  nc_tag("b")

save_nc(p, file.path(OUT, "Fig5b_aminoacid_boxplots"), width_mm = NC_W2, height_mm = 56)

cat("\n-- panel statistics (two-tailed Welch t on mice; Storey q over the platform) --\n")
for (i in seq_len(nrow(PANELS)))
  cat(sprintf("  %-14s P = %-9s q = %-8s %s\n", PANELS$name[i],
              fmt_p(PANELS$P[i]), fmt_p(PANELS$q[i]),
              if (PANELS$P[i] < 0.05 && PANELS$q[i] < 0.10) "clears the Fig 1 bar" else "DOES NOT clear"))

src <- long
src$P <- PANELS$P[match(as.character(src$name), PANELS$name)]
src$q <- PANELS$q[match(as.character(src$name), PANELS$name)]
utils::write.csv(src, file.path(OUT, "SourceData_Fig5b.csv"), row.names = FALSE)
