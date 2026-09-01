# =============================================================================
#  Figure 4 panel b -- gene arm against metabolite arm: the concordance scatter
#
#  Every pathway measured in BOTH arms is a point: gene-arm Z on x, metabolite-
#  arm Z on y. This is the standard idiom for two-modality concordance, and it
#  carries in one panel what the bars and the heatmap needed three to say:
#
#    - agreement is the diagonal. Lower-left and upper-right quadrants are
#      concordant; upper-left and lower-right are discordant, by construction.
#    - colour is pathway category, so the amino-acid cluster reads as a cluster
#      rather than as adjacent rows in a ranked list.
#    - SIZE is joint magnitude, sqrt(Z_gene^2 + Z_met^2) -- how far the pathway
#      sits from the origin, made legible at a glance. Position already carries
#      this, so size is deliberately redundant reinforcement, not new
#      information; it exists so the strong hits pop out of the 90 grey points
#      in a 72 mm-wide panel where distance is hard to judge by eye.
#
#  Quadrants are labelled rather than shaded -- shading would add two more
#  colours to a figure that already has enough.
#
#  NOTE no y = x line. The two Z scales are not comparable: the gene arm is a
#  CT-weighted Stouffer combination over cell types (spans +/-6.7), the
#  metabolite arm a single limma moderated t (spans +/-2.8). The diagonal would
#  invite a reader to compare magnitudes that are not on one scale. Sign
#  agreement, i.e. the quadrant, is the claim -- not slope.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING -- label legibility (revised 31 Aug)
## ---------------------------------------------------------------------------
## The 14 labels used to be set at NC_DENSE (5 pt) because they would not fit
## at the 6 pt floor. They fit now, and the fix was not type size or panel
## width -- it was giving ggrepel somewhere to put them:
##
##   HEADROOM_Y  the seven amino-acid pathways all land in one corner at
##               y ~ 2.2-2.5. At the old 1.26 multiplier the top of the panel
##               ended just above them and repel had no room, so it stacked
##               labels on top of each other. 1.55 opens a band the labels can
##               occupy without covering points.
##   WRAP_W      13 characters, not 17. Narrow tall blocks pack into a crowded
##               horizontal band; wide short ones do not.
##
## Checked by rendering, not by arithmetic: at 74 x 86 with these values all 14
## labels place at 6 pt with no overlaps. 86 mm is cleaner still -- fewer labels
## sit over background circles -- if panel a's diagram can give up the width.
LAB_PT     <- NC_MIN   # 6 pt, the manuscript floor. Was NC_DENSE (5 pt).
HEADROOM_Y <- 1.55
WRAP_W     <- 13
PANEL_W    <- 74       # mm; 86 reduces label-over-point overprinting
## The shared category key, emitted as its own strip. Sized to sit under panels
## b and c together: 74 + 4 mm gutter + 54 = 132 mm.
KEY_W      <- 132
KEY_H      <- 9    # two rows of 8 pt; anything taller is empty space
KEY_NROW   <- 2

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("FIG4_IN",    file.path(ROOT, "GB_metabolomics_20260628/Q3_StoufferPMAX_pooled"))
JOINT   <- Sys.getenv("FIG4_JOINT", file.path(FIGROOT, "Fig 4/C1_joint_hits.csv"))
CAT     <- Sys.getenv("FIG4_CATS",  file.path(FIGROOT, "Fig 4/pathway_categories.csv"))
OUT     <- Sys.getenv("FIG4_OUT",   file.path(FIGROOT, "Fig 4/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

combo <- utils::read.csv(file.path(IN, "modalities_pmax_gene_plus_metDec12.csv"), stringsAsFactors = FALSE)
cats  <- utils::read.csv(CAT, stringsAsFactors = FALSE)
catof <- stats::setNames(cats$category, cats$KEGG_id)
J     <- utils::read.csv(JOINT, stringsAsFactors = FALSE)

d <- combo[!is.na(combo$Z_gene) & !is.na(combo$Z_met), ]
d$title <- ifelse(is.na(d$KEGG_title) | d$KEGG_title == "", d$KEGG_id, d$KEGG_title)
d$joint <- d$KEGG_id %in% J$KEGG_id
d$cat   <- unname(catof[d$KEGG_id]); d$cat[is.na(d$cat)] <- "Other"
d$cat   <- factor(d$cat, levels = names(NC_CAT))
d$mag   <- sqrt(d$Z_gene^2 + d$Z_met^2)          # joint magnitude -> point size

hit <- d[d$joint, ]
## Generous corner whitespace: the quadrant labels live out there, and
## ggrepel is fenced out of it (xlim/ylim below) so labels never collide
## with them.
lim_x <- max(abs(d$Z_gene)) * 1.10
lim_y <- max(abs(d$Z_met))  * HEADROOM_Y

## Labels are long. Wrap them so ggrepel has compact blocks to place rather
## than 50-character ribbons that cannot fit anywhere in a 72 mm panel.
wrap_lab <- function(x, width = WRAP_W) vapply(x, function(s)
  paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
hit$lab <- wrap_lab(hit$title)

cat(sprintf("  %d pathways measured in both arms; %d are joint hits\n", nrow(d), nrow(hit)))
q <- c(sum(hit$Z_gene > 0 & hit$Z_met > 0), sum(hit$Z_gene < 0 & hit$Z_met < 0),
       sum(hit$Z_gene < 0 & hit$Z_met > 0), sum(hit$Z_gene > 0 & hit$Z_met < 0))
cat(sprintf("  concordant up %d, concordant down %d, discordant %d\n", q[1], q[2], q[3] + q[4]))
cat(sprintf("  joint magnitude range %.2f - %.2f\n", min(hit$mag), max(hit$mag)))

p <- ggplot(d, aes(Z_gene, Z_met)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#BBBBBB") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BBBBBB") +
  ## Background pathways are drawn HOLLOW, not grey-filled. "Other" is itself a
  ## grey category (NC_CAT), so a grey fill had to mean one thing only; fill vs
  ## no fill separates "not a joint hit" from "joint hit, uncategorised"
  ## without spending another hue.
  geom_point(data = d[!d$joint, ], aes(size = mag), shape = 21, fill = NA,
             colour = "#D0D0D0", stroke = 0.3) +
  geom_point(data = hit, aes(fill = cat, size = mag), shape = 21,
             colour = NC_INK, stroke = 0.2) +
  ## LAB_PT, HEADROOM_Y and WRAP_W are set together at the top of this file --
  ## changing one without the others is what put labels on top of each other.
  ggrepel::geom_text_repel(data = hit, aes(label = lab, colour = cat),
                           size = LAB_PT / .pt, seed = 1, max.overlaps = Inf,
                           lineheight = 0.85,
                           min.segment.length = 0, segment.size = 0.15,
                           segment.colour = "#AAAAAA", box.padding = 0.28,
                           point.padding = 0.25, force = 6, force_pull = 0.4,
                           max.iter = 30000, show.legend = FALSE,
                           xlim = c(-lim_x * 0.98, lim_x * 0.98),
                           ylim = c(-lim_y * 0.95, lim_y * 0.95)) +
  annotate("text", x = lim_x, y = lim_y, hjust = 1, vjust = 1,
           label = "concordant\nup", size = NC_MIN / .pt, colour = NC_INK2, lineheight = 0.9) +
  annotate("text", x = -lim_x, y = -lim_y, hjust = 0, vjust = 0,
           label = "concordant\ndown", size = NC_MIN / .pt, colour = NC_INK2, lineheight = 0.9) +
  annotate("text", x = lim_x, y = -lim_y, hjust = 1, vjust = 0,
           label = "discordant", size = NC_MIN / .pt, colour = NC_INK2) +
  annotate("text", x = -lim_x, y = lim_y, hjust = 0, vjust = 1,
           label = "discordant", size = NC_MIN / .pt, colour = NC_INK2) +
  scale_fill_manual(values = NC_CAT, name = NULL, drop = FALSE) +
  scale_colour_manual(values = NC_CAT, guide = "none") +
  scale_size_area(name = NULL, max_size = 2.6, breaks = c(2, 4, 6),
                  labels = c("2", "4", "6")) +
  coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
  labs(x = "Gene arm Z   (Stouffer, snRNA-seq)",
       y = "Metabolite arm Z   (ssGSEA)") +
  ## Category colours are shared with panel c, so the category key is emitted
  ## ONCE as its own strip (below) and placed under the row in Illustrator.
  ## Only the size key, which is specific to this panel, rides along here.
  guides(fill = "none",
         size = guide_legend(nrow = 1, title = "Joint magnitude",
                             override.aes = list(fill = "#BDBDBD", shape = 21,
                                                 colour = NC_INK, stroke = 0.2))) +
  theme_nc(grid_y = TRUE) +
  theme(legend.position = "top",
        legend.title = element_text(size = NC_MIN, colour = NC_INK2),
        legend.box.spacing = grid::unit(2, "pt"),
        legend.key.size = grid::unit(6, "pt")) +
  nc_tag("b")

save_nc(p, file.path(OUT, "Fig4b_concordance_scatter"), width_mm = PANEL_W, height_mm = 86)

## ---------------------------------------------------------------------------
## Shared category key for the whole figure row. One strip, placed once.
## ---------------------------------------------------------------------------
## "not a joint hit" is carried as a sixth swatch with a transparent fill and a
## pale outline, so the open circles in panel b are defined in the same strip.
KEY_V <- c(NC_CAT, "Not a joint hit" = NA)
key <- data.frame(cat = factor(names(KEY_V), levels = names(KEY_V)), x = 1, y = 1)
## alpha = 0 on the layer: the six key points are all plotted at (1, 1) and
## overprint into a single grey dot floating under the strip, which is what the
## placed key showed on 30 Aug. The legend needs the layer to exist, not to be
## visible, so it is drawn transparent and override.aes puts the swatches back.
pk <- ggplot(key, aes(x, y, fill = cat, colour = cat)) +
  geom_point(shape = 21, size = 2, stroke = 0.25, alpha = 0) +
  scale_fill_manual(values = KEY_V, name = NULL, drop = FALSE, na.value = NA) +
  scale_colour_manual(values = c(stats::setNames(rep(NC_INK, length(NC_CAT)), names(NC_CAT)),
                                 "Not a joint hit" = "#D0D0D0"), guide = "none") +
  ## Two rows: six entries totalling ~95 characters do not fit on one row at
  ## 132 mm and 8 pt, and shortening the category names to make them fit would
  ## make the key disagree with the pathway_categories.csv it describes.
  guides(fill = guide_legend(nrow = KEY_NROW, byrow = TRUE, override.aes = list(
           alpha = 1,
           colour = c(rep(NC_INK, length(NC_CAT)), "#D0D0D0")))) +
  theme_void(base_size = NC_BASE) +
  theme(text = element_text(family = FONT, size = NC_BASE, colour = NC_INK),
        ## Left-justified with the margins zeroed so the strip's first swatch
        ## aligns with panel b's left edge when it is placed under the row.
        ## Default justification centres the key and it floats off-axis.
        legend.position = "top", legend.justification = "left",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(0, "pt"),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "pt"),
        legend.key.size = grid::unit(6, "pt"),
        legend.text = element_text(size = NC_BASE, colour = NC_INK))
save_nc(pk, file.path(OUT, "Fig4_category_key"), width_mm = KEY_W, height_mm = KEY_H)

utils::write.csv(d[, c("KEGG_id","title","cat","Z_gene","FDR_gene","Z_met","P_met","mag","joint")],
                 file.path(OUT, "SourceData_Fig4b_scatter.csv"), row.names = FALSE)
