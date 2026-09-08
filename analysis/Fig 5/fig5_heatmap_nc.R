# =============================================================================
#  Figure 5 — integrated pathways across the four comparisons, NC style
#
#  Replaces fig5_concordance_heatmap.R (heatmap section only; the graphical-
#  abstract section of that file is untouched and still lives there).
#
#  RESTYLE ONLY. The significance rule, the input run, the row set and the row
#  order are all unchanged, so the cells that were filled before are filled now.
#  What changed is how it is drawn:
#
#    - COLOUR NOW MEANS DIRECTION, as in Figures 1-4. Previously blue meant
#      "gene arm" and red meant "metabolite arm", so a blue triangle pointing up
#      was a gene-arm increase -- and a reader arriving from four figures where
#      blue has meant DECREASED reads it backwards. Arm is now carried by
#      position (left glyph = gene, right glyph = metabolite), which was already
#      true, and the freed colour channel carries direction: red up, blue down.
#      Discordance becomes a colour mismatch inside one cell, which is far
#      easier to see than comparing two triangle orientations.
#    - Cell shading moved off the red ramp onto a neutral one. Red now means
#      increased; a red cell behind red and blue glyphs was a third meaning for
#      the same hue. The shading still encodes -log10(P) from the metabolite arm,
#      the limiting arm, so nothing is lost.
#    - 180 mm, 5-7 pt type, vector SVG through nc_style.R.
#    - Contrast labels in plain language (task A2), from CONTRAST_LABEL, with a
#      "gene / metab." sub-header so the left-right convention is readable at the
#      cell rather than only in the caption.
#    - ggtext dropped, and with it the install.packages() call that ran on load.
#      A script that installs packages as a side effect of being sourced is not
#      reproducible.
#    - The triangles are ggplot POINT SHAPES (24 up, 25 down), not the text
#      characters U+25B2/U+25BC. Text glyphs depend on the font having them: in
#      a font that does not, they silently fall back to dots, and the same class
#      of failure has already bitten this pipeline twice (the missing-glyph
#      boxes in Figure 2). Point shapes are drawn as vector paths and cannot go
#      missing, in any font or in Illustrator.
#
#  KNOWN LIMITATION, deliberately left as it is pending a decision:
#  a white cell means "did not pass the joint threshold", NOT "no effect". Some
#  rows carry the same direction in a comparison where they miss significance --
#  arachidonic acid metabolism is down in both arms in C1, C2 AND C3 but only
#  reaches the threshold in C1. Reading the blanks as background specificity is
#  the difference-in-significance fallacy. Shading every cell by Z and marking
#  significance with a border instead would remove the trap; that is a content
#  change, not a restyle, so it is not done here. Say it in the legend meanwhile.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN_BASE <- Sys.getenv("FIG5_IN",    file.path(ROOT, "GB_metabolomics_20260628"))
OUT     <- Sys.getenv("FIG5_OUT",   file.path(FIGROOT, "Fig 5/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

GENE_FDR <- 0.10
MET_P    <- 0.05
qmap  <- c(C1 = "Q3", C2 = "Q4", C3 = "Q2", C4 = "Q5")
## The columns are keyed C1..C4 internally. ckey turns each into the internal
## contrast key, and nc_label_expr turns THAT into what the reader sees --
## nc_style.R's display layer, two-line form, with the I148M superscript.
## Nothing here is ever drawn as C1..C4 (revision task A2).
ckey  <- c(C1 = "Alcohol on variant", C2 = "Alcohol on wild type",
           C3 = "Variant alone",      C4 = "Variant on alcohol")
X_LAB <- function(x) nc_label_expr(ckey[as.character(x)], two_line = TRUE)

long <- do.call(rbind, lapply(names(qmap), function(cc) {
  f <- file.path(IN_BASE, paste0(qmap[[cc]], "_StoufferPMAX_pooled"),
                 "modalities_pmax_gene_plus_metDec12.csv")
  if (!file.exists(f)) stop("Missing ", f)
  z <- utils::read.csv(f, stringsAsFactors = FALSE)
  data.frame(KEGG_id = z$KEGG_id,
             pw_label = ifelse(is.na(z$KEGG_title) | z$KEGG_title == "", z$KEGG_id, z$KEGG_title),
             contrast = cc, Z_gene = z$Z_gene, FDR_gene = z$FDR_gene,
             Z_met = z$Z_met, P_met = z$P_met, stringsAsFactors = FALSE)
}))

long$pw_label <- gsub(" metabolism",  " metab.",    long$pw_label, fixed = TRUE)
long$pw_label <- gsub(" biosynthesis", " biosynth.", long$pw_label, fixed = TRUE)
long$sig <- !is.na(long$FDR_gene) & !is.na(long$P_met) &
            long$FDR_gene < GENE_FDR & long$P_met < MET_P
long$contrast <- factor(long$contrast, levels = names(qmap))

ids <- unique(long$KEGG_id[long$sig])
hd  <- long[long$KEGG_id %in% ids, ]

## Shading: strength of the metabolite-arm P, the limiting arm. The original
## mapped -log10(P) onto a continuous ramp whose legend read "0 / 2 / 4", which
## tells a reader nothing -- nobody converts -log10 in their head. Same
## information, now in three DISCRETE bands labelled with the actual P values.
## Non-significant cells stay white, unchanged.
## "not significant" is a REAL level, not NA. Coding the white cells as NA left
## the reader's key silent about the most common cell in the figure, and on some
## ggplot2 versions it also puts an "NA" entry in the legend -- which would imply
## missing measurements. There are none: all cells here are measured in both
## arms, and the stopifnot below fails loudly if that ever stops being true
## rather than drawing an unmeasured pathway as if it were merely non-significant.
BANDS <- c("not significant", "P \u2264 0.05", "P \u2264 0.01", "P \u2264 0.001")
stopifnot("Some cells are not measured in both arms -- white would be ambiguous" =
          all(!is.na(hd$FDR_gene) & !is.na(hd$P_met)))
hd$band <- ifelse(!hd$sig, BANDS[1],
           ifelse(hd$P_met <= 0.001, BANDS[4],
           ifelse(hd$P_met <= 0.01,  BANDS[3], BANDS[2])))
hd$band <- factor(hd$band, levels = BANDS)
cat(sprintf("  all %d cells measured in both arms; no missing data\n", nrow(hd)))

## Glyphs. Position is the arm, colour is the direction.
gl <- function(z, dx, arm) {
  s <- hd[hd$sig & !is.na(z), ]
  zz <- z[hd$sig & !is.na(z)]
  data.frame(contrast = s$contrast, pw_label = s$pw_label, arm = arm,
             x = as.numeric(s$contrast) + dx,
             shp  = ifelse(zz > 0, 24L, 25L),          # 24 = up, 25 = down
             dir  = ifelse(zz > 0, "increased", "decreased"), stringsAsFactors = FALSE)
}
G <- rbind(gl(hd$Z_gene, -0.19, "gene"), gl(hd$Z_met, +0.19, "metabolite"))

## Row order: lowest gene-arm FDR among that row's significant cells at the top.
## Unchanged from the original script.
o   <- hd[hd$sig, ]
ord <- tapply(o$FDR_gene, o$pw_label, min)
lev <- names(sort(ord, decreasing = TRUE))          # rev() of ascending -> top = smallest
hd$pw_label <- factor(hd$pw_label, levels = lev)
G$pw_label  <- factor(G$pw_label,  levels = lev)

n <- length(lev)
cat(sprintf("  %d pathways shown; %d of %d cells significant (%.0f%%)\n",
            n, sum(hd$sig), nrow(hd), 100 * sum(hd$sig) / nrow(hd)))
print(table(contrast = hd$contrast[hd$sig]))

## Arm sub-header. The left/right convention was only in the caption, which is
## the wrong place for something a reader needs while looking at a cell. These
## sit at the exact x offsets the glyphs use, directly above the matrix.
armhdr <- data.frame(
  x = rep(seq_along(levels(hd$contrast)), each = 2) + c(-0.19, 0.19),
  lab = rep(c("gene", "metab."), times = nlevels(hd$contrast)),
  stringsAsFactors = FALSE)

p <- ggplot(hd, aes(contrast, pw_label)) +
  geom_tile(aes(alpha = band), fill = "#757575", colour = "#BDBDBD", linewidth = 0.15) +
  geom_point(data = G, aes(x = x, y = pw_label, shape = shp, fill = dir),
             colour = NC_INK, stroke = 0.2, size = 2.4, inherit.aes = FALSE) +
  scale_shape_identity() +
  ## Two fill scales are not possible in one ggplot, so the tile bands are drawn
  ## through `alpha` on a fixed grey and `fill` is reserved for the glyphs.
  scale_alpha_manual(values = stats::setNames(c(0, 0.16, 0.34, 0.55), BANDS),
                     drop = FALSE, na.translate = FALSE, name = "Metabolite arm") +
  scale_fill_manual(values = c(increased = NC_RED, decreased = NC_BLUE), name = NULL,
                    labels = c(increased = "Increased", decreased = "Decreased")) +
  geom_text(data = armhdr, aes(x = x, y = n + 0.72, label = lab), inherit.aes = FALSE,
            size = NC_MIN / .pt, colour = NC_INK2, vjust = 0) +
  scale_x_discrete(labels = X_LAB, position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL,
       caption = paste0("Not significant means the cell did not pass BOTH arms\n",
                        "(gene FDR < 0.10 and metabolite P < 0.05). It is not the\n",
                        "same as no effect: every cell here was measured.")) +
  guides(fill  = guide_legend(order = 1, override.aes = list(shape = 22, size = 2.4)),
         alpha = guide_legend(order = 2, override.aes = list(fill = "#757575", colour = "#BDBDBD"),
                              keywidth = grid::unit(8, "pt"), keyheight = grid::unit(8, "pt"))) +
  theme_nc() +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = NC_BASE, lineheight = 0.9,
                                       margin = margin(b = 2)),
        axis.text.y = element_text(size = NC_BASE),
        legend.position = "right", legend.key.size = grid::unit(7, "pt"),
        plot.caption = element_text(size = NC_MIN, colour = NC_INK2, hjust = 0,
                                    margin = margin(t = 3))) +
  nc_tag("a")

## Height is set so that a + b fit the 225 mm page together: this panel at
## 135 mm plus the amino-acid panel at 72 mm leaves ~18 mm for the gutter
## and the panel letters.
H <- min(NC_HMAX, 26 + 5.2 * n)
save_nc(p, file.path(OUT, "Fig5_integrated_pathways"), width_mm = NC_W2, height_mm = H)

src <- hd[, c("KEGG_id", "pw_label", "contrast", "Z_gene", "FDR_gene", "Z_met", "P_met", "sig")]
src <- src[order(src$pw_label, src$contrast), ]
utils::write.csv(src, file.path(OUT, "SourceData_Fig5.csv"), row.names = FALSE)
