## ════════════════════════════════════════════════════════════════════
##  LEGENDS — two standalone pieces for the volcano panels
##    1) class-color legend (colorful names) + a grey "set aside" entry
##    2) -log10(p) size legend  (dots match the plot because the size
##       scale below is IDENTICAL to the volcano script)
##  No data needed — built straight from the palette + size scale.
##  Hit Source; crop the white space in Illustrator.
##
##  IMPORTANT: keep PT_RANGE / SIZE_LIM / SIZE_BREAKS / SIZE_LABELS the
##  same here and in the volcano script, or the legend dots won't match.
## ════════════════════════════════════════════════════════════════════

OUT_DIR <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/F1"

## ── Size scale (MUST match the volcano script) ──────────────────────
PT_RANGE    <- c(2, 7)
SIZE_LIM    <- c(1.3, 4.6)
SIZE_BREAKS <- c(1.3, 2, 3, 4)
SIZE_LABELS <- c("0.05","0.01","0.001","0.0001")

## ── Class-color legend ──────────────────────────────────────────────
## List the colored (non-set-aside) classes you want shown, top order.
CLASS_LIST <- c(
  "Amino acids & derivatives", "Bile acids", "Disaccharides", "Fatty acids",
  "Fatty esters", "Flavin Nucleotides", "Glycerolipids", "Purine Nucleotides",
  "Sphingolipids", "Sterols", "Sulfinic Acids and Derivatives")
SHOW_GRAY_ENTRY <- TRUE                              # explain Option B's grey dots
GRAY_LABEL      <- "Glycerophospholipids (set aside)"

CLASS_NROW <- 3        # rows in the (horizontal) class legend
FONT_LEG   <- 13       # legend text size
CLASS_W    <- 9.0; CLASS_H <- 2.2     # class legend canvas (inches)
SIZE_W     <- 1.9; SIZE_H  <- 3.0     # size  legend canvas (inches)
## ════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({ library(ggplot2) })
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_NS <- "#BBBBBB"; MUTE_COL <- "#9E9E9E"
class_pal <- c(
  "Amino acids & derivatives"="#E53935","Bile acids"="#F57F17","Disaccharides"="#FF8F00",
  "Fatty acids"="#558B2F","Fatty esters"="#33691E","Flavin Nucleotides"="#00695C",
  "Glycerolipids"="#42A5F5","Glycerophospholipids"="#1565C0","Purine Nucleotides"="#6A1B9A",
  "Sphingolipids"="#9C27B0","Sterols"="#880E4F","Sulfinic Acids and Derivatives"="#BF360C")

## ── 1) Class-color legend ───────────────────────────────────────────
items <- CLASS_LIST
cols  <- class_pal[CLASS_LIST]
if (SHOW_GRAY_ENTRY) { items <- c(items, GRAY_LABEL); cols <- c(cols, setNames(MUTE_COL, GRAY_LABEL)) }
class_legend <- ggplot(data.frame(cls = factor(items, levels = items), x = 1),
                       aes(x, cls, colour = cls)) +
  geom_point(alpha = 0) +
  scale_colour_manual(values = cols, limits = items, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1),
                               nrow = CLASS_NROW, byrow = TRUE)) +
  theme_void(base_size = 14) +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.text = element_text(size = FONT_LEG))
ggsave(file.path(OUT_DIR, "legend_classes.png"), class_legend, width = CLASS_W, height = CLASS_H, dpi = 300)
ggsave(file.path(OUT_DIR, "legend_classes.svg"), class_legend, width = CLASS_W, height = CLASS_H, dpi = 300)

## ── 2) -log10(p) size legend (dots match the plot) ──────────────────
size_legend <- ggplot(data.frame(x = 1, y = SIZE_BREAKS, s = SIZE_BREAKS),
                      aes(x, y, size = s)) +
  geom_point(alpha = 0) +
  scale_size_continuous(range = PT_RANGE, limits = SIZE_LIM,
                        breaks = SIZE_BREAKS, labels = SIZE_LABELS, name = "-log10(p)") +
  guides(size = guide_legend(override.aes = list(alpha = 1, colour = "#555555"), ncol = 1)) +
  theme_void(base_size = 14) +
  theme(legend.position = "right",
        legend.text = element_text(size = FONT_LEG),
        legend.title = element_text(size = FONT_LEG + 1, face = "bold"))
ggsave(file.path(OUT_DIR, "legend_size_logp.png"), size_legend, width = SIZE_W, height = SIZE_H, dpi = 300)
ggsave(file.path(OUT_DIR, "legend_size_logp.svg"), size_legend, width = SIZE_W, height = SIZE_H, dpi = 300)

cat("Saved legend_classes.png/.svg and legend_size_logp.png/.svg to", OUT_DIR, "\n")
