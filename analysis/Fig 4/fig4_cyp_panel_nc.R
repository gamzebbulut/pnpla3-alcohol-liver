# =============================================================================
#  Figure 4 panel e -- the Cyp resolution
#
#  Panels b and c leave three pathways pointing the same way: arachidonic acid
#  metabolism and linoleic acid metabolism concordantly down in both arms, PPAR
#  signalling the most decreased pathway in the gene arm. They look like three
#  findings. They are largely ONE finding seen three times, because their KEGG
#  gene sets overlap heavily in the cytochrome P450s -- and it is the Cyps that
#  are down.
#
#  This panel is that claim, gene by gene: every Cyp belonging to any of the
#  three sets that was detected in the six hepatocyte populations, ordered by effect.
#
#  Layout is one ggplot, not a montage. The membership matrix is drawn in the
#  right-hand margin of the SAME x scale (x = 0.7/1.0/1.3, separated by a rule),
#  because patchwork/cowplot are not dependencies of this pipeline and a montage
#  assembled in Illustrator would drift out of row alignment on any re-render.
#
#  Colour is DIRECTION here, not category -- these are genes, not pathways, and
#  direction is the same red/blue that Figures 1-3 use. Saturation carries
#  significance: full tone for FDR < 0.05, tint for not. That avoids a separate
#  asterisk column and keeps the panel from screaming.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
CYP     <- Sys.getenv("FIG4_CYP",   file.path(FIGROOT, "Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv"))
OUT     <- Sys.getenv("FIG4_OUT",   file.path(FIGROOT, "Fig 4/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

FDR_CUT <- 0.05
## Which FDR shades the bars. "animal" is the replication unit the experiment
## randomised and must match panel e; "zone" is the stacked mouse x population model,
## kept only as a sensitivity view. See fig4_animal_level_limma.R.
FDR_SOURCE <- "animal"
OUTS <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
QUESTION <- "Q3"

cyp <- utils::read.csv(CYP, stringsAsFactors = FALSE)
cyp$AA <- as.logical(cyp$AA); cyp$LA <- as.logical(cyp$LA); cyp$PPAR <- as.logical(cyp$PPAR)

## The FDR column in shared_Cyps_AA_LA_PPAR_C1.csv came from the population-stacked
## model. Swap it for the animal-level one so panels d and e cannot disagree
## about which Cyps are significant.
if (FDR_SOURCE == "animal") {
  af <- file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION))
  if (!file.exists(af)) stop("Run fig4_animal_level_limma.R first, or set FDR_SOURCE <- \"zone\".")
  an <- utils::read.csv(af, stringsAsFactors = FALSE)
  cyp$FDR <- unname(stats::setNames(an$FDR_animal, an$SYMBOL)[cyp$symbol])
  cat(sprintf("  FDR source: animal level (n = mice), from %s\n", basename(af)))
} else {
  cat("  FDR source: population-stacked model -- SENSITIVITY VIEW, not the reported number\n")
}

## Undetected genes carry no effect to plot. They are counted in the caption
## line below rather than drawn as 27 empty rows.
n_all <- nrow(cyp)
d <- cyp[!is.na(cyp$logFC), ]
d <- d[order(d$logFC), ]
d$symbol <- factor(d$symbol, levels = d$symbol)
d$sig    <- !is.na(d$FDR) & d$FDR < FDR_CUT
d$key    <- ifelse(d$logFC < 0,
                   ifelse(d$sig, "down_sig", "down_ns"),
                   ifelse(d$sig, "up_sig",   "up_ns"))

BAR <- c(down_sig = NC_BLUE, down_ns = NC_BLUE_TINT,
         up_sig   = NC_RED,  up_ns   = NC_RED_TINT)

## ---- membership matrix, drawn into the right margin of the same x scale -----
MX  <- c(AA = 0.72, LA = 1.02, PPAR = 1.32)
RULE <- 0.50
mem <- do.call(rbind, lapply(names(MX), function(p) {
  data.frame(symbol = d$symbol, set = p, x = MX[[p]], present = d[[p]],
             stringsAsFactors = FALSE)
}))
mem <- mem[mem$present, ]

## The left half of the panel is empty by construction -- it is a ranked bar
## chart, so only the deepest bars reach it. That space carries the key, which
## the panel badly needed: light vs dark blue was undefined anywhere.
xlo <- min(d$logFC) * 1.10
xhi <- 1.52
KX  <- xlo * 0.97; KW <- abs(xlo) * 0.09; KY <- nrow(d) - 0.2   # key anchor
brk <- pretty(c(xlo, 0), n = 4); brk <- brk[brk >= xlo & brk <= 0.4]

p <- ggplot(d, aes(logFC, symbol)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = NC_INK) +
  geom_col(aes(fill = key), width = 0.68, colour = NC_INK, linewidth = 0.1) +
  ## separator between the effect axis and the membership block
  annotate("segment", x = RULE, xend = RULE, y = 0.4, yend = nrow(d) + 1.2,
           linewidth = 0.25, colour = "#CCCCCC") +
  ## empty slots, so an absent membership reads as an empty cell not a gap
  geom_point(data = expand.grid(symbol = d$symbol, x = unname(MX)),
             aes(x = x, y = symbol), inherit.aes = FALSE,
             shape = 21, size = 1.5, fill = NA, colour = "#E0E0E0", stroke = 0.25) +
  geom_point(data = mem, aes(x = x, y = symbol), inherit.aes = FALSE,
             shape = 21, size = 1.5, fill = NC_INK2, colour = NC_INK, stroke = 0.2) +
  annotate("text", x = unname(MX), y = nrow(d) + 1.6, label = names(MX),
           size = NC_MIN / .pt, colour = NC_INK2, vjust = 0) +
  ## ---- key, drawn into the empty left half in data coordinates -----------
  annotate("rect", xmin = KX, xmax = KX + KW, ymin = KY - 0.30, ymax = KY + 0.30,
           fill = NC_BLUE, colour = NC_INK, linewidth = 0.1) +
  annotate("rect", xmin = KX, xmax = KX + KW, ymin = KY - 1.60, ymax = KY - 1.00,
           fill = NC_BLUE_TINT, colour = NC_INK, linewidth = 0.1) +
  annotate("rect", xmin = KX, xmax = KX + KW, ymin = KY - 2.90, ymax = KY - 2.30,
           fill = NC_RED_TINT, colour = NC_INK, linewidth = 0.1) +
  annotate("text", x = KX + KW + 0.10, y = c(KY, KY - 1.30, KY - 2.60), hjust = 0,
           size = NC_MIN / .pt, colour = NC_INK,
           label = c(sprintf("decreased, FDR < %.2f", FDR_CUT),
                     "decreased, not significant",
                     "increased, not significant")) +
  annotate("text", x = KX, y = KY - 4.1, hjust = 0, vjust = 1,
           size = NC_MIN / .pt, colour = NC_INK2, lineheight = 0.9,
           label = sprintf("%d of %d Cyps reach FDR < %.2f;\nall of them decreased.\nNo Cyp increased significantly.",
                           sum(d$sig), nrow(d), FDR_CUT)) +
  annotate("text", x = KX, y = KY - 8.0, hjust = 0, vjust = 1,
           size = NC_MIN / .pt, colour = NC_INK2, lineheight = 0.9, fontface = "italic",
           label = sprintf("%d further Cyps in these sets\nwere not detected.", n_all - nrow(d))) +
  scale_fill_manual(values = BAR, guide = "none") +
  scale_x_continuous(breaks = brk, limits = c(xlo, xhi), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression("Log"[2]~"fold-change, alcohol on variant"), y = NULL) +
  theme_nc(grid_y = FALSE) +
  theme(axis.text.y = element_text(size = NC_MIN, face = "italic"),
        plot.margin = margin(8, 2, 2, 2)) +
  nc_tag("e")

save_nc(p, file.path(OUT, "Fig4e_cyp_panel"), width_mm = 72, height_mm = 96)

utils::write.csv(d[, c("symbol","AA","LA","PPAR","sets","n_pathway","logFC","FDR","direction")],
                 file.path(OUT, "SourceData_Fig4e_cyps.csv"), row.names = FALSE)

## ---- the numbers the legend and Results need -------------------------------
cat(sprintf("\n  Cyps in AA / LA / PPAR gene sets : %d (union)\n", n_all))
cat(sprintf("  detected in six hepatocyte populations    : %d  (%d undetected, not drawn)\n",
            nrow(d), n_all - nrow(d)))
cat(sprintf("  decreased                        : %d of %d\n", sum(d$logFC < 0), nrow(d)))
cat(sprintf("  FDR < %.2f (%s)              : %d, all %s\n", FDR_CUT, FDR_SOURCE, sum(d$sig),
            if (all(d$logFC[d$sig] < 0)) "decreased" else "mixed"))
cat(sprintf("  in all three sets                : %d\n", sum(d$n_pathway == 3)))
cat("  strongest: ")
cat(paste(sprintf("%s (%.2f)", d$symbol[order(d$logFC)][1:3], sort(d$logFC)[1:3]), collapse = ", "), "\n")
for (s in unique(cyp$sets)) cat(sprintf("    %-8s n = %2d (detected %2d)\n", s,
    sum(cyp$sets == s), sum(cyp$sets == s & !is.na(cyp$logFC))))
