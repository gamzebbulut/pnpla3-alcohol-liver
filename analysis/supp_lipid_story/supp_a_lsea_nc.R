# =============================================================================
#  Combined lipid supplement, panel a — LSEA by chain length and unsaturation
#
#  WHY NOT THE CLASS SETS. Figure 1d already carries the class sets for the two
#  alcohol contrasts. Repeating them here would put the same result in two
#  figures; what Figure 1d deliberately left alone is the chain-length and
#  unsaturation half, and that is this panel. CLASS_TOO <- TRUE puts the class
#  block back if the supplement is ever wanted as a complete superset -- in
#  which case say so in the legend.
#
#  ORDER IS NUMERIC, NOT BY EFFECT. Chain length and double-bond count are
#  ordinal, so sorting them by enrichment score would destroy the only thing the
#  panel is for: reading the response as a profile along the axis. Ordered
#  numerically, the alcohol rows show short chains up and long chains down, and
#  the variant rows show the mirror -- which is the claim panels b and c then
#  quantify.
#
#  NO C1-C4 (revision task A2). Rows are named in plain language, in
#  nc_style.R's CONTRAST_ORDER.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
SIG        <- 0.10       # LSEA padj; the alcohol-on-variant arm clears nothing at 0.05
CLASS_TOO  <- FALSE
ARMS       <- c(Q4 = "Alcohol on wild type", Q3 = "Alcohol on variant",
                Q2 = "Variant alone",       Q5 = "Variant on alcohol")
W_MM       <- 180
H_MM       <- 46
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("SUPP_LSEA", file.path(FIGROOT, "Fig 3/Lipid_LSEA_Q3_vs_Q4",
                                             "lsea_ALLQ_comparison.csv"))
OUT     <- Sys.getenv("SUPP_OUT",  file.path(FIGROOT, "supp_lipid_story/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

if (!file.exists(IN)) stop("LSEA table not found: ", IN)
z <- utils::read.csv(IN, stringsAsFactors = FALSE)

## Set names arrive as Class_TG / total_cl_56 / total_cs_8. Split them into a
## block and a label; the numeric part is kept as a number so the axis sorts.
block <- ifelse(grepl("^total_cl_", z$set), "Acyl chain length (carbons)",
         ifelse(grepl("^total_cs_", z$set), "Double bonds",
                "Lipid class"))
value <- suppressWarnings(as.numeric(sub("^total_c[ls]_", "", z$set)))
label <- ifelse(block == "Lipid class", sub("^Class_", "", z$set), as.character(value))

D <- do.call(rbind, lapply(names(ARMS), function(q) data.frame(
  block = block, label = label, value = value, arm = unname(ARMS[q]),
  NES = z[[paste0("NES_", q)]], padj = z[[paste0("padj_", q)]],
  stringsAsFactors = FALSE)))
if (!CLASS_TOO) D <- D[D$block != "Lipid class", ]
D <- D[!is.na(D$NES), ]
stopifnot("No sets left" = nrow(D) > 0)
D$sig <- !is.na(D$padj) & D$padj < SIG

BLOCKS <- c("Acyl chain length (carbons)", "Double bonds")
if (CLASS_TOO) BLOCKS <- c("Lipid class", BLOCKS)
D$block <- factor(D$block, levels = BLOCKS)
## numeric order within each block; classes (if drawn) keep alphabetical
lev <- unique(D[order(D$block, D$value, D$label), c("block", "label")])
D$label <- factor(D$label, levels = unique(lev$label))
D$arm   <- factor(D$arm, levels = rev(unname(ARMS)))

LIM <- max(abs(D$NES), na.rm = TRUE)
p <- ggplot(D, aes(label, arm)) +
  geom_tile(aes(fill = NES), colour = "white", linewidth = 0.35) +
  geom_tile(data = D[D$sig, ], fill = NA, colour = NC_INK, linewidth = 0.5) +
  facet_grid(~ block, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                       limits = c(-LIM, LIM), na.value = "#F2F2F2",
                       name = "Normalized enrichment score",
                       breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(expand = expansion(add = 0.04)) +
  ## Reader-facing contrast names from nc_style.R's display layer, parsed for
  ## the I148M superscript. The data keeps the internal keys.
  scale_y_discrete(expand = expansion(add = 0.04), labels = nc_label_expr) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(34, "pt"),
                                barheight = grid::unit(5.5, "pt"),
                                title.position = "left", ticks.colour = NC_INK)) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = NC_BIG),
        axis.text.y = element_text(size = NC_BIG),
        strip.text = element_text(size = NC_BIG),
        panel.spacing = grid::unit(6, "pt"),
        legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = NC_BIG),
        legend.text = element_text(size = NC_BIG),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "pt")) +
  nc_tag("a")

save_nc(p, file.path(OUT, "Supp_lipid_a_lsea"), W_MM, H_MM)
utils::write.csv(D, file.path(OUT, "SourceData_Supp_lipid_a_lsea.csv"), row.names = FALSE)

cat(sprintf("\n  %d sets x %d contrasts drawn (%s)\n", nlevels(D$label), nlevels(D$arm),
            paste(BLOCKS, collapse = " + ")))
for (a in rev(levels(D$arm)))
  cat(sprintf("  %-22s %2d of %2d sets at padj < %.2f\n", a,
              sum(D$sig[D$arm == a]), sum(D$arm == a), SIG))
cat("\n  LEGEND MATERIAL: nothing here clears padj < 0.05 in the variant arm, so\n")
cat(sprintf("  the threshold drawn is %.2f and the legend must say so.\n", SIG))
