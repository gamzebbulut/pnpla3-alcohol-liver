# =============================================================================
#  EXPLORATORY — the two alcohol arms side by side
#
#  Not a manuscript figure. This exists to look at something the 30 Aug run
#  turned up that the manuscript does not currently say: in the transcriptome,
#  alcohol on the WILD-TYPE background moves far more genes than alcohol on the
#  variant background, in every cell type tested. The metabolome says the
#  opposite (107 significant on the variant against 32 on the wild type). Before
#  any of that reaches the text it needs looking at, and this is the looking.
#
#  Panel a  every cell type, the three contrasts, counted at FDR < 0.05.
#  Panel b  the Cyps across the six hepatocyte populations, one block per arm,
#           on ONE colour scale so the two blocks are directly comparable.
#
#  INPUT: fig4_outs_ALLCELLS/cyp_specificity_qlf_long.csv and
#         cyp_specificity_qlf_celltype_summary.csv, from fig4_cyp_specificity_qlf.R
# =============================================================================

FDR_CUT <- 0.05
MIN_POP <- 4     # a Cyp is drawn if it was detected in this many populations

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs_ALLCELLS"))
OUT     <- Sys.getenv("FIG4_OUT",  file.path(FIGROOT, "Fig 4/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S <- utils::read.csv(file.path(OUTS, "cyp_specificity_qlf_celltype_summary.csv"),
                     stringsAsFactors = FALSE)
L <- utils::read.csv(file.path(OUTS, "cyp_specificity_qlf_long.csv"),
                     stringsAsFactors = FALSE)
L <- L[L$detected %in% c(TRUE, "TRUE"), ]

ARM <- c(variant = "alcohol on variant", wild_type = "alcohol on wild type",
         extra_hit = "difference of the two", baseline = "variant vs wild type,\nno alcohol")

## ---- panel a : how many genes each contrast moves --------------------------
COLS <- c("variant", "wild_type", "extra_hit")
if ("DE_baseline" %in% names(S)) COLS <- c(COLS, "baseline")
cnt <- do.call(rbind, lapply(COLS, function(k) data.frame(
  cell = S$cell_type, arm = unname(ARM[k]),
  n = S[[paste0("DE_", k)]], stringsAsFactors = FALSE)))
cnt$arm  <- factor(cnt$arm, levels = unname(ARM[COLS]))
cnt$cell <- factor(cnt$cell, levels = S$cell_type[order(S$DE_wild_type)])

pa <- ggplot(cnt, aes(arm, cell, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n), size = NC_MIN / .pt, colour = NC_INK) +
  ## Counts span 0 to ~500, so a linear ramp would show one dark column and
  ## nothing else. sqrt keeps the small numbers legible without pretending
  ## the difference is smaller than it is; the printed numbers carry the truth.
  scale_fill_gradient(low = "white", high = NC_BLUE_L, trans = "sqrt",
                      breaks = c(0, 50, 150, 300, 500),
                      name = "genes at FDR < 0.05") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(30, "pt"),
                                barheight = grid::unit(5, "pt"),
                                title.position = "top")) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = NC_MIN, angle = 20, hjust = 0),
        axis.text.y = element_text(size = NC_MIN),
        legend.position = "bottom", legend.justification = "left",
        legend.title = element_text(size = NC_MIN),
        legend.text = element_text(size = NC_MIN)) +
  nc_tag("a")
save_nc(pa, file.path(OUT, "Explore_a_DE_counts_by_arm"), 112, 96)

## ---- panel b : the Cyps, both arms, one scale ------------------------------
HEP <- c("Periportal hepatocyte", "Mid-zone hepatocyte", "Pericentral hepatocyte",
         "ISG expressing hepatocyte", "Migratory_IR hepatocyte", "Remodeling hepatocyte")
SHORT <- c("periportal", "mid-zone", "pericentral", "ISG", "migratory/IR", "remodeling")
names(SHORT) <- HEP

d <- L[L$cell_type %in% HEP & L$contrast %in% c("variant", "wild_type"), ]
np <- tapply(d$symbol[d$contrast == "variant"], d$symbol[d$contrast == "variant"], length)
keep <- names(np)[np >= MIN_POP]
d <- d[d$symbol %in% keep, ]

## Row order by the variant-arm mean, so the panel reads top to bottom as the
## reported story does; the wild-type block is then directly comparable.
mv <- tapply(d$logFC[d$contrast == "variant"], d$symbol[d$contrast == "variant"],
             mean, na.rm = TRUE)
d$symbol    <- factor(d$symbol, levels = names(sort(mv, decreasing = TRUE)))
d$cell_type <- factor(unname(SHORT[d$cell_type]), levels = unname(SHORT[HEP]))
d$arm       <- factor(unname(ARM[d$contrast]), levels = unname(ARM[c("variant", "wild_type")]))
d$sig       <- !is.na(d$FDR) & d$FDR < FDR_CUT

LIM <- max(abs(d$logFC), na.rm = TRUE)
pb <- ggplot(d, aes(cell_type, symbol, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = d[d$sig, ], fill = NA, colour = NC_INK, linewidth = 0.5) +
  facet_wrap(~ arm, nrow = 1) +
  scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                       limits = c(-LIM, LIM), na.value = "#F2F2F2",
                       breaks = c(-3, -1.5, 0, 1.5, 3),
                       name = expression("log"[2]~"fold-change")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(30, "pt"),
                                barheight = grid::unit(5, "pt"),
                                title.position = "top")) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = NC_MIN, angle = 45, hjust = 0),
        axis.text.y = element_text(size = NC_MIN, face = "italic"),
        strip.text = element_text(size = NC_MIN, face = "bold"),
        panel.spacing = grid::unit(6, "pt"),
        legend.position = "bottom", legend.justification = "left",
        legend.title = element_text(size = NC_MIN),
        legend.text = element_text(size = NC_MIN)) +
  nc_tag("b")
save_nc(pb, file.path(OUT, "Explore_b_cyps_both_arms"), 150, 112)

## ---- what the panels say ---------------------------------------------------
cat("\n  DE genes at FDR < 0.05, summed over cell types\n")
cat(sprintf("    alcohol on variant    : %5d\n", sum(S$DE_variant)))
cat(sprintf("    alcohol on wild type  : %5d\n", sum(S$DE_wild_type)))
cat(sprintf("    difference of the two : %5d\n", sum(S$DE_extra_hit)))
cat(sprintf("    cell types where wild type > variant : %d of %d\n",
            sum(S$DE_wild_type > S$DE_variant), nrow(S)))
mw <- tapply(d$logFC[d$arm == ARM["wild_type"]], d$symbol[d$arm == ARM["wild_type"]],
             mean, na.rm = TRUE)
mv2 <- mv[names(mw)]
cat(sprintf("\n  Cyps in hepatocyte populations (n = %d drawn)\n", length(mw)))
cat(sprintf("    mean over genes, variant arm   : %+.2f\n", mean(mv2, na.rm = TRUE)))
cat(sprintf("    mean over genes, wild-type arm : %+.2f\n", mean(mw,  na.rm = TRUE)))
cat(sprintf("    more decreased on the variant  : %d of %d genes\n",
            sum(mv2 < mw, na.rm = TRUE), sum(!is.na(mv2 - mw))))
