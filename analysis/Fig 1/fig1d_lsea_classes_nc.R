# =============================================================================
#  Figure 1 panel d — lipid classes, set-level, both alcohol contrasts
#
#  WHY THIS IS NOT PANEL c AGAIN
#  Panel c counts how many species in a class changed. This asks a different
#  question: did the class move as a SET. The two differ most exactly where
#  panel c is weakest -- its glycerophospholipid bar is 48 of the 102 species
#  and pools PC, PE, LPC, LPE, LPG and LPI into one row. LSEA separates them,
#  and they do not behave alike: phosphatidylcholines and the lysophospholipids
#  go up as sets while triglycerides go down. A count cannot show that, and a
#  reader looking at panel c would reasonably assume the 48 were one thing.
#
#  CLASS SETS ONLY. The LSEA also carries chain-length and unsaturation sets.
#  Those are Figure 3's subject and drawing them here would put the same result
#  in two figures. CHAIN_AND_UNSAT flips them on if that changes.
#
#  ORIENTATION matches Fig 4e: a wide short strip, categories across x, the two
#  contrasts as rows. Two columns of a tall heatmap would waste a page width.
#
#  WHAT IT MAY NOT SAY
#  Two rows side by side invite "this class responds only on the wild type".
#  Phosphatidylcholine is the trap here: NES +2.55 (q 0.001) on the wild type
#  against +1.28 (q 0.25) on the variant, which is a difference in significance
#  and not a demonstrated difference. The honest reading is the shared one --
#  the same classes move the same way on both backgrounds, triglycerides down
#  and phosphatidylcholines and lysophospholipids up -- and the script prints
#  which classes agree so the legend can say it in numbers.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
SIG            <- 0.10     # LSEA padj; nothing in the variant arm clears 0.05
CHAIN_AND_UNSAT <- FALSE   # TRUE also draws the total_cl_* / total_cs_* sets
ARMS           <- c(Q4 = "Alcohol on wild type", Q3 = "Alcohol on variant")
## Half width, sharing the bottom row with panel c. The class abbreviations
## on the x axis are short (CAR, TG, PC ...), so 12 columns fit 88 mm
## comfortably; it was never this panel that needed the page width.
W_MM           <- 88
H_MM           <- 50   # matched to BAR_H_MM in fig1_volcano_bars_nc.R
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("FIG1_LSEA", file.path(FIGROOT, "Fig 3/Lipid_LSEA_Q3_vs_Q4",
                                             "lsea_ALLQ_comparison.csv"))
OUT     <- Sys.getenv("FIG1_OUT",  file.path(FIGROOT, "Fig 1/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

if (!file.exists(IN)) stop("LSEA table not found: ", IN)
z <- utils::read.csv(IN, stringsAsFactors = FALSE)

keep <- if (CHAIN_AND_UNSAT) rep(TRUE, nrow(z)) else grepl("^Class_", z$set)
z <- z[keep, ]
stopifnot("No sets left after filtering" = nrow(z) > 0)

D <- do.call(rbind, lapply(names(ARMS), function(q) data.frame(
  set  = sub("^Class_", "", z$set),
  arm  = unname(ARMS[q]),
  NES  = z[[paste0("NES_",  q)]],
  padj = z[[paste0("padj_", q)]],
  stringsAsFactors = FALSE)))
D <- D[!is.na(D$NES), ]
D$sig <- !is.na(D$padj) & D$padj < SIG

## Order by the variant arm, the same arm panel c is ordered by, so the two
## panels stack in a consistent order rather than each with its own.
ord <- D[D$arm == unname(ARMS["Q3"]), ]
LEV <- ord$set[order(ord$NES)]
D$set <- factor(D$set, levels = LEV)
D$arm <- factor(D$arm, levels = rev(unname(ARMS)))   # first level plots at bottom

LIM <- max(abs(D$NES), na.rm = TRUE)

p <- ggplot(D, aes(set, arm)) +
  geom_tile(aes(fill = NES), colour = "white", linewidth = 0.4) +
  geom_tile(data = D[D$sig, ], fill = NA, colour = NC_INK, linewidth = 0.55) +
  scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                       limits = c(-LIM, LIM), na.value = "#F2F2F2",
                       name = "Normalized enrichment score",
                       breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(expand = expansion(add = 0.04)) +
  ## Reader-facing arm names from nc_style.R's display layer (nc_label_expr),
  ## parsed for the I148M superscript. The data keeps the internal keys.
  scale_y_discrete(expand = expansion(add = 0.04), labels = nc_label_expr) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(34, "pt"),
                                barheight = grid::unit(5.5, "pt"),
                                title.position = "left", ticks.colour = NC_INK)) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        ## Rotated: at 88 mm the twelve class abbreviations are about 5 mm
        ## apart and run together horizontally ("CARTG", "PCLPCLPILPG").
        axis.text.x = element_text(size = NC_BIG, angle = 90, hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(size = NC_BIG),
        legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = NC_BIG),
        legend.text = element_text(size = NC_BIG),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "pt")) +
  nc_tag("d")

save_nc(p, file.path(OUT, "Fig1d_lsea_classes"), W_MM, H_MM)

## ---- what it says ----------------------------------------------------------
w <- D[D$arm == unname(ARMS["Q4"]), ]; v <- D[D$arm == unname(ARMS["Q3"]), ]
m <- merge(w[, c("set","NES","padj")], v[, c("set","NES","padj")], by = "set",
           suffixes = c("_wt", "_var"))
m <- m[order(m$NES_var), ]
cat(sprintf("\n  %-8s %9s %8s %9s %8s %8s\n", "class", "NES wt", "q wt", "NES var", "q var", "agree"))
for (i in seq_len(nrow(m)))
  cat(sprintf("  %-8s %9.2f %8.2f %9.2f %8.2f %8s\n", m$set[i], m$NES_wt[i],
              m$padj_wt[i], m$NES_var[i], m$padj_var[i],
              if (sign(m$NES_wt[i]) == sign(m$NES_var[i])) "yes" else "-"))
ag <- sign(m$NES_wt) == sign(m$NES_var)
cat(sprintf("\n  same direction in both arms : %d of %d classes\n", sum(ag), nrow(m)))
cat(sprintf("  significant (q < %.2f): %d on wild type, %d on variant\n", SIG,
            sum(w$sig), sum(v$sig)))
cat("\n  LEGEND MATERIAL: say the shared direction, not the difference. A class\n")
cat("  significant in one arm and not the other has NOT been shown to differ\n")
cat("  between arms -- phosphatidylcholine is the one that will tempt you.\n")
