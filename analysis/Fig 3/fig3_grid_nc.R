# =============================================================================
#  Figure 3, panels a-c — Nature Communications style
#    a  species log2 fold-change by total acyl-chain unsaturation
#    b  mean species log2 fold-change by total chain length
#    c  per-species view of triglycerides: chain length vs unsaturation
#  Each panel shows all four comparisons side by side.
#
#  Replaces fig3_2x4_bars.R and fig3_banini_arrows.R.
#
#  CHANGED:
#    - Q2/Q3/Q4/Q5 and C1-C4 gone; comparisons named in plain language (task A2),
#      ordered alcohol first (task A1).
#    - Y AXES ARE NOW SHARED WITHIN EACH ROW. They were free before, so the four
#      panels of a ran -2..1.7, -1.6..2.4, -1.5..1.5 and -1..1.3: a bar that
#      looked tall in one panel was smaller than a short-looking bar in another,
#      even though the figure exists to compare them.
#    - COLOUR CARRIES WHATEVER POSITION DOES NOT ALREADY SHOW.
#      In a and b, direction is already visible: a bar above zero is increased,
#      a box above the dashed line is up. Colouring by direction there would
#      repeat the axis and leave significance to a 7 pt asterisk. So in a and b
#      red marks the enriched sets, as before.
#      In c, position encodes chain length and unsaturation, and nothing shows
#      direction but the triangle's orientation at 1.5 pt. So colour carries
#      direction there: red increased, blue decreased.
#    - Panel c drops its significance encoding. Two of 264 triglyceride species
#      reach species-level significance; that belongs in one sentence of the
#      legend, not in a colour channel.
#      NOTE the cost of this choice: red means "enriched" in a and b and
#      "increased" in c. The legend must say so explicitly.
#    - Computed straight from the two result CSVs rather than by reaching into
#      lipidr's ggplot objects, which broke whenever lipidr changed its columns.
#    - 180 mm rows, 5-7 pt type, vector SVG.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("FIG3_IN",    file.path(FIGROOT, "Fig 3/outs"))
OUT     <- Sys.getenv("FIG3_OUT",   file.path(FIGROOT, "Fig 3/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

FDR <- 0.10

## ---- provenance ------------------------------------------------------------
## This script does NOT run the enrichment. It reads the stored results, so its
## output is deterministic: the same bins are marked every time. The permutation
## step lives upstream in lipid_LSEA_Q3_vs_Q4.R, where set.seed(1) sits at the
## top of the file rather than beside the lsea() call ~127 lines below it. That
## is reproducible for a clean source(), but not if that script is run in pieces
## interactively or after other code has drawn from the RNG. Deliberately left
## as is; the fingerprint below is what makes any drift visible.
prov <- function(f) {
  p <- file.path(IN, f)
  sprintf("  %-26s %s  md5 %s", f,
          format(file.info(p)$mtime, "%Y-%m-%d %H:%M"),
          tools::md5sum(p)[[1]])
}
cat("\n-- input fingerprint --\n", prov("de_lipids_Q2_Q5.csv"), "\n",
    prov("lsea_results_Q2_Q5.csv"), "\n", sep = "")

de <- utils::read.csv(file.path(IN, "de_lipids_Q2_Q5.csv"), stringsAsFactors = FALSE)
en <- utils::read.csv(file.path(IN, "lsea_results_Q2_Q5.csv"), stringsAsFactors = FALSE)

## plain language, alcohol first
MAP <- c("Q4 (alcohol | WT bg)"       = "Alcohol on wild type",
         "Q3 (alcohol | mutant bg)"   = "Alcohol on variant",
         "Q2 (mutation | no-alcohol)" = "Variant alone",
         "Q5 (mutation | alcohol bg)" = "Variant on alcohol")
ORDER <- unname(MAP)
de$cmp <- factor(unname(MAP[de$contrast]), levels = ORDER)
en$cmp <- factor(unname(MAP[en$contrast]), levels = ORDER)
de <- de[!is.na(de$cmp), ]; en <- en[!is.na(en$cmp), ]

## enriched sets, by family
enriched <- function(prefix) {
  e <- en[startsWith(en$set, prefix) & !is.na(en$padj) & en$padj < FDR, ]
  paste(e$cmp, sub(prefix, "", e$set))
}
SIG_CS <- enriched("total_cs_"); SIG_CL <- enriched("total_cl_")

DIRCOL <- c(increased = NC_RED, decreased = NC_BLUE)          # panel c only
ENRCOL <- c("Enriched (FDR < 0.10)" = NC_RED, "Not enriched" = "#9E9E9E")

facets <- function() facet_wrap(~ cmp, nrow = 1)
legend_top <- function() theme(legend.position = "top",
                               legend.key.size = grid::unit(6, "pt"),
                               strip.text = element_text(size = NC_BASE))

## ---- a : distribution of species fold-changes by unsaturation --------------
A <- de[!is.na(de$total_cs), ]
A$bin <- factor(A$total_cs, levels = sort(unique(A$total_cs)))
A$enr <- ifelse(paste(A$cmp, A$total_cs) %in% SIG_CS,
                "Enriched (FDR < 0.10)", "Not enriched")

a <- ggplot(A, aes(bin, logFC, colour = enr)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "#888888") +
  geom_boxplot(outlier.size = 0.25, outlier.colour = "#666666",
               linewidth = 0.3, fatten = 1.1, fill = NA) +
  scale_colour_manual(values = ENRCOL, name = NULL, drop = FALSE,
                      breaks = names(ENRCOL)) +
  facets() +
  labs(x = "Total acyl-chain double bonds", y = expression("Log"[2]~"fold-change")) +
  theme_nc(grid_y = TRUE) + legend_top() +
  nc_tag("a")
save_nc(a, file.path(OUT, "Fig3a_unsaturation"), width_mm = 180, height_mm = 52)

## ---- b : mean fold-change by chain length ----------------------------------
B <- de[!is.na(de$total_cl), ]
agg <- do.call(rbind, lapply(split(B, list(B$cmp, B$total_cl), drop = TRUE), function(z)
  data.frame(cmp = z$cmp[1], cl = z$total_cl[1], mean = mean(z$logFC, na.rm = TRUE),
             n = nrow(z), stringsAsFactors = FALSE)))
agg$cmp <- factor(agg$cmp, levels = ORDER)
agg$enr <- ifelse(paste(agg$cmp, agg$cl) %in% SIG_CL,
                  "Enriched (FDR < 0.10)", "Not enriched")
brk <- sort(unique(agg$cl)); brk <- brk[seq(1, length(brk), by = 4)]

b <- ggplot(agg, aes(cl, mean, fill = enr)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NC_INK) +
  geom_col(width = 1.6, colour = NC_INK, linewidth = 0.08) +
  scale_fill_manual(values = ENRCOL, name = NULL, drop = FALSE,
                    breaks = names(ENRCOL)) +
  scale_x_continuous(breaks = brk) +
  facets() +
  labs(x = "Total chain length (carbons)", y = expression("Mean log"[2]~"fold-change")) +
  theme_nc(grid_y = TRUE) + legend_top() +
  nc_tag("b")
save_nc(b, file.path(OUT, "Fig3b_chainlength"), width_mm = 180, height_mm = 52)

## ---- c : per-species triglycerides ----------------------------------------
C <- de[de$Class == "TG" & !is.na(de$total_cl) & !is.na(de$total_cs), ]
C$dir <- ifelse(C$logFC > 0, "increased", "decreased")

cc <- ggplot(C, aes(total_cs, total_cl, shape = dir, fill = dir)) +
  geom_point(size = 1.5, colour = NC_INK, stroke = 0.15) +
  ## The data values stay lower case -- only what the reader sees is relabelled,
  ## so nothing that keys off C$dir has to change with it.
  scale_shape_manual(values = c(increased = 24, decreased = 25), name = NULL,
                     breaks = c("increased", "decreased"),
                     labels = c(increased = "Increased", decreased = "Decreased")) +
  scale_fill_manual(values = DIRCOL, name = NULL,
                    breaks = c("increased", "decreased"),
                    labels = c(increased = "Increased", decreased = "Decreased")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  facets() +
  labs(x = "Acyl-chain double bonds", y = "Carbon number") +
  guides(shape = guide_legend(order = 1), fill = guide_legend(order = 1)) +
  theme_nc(grid_y = TRUE) + legend_top() +
  nc_tag("c")
save_nc(cc, file.path(OUT, "Fig3c_TG_species"), width_mm = 180, height_mm = 58)

## ---- Source Data + what the panels contain ---------------------------------
utils::write.csv(de[, c("cmp","Molecule","Class","total_cl","total_cs",
                        "logFC","P.Value","adj.P.Val")],
                 file.path(OUT, "SourceData_Fig3_species.csv"), row.names = FALSE)
utils::write.csv(en[, c("cmp","set","NES","pval","padj","size")],
                 file.path(OUT, "SourceData_Fig3_enrichment.csv"), row.names = FALSE)

cat("\n-- enriched sets at FDR < 0.10 --\n")
for (lv in ORDER) {
  s <- en[en$cmp == lv & !is.na(en$padj) & en$padj < FDR, ]
  cat(sprintf("  %-22s %2d sets: %s\n", lv, nrow(s),
              if (nrow(s)) paste(s$set, collapse = ", ") else "none"))
}
cat(sprintf("\n  panel c: %d triglyceride species; %d reach species-level FDR < 0.10\n",
            nrow(C), sum(!is.na(C$adj.P.Val) & C$adj.P.Val < FDR)))
cat("           (stated in the legend, no longer encoded in the artwork)\n")
