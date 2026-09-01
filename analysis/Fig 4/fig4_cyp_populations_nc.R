# =============================================================================
#  Figure 4 — the Cyp panel as a population grid  (ALTERNATIVE to Fig4e)
#
#  WHY THIS EXISTS
#  The reported statistic is one number per gene, computed on one value per
#  mouse. That is the right unit, and it is also the reason the current panel e
#  now shows only three full-tone bars: at n = 5 vs 4 animals, a transcriptome-
#  wide FDR is a high bar and most of the programme sits just under it.
#
#  Steve's sanctioned route (edgeR glmQLFTest, one test per hepatocyte
#  population, mouse as the unit throughout — fig4_steve_edgeR_qlf.R) gives
#  something the single number cannot: the same contrast measured six times, in
#  six different hepatocyte populations. Every Cyp with a substantial effect
#  moves the same way in all six. That is worth showing, because it is the
#  strongest honest statement available about this programme.
#
#  WHAT THIS PANEL CLAIMS, AND WHAT IT MUST NOT
#  It claims PERVASIVENESS: the decrease is present in every hepatocyte
#  population, not carried by one. It does NOT claim replication. The six
#  populations come from the SAME nine animals, so six consistent columns are
#  one experiment seen six ways, not six experiments. Never quote a P value for
#  "same sign in 6 of 6" — the columns are not independent, which is the same
#  error that made the stacked model return 1e-13. The inferential claim stays
#  the animal-level one: 18 of 24 down, sign test P = 0.023.
#
#  ENCODING
#  Fill is the effect size, red increased / blue decreased, on one continuous
#  scale. Significance is a BORDER, not a saturation step: with three or four
#  significant cells in the grid, tinting everything else would read as missing
#  data. This is the same grammar proposed for Figure 5a — shade by effect,
#  outline for significance — so the two figures teach the reader one rule.
#
#  INPUT: fig4_outs/steve_edgeR_qlf_Cyps_Q3.csv, written by fig4_steve_edgeR_qlf.R
#         (which needs fig4_animal_level_limma.R to have run first).
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
QUESTION   <- "Q3"
FDR_CUT    <- 0.05
SHOW_SETS  <- TRUE     # the AA / LA / PPAR membership string in the right margin
W_MM       <- 108
H_MM       <- 96
## Populations in reading order, and the short labels the panel prints.
POP_LABEL  <- c("Periportal hepatocyte"     = "periportal",
                "Mid-zone hepatocyte"       = "mid-zone",
                "Pericentral hepatocyte"    = "pericentral",
                "ISG expressing hepatocyte" = "ISG",
                "Migratory_IR hepatocyte"   = "migratory/IR",
                "Remodeling hepatocyte"     = "remodeling")
SUMMARY_LAB <- "per mouse"   # the animal-level column, drawn first and ruled off
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
OUTS    <- Sys.getenv("FIG4_OUTS",  file.path(FIGROOT, "Fig 4/fig4_outs"))
OUT     <- Sys.getenv("FIG4_OUT",   file.path(FIGROOT, "Fig 4/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

f <- file.path(OUTS, sprintf("steve_edgeR_qlf_Cyps_%s.csv", QUESTION))
if (!file.exists(f)) stop("Run fig4_steve_edgeR_qlf.R first -- missing ", f)
z <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

pops <- sub("^logFC\\.", "", grep("^logFC\\.", names(z), value = TRUE))
pops <- names(POP_LABEL)[names(POP_LABEL) %in% pops]     # reading order
if (!length(pops)) stop("No per-population columns found in ", basename(f))

## ---- long form: one row per gene x column ----------------------------------
mk <- function(col_fc, col_q, label) data.frame(
  symbol = z$symbol, sets = z$sets, col = label,
  fc = z[[col_fc]], q = z[[col_q]], stringsAsFactors = FALSE)

D <- rbind(mk("logFC_animal", "FDR_animal", SUMMARY_LAB),
           do.call(rbind, lapply(pops, function(p)
             mk(paste0("logFC.", p), paste0("FDR.", p), unname(POP_LABEL[p])))))
D$sig <- !is.na(D$q) & D$q < FDR_CUT

## Rows ordered by the reported effect, most decreased at the top. ggplot puts
## the first level at the bottom, so the levels run the other way.
ord <- z$symbol[order(z$logFC_animal, decreasing = TRUE)]
D$symbol <- factor(D$symbol, levels = ord)
D$col    <- factor(D$col, levels = c(SUMMARY_LAB, unname(POP_LABEL[pops])))

LIM <- max(abs(D$fc), na.rm = TRUE)

## ---- consistency, for the console and the legend ---------------------------
same6 <- vapply(z$symbol, function(s) {
  v <- D$fc[D$symbol == s & D$col != SUMMARY_LAB]
  v <- v[!is.na(v)]
  length(v) > 0 && (all(v < 0) || all(v > 0)) }, logical(1))
down6 <- vapply(z$symbol, function(s) {
  v <- D$fc[D$symbol == s & D$col != SUMMARY_LAB]; v <- v[!is.na(v)]
  length(v) > 0 && all(v < 0) }, logical(1))

## ---- the panel -------------------------------------------------------------
NCOL <- nlevels(D$col)
p <- ggplot(D, aes(col, symbol)) +
  geom_tile(aes(fill = fc), colour = "white", linewidth = 0.4) +
  ## Significance is an outline drawn ON TOP, so it survives any fill.
  geom_tile(data = D[D$sig, ], fill = NA, colour = NC_INK, linewidth = 0.55) +
  ## Rule between the reported number and the six populations that support it.
  ## A segment, not geom_vline: with clip = "off" a vline runs past the grid
  ## and into the legend.
  annotate("segment", x = 1.5, xend = 1.5,
           y = 0.5, yend = length(unique(D$symbol)) + 0.5,
           linewidth = 0.4, colour = NC_INK) +
  scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                       limits = c(-LIM, LIM), na.value = "#F2F2F2",
                       name = expression("Log"[2]~"fold-change"),
                       breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(28, "pt"),
                                barheight = grid::unit(5, "pt"),
                                title.position = "top", ticks.colour = NC_INK)) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = NC_MIN, angle = 45, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = NC_MIN, face = "italic"),
        legend.position = "bottom", legend.justification = "left",
        legend.title = element_text(size = NC_MIN),
        legend.text = element_text(size = NC_MIN),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "pt")) +
  nc_tag("e")

if (SHOW_SETS) {
  sets <- data.frame(symbol = factor(z$symbol, levels = ord),
                     col = NCOL + 0.75, lab = z$sets, stringsAsFactors = FALSE)
  p <- p + geom_text(data = sets, aes(x = col, y = symbol, label = lab),
                     inherit.aes = FALSE, hjust = 0, size = NC_MIN / .pt,
                     colour = NC_INK2) +
    coord_cartesian(xlim = c(0.5, NCOL + 2.9), clip = "off")
}

save_nc(p, file.path(OUT, "Fig4e_cyp_population_grid"), W_MM, H_MM)

## ---- report ----------------------------------------------------------------
cat(sprintf("\n  %d Cyps x %d populations, plus the reported per-mouse column\n",
            nrow(z), length(pops)))
cat(sprintf("  same direction in all %d populations : %d of %d genes\n",
            length(pops), sum(same6), nrow(z)))
cat(sprintf("  DOWN in all %d populations           : %d of %d genes\n",
            length(pops), sum(down6), nrow(z)))
cat(sprintf("  cells outlined (FDR < %.2f)          : %d of %d\n",
            FDR_CUT, sum(D$sig), nrow(D)))
cat("\n  down in every population:\n     ")
cat(paste(z$symbol[down6], collapse = ", "), "\n")
cat("\n  CAPTION MATERIAL, worded so it claims only what it can:\n")
cat("  \"the decrease is present in every hepatocyte population\" -- pervasiveness.\n")
cat("  NOT \"confirmed in six independent tests\": same nine animals, six views.\n")
