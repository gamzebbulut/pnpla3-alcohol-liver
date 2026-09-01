# =============================================================================
# Figure 3 — LSEA, 2 rows x 4 columns  (C1 | C2 | C3 | C4  =  Q3 | Q4 | Q2 | Q5)
#   Row 1 (unsaturation): YOUR boxplots, unchanged look (plot_enrichment).
#   Row 2 (chain length): BARS = mean log2FC per bin, RED when the bin is
#                         enriched (FDR < SIG_FDR). x-labels rotated + sized
#                         so they no longer overlap.
#
# Significance + binning for the bars are taken from lipidr's own
# plot_enrichment() output, so the red in row 2 matches the red in row 1.
# Reuses `de`/`enr`/`short`/plot_enrichment()/significant_lipidsets() exactly
# like fig2_mountains_2x3.R. Run after lipid_LSEA_Q3_vs_Q4.R, or it sources it.
# =============================================================================

## ------------------------------- CONFIG -------------------------------------
LSEA_SCRIPT <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 3/lipid_LSEA_Q3_vs_Q4.R")
OUT_DIR <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 3/outs")   # always set (avoid stale OUT_DIR from a prior script in the same session)
if (!exists("SIG_FDR")) SIG_FDR <- 0.10

CONTRASTS   <- c("Q3", "Q4", "Q2", "Q5")   # display order C1, C2, C3, C4 (keys stay Q for data lookup)
FONT        <- "Arial"     # set "" if you get a font warning
BASE        <- 14
XTICK       <- 7           # row-2 x-axis tick font (small -> no overlap)
LABEL_EVERY <- 1           # label every Nth chain-length bin (2 = every other)
COL_SIG     <- "#D32F2F"; COL_NS <- "grey80"
FIG_W       <- 17; FIG_H <- 8.5
VERBOSE     <- TRUE        # print which column drives the red + which bins are sig

OUT_PNG <- file.path(OUT_DIR, "fig2_2x4_bars_Q2_Q4_Q3_Q5.png")
OUT_SVG <- file.path(OUT_DIR, "fig2_2x4_bars_Q2_Q4_Q3_Q5.svg")

## ----------------------------- dependencies ---------------------------------
if (!requireNamespace("patchwork", quietly = TRUE))
  stop("install.packages('patchwork')")
suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork) }))

if (!exists("de") || !exists("enr")) {
  message("`de`/`enr` not found -- sourcing LSEA script:\n  ", LSEA_SCRIPT)
  if (!file.exists(LSEA_SCRIPT)) stop("Set LSEA_SCRIPT to the real path.")
  source(LSEA_SCRIPT)
}
if (!exists("short"))
  short <- c("Q2 (mutation | no-alcohol)" = "Q2", "Q3 (alcohol | mutant bg)" = "Q3",
             "Q4 (alcohol | WT bg)" = "Q4", "Q5 (mutation | alcohol bg)" = "Q5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------- styling ------------------------------------
polished <- theme_bw(base_size = BASE) +
  theme(text             = element_text(family = FONT),
        plot.title       = element_text(size = BASE + 3, face = "bold", hjust = 0.5),
        axis.title       = element_text(size = BASE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text       = element_blank(),
        strip.background = element_blank(),
        legend.position  = "bottom")

## --------- hide the unparsed "NA" bin from the DISPLAY (not the stats) -------
drop_na_bin <- function(p) {
  xv <- tryCatch(rlang::as_label(rlang::quo_squash(p$mapping$x)), error = function(e) NA)
  if (is.na(xv) || !xv %in% names(p$data)) return(p)
  keep <- !is.na(p$data[[xv]]) & as.character(p$data[[xv]]) != "NA"
  p$data <- p$data[keep, , drop = FALSE]
  if (is.factor(p$data[[xv]])) p$data[[xv]] <- droplevels(p$data[[xv]])
  p
}

## --------- find which data column an aesthetic points at ---------------------
aes_col <- function(p, nm, d) {
  q <- p$mapping[[nm]]
  if (is.null(q)) return(NA_character_)
  v <- tryCatch(rlang::as_label(rlang::quo_squash(q)), error = function(e) NA_character_)
  if (!is.na(v) && v %in% names(d)) v else NA_character_
}
to_sig <- function(v) {
  if (is.logical(v)) return(v)
  vc <- tolower(trimws(as.character(v))); out <- rep(FALSE, length(vc))
  out[vc %in% c("true", "t", "1", "yes")] <- TRUE
  out[grepl("^sig|significant", vc) & !grepl("not|non|^ns$", vc)] <- TRUE
  out
}

## --------------------- the plot_enrichment for one contrast -----------------
enr_plot <- function(qcode, annot) {
  cn  <- names(short)[match(qcode, short)]
  sub <- de[de$contrast == cn, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  sig <- tryCatch(significant_lipidsets(enr[enr$contrast == cn, ], p.cutoff = SIG_FDR),
                  error = function(e) NULL)
  p <- tryCatch(plot_enrichment(sub, sig, annotation = annot),
                error = function(e) { message("plot_enrichment ", qcode, "/", annot,
                                               ": ", conditionMessage(e)); NULL })
  if (is.null(p)) NULL else drop_na_bin(p)
}
blank_panel <- function(title) ggplot() + theme_void() +
  labs(title = title) + theme(plot.title = element_text(face = "bold", hjust = 0.5))

## -------------------------- ROW 1: boxplots (unsat) -------------------------
box_panel <- function(qcode, title, xlab) {
  p <- enr_plot(qcode, "unsat")
  if (is.null(p)) return(blank_panel(title))
  p + labs(title = title, x = xlab, y = expression(log[2] ~ "fold-change")) +
    polished + theme(legend.position = "none")
}

## --------------------------- ROW 2: bars (length) ---------------------------
bar_panel <- function(qcode, title, xlab) {
  p <- enr_plot(qcode, "length")
  if (is.null(p)) return(blank_panel(title))
  d  <- p$data
  xc <- aes_col(p, "x", d); yc <- aes_col(p, "y", d)
  sc <- aes_col(p, "colour", d); if (is.na(sc)) sc <- aes_col(p, "color", d)
  if (is.na(sc)) sc <- aes_col(p, "fill", d)
  if (is.na(sc)) { cand <- intersect(c("significant","Significant","sig","is_significant","enriched"),
                                      names(d)); if (length(cand)) sc <- cand[1] }
  if (is.na(xc) || is.na(yc)) {
    message("bar ", qcode, ": x/y not found. names(p$data): ", paste(names(d), collapse = ", "))
    return(blank_panel(title))
  }
  d$.x <- as.character(d[[xc]]); d$.y <- suppressWarnings(as.numeric(d[[yc]]))
  d$.s <- if (!is.na(sc)) to_sig(d[[sc]]) else FALSE
  xlev <- if (is.factor(d[[xc]])) levels(droplevels(d[[xc]]))
          else { u <- unique(d$.x); u[order(suppressWarnings(as.numeric(u)))] }
  if (VERBOSE) message("bar ", qcode, " | x=", xc, " y=", yc,
                       " sig=", ifelse(is.na(sc), "<none>", sc),
                       " | sig bins: ", paste(sort(unique(d$.x[d$.s])), collapse = ", "))
  m  <- tapply(d$.y, d$.x, mean, na.rm = TRUE)
  sg <- tapply(d$.s, d$.x, function(v) any(v, na.rm = TRUE))
  agg <- data.frame(.x = factor(names(m), levels = xlev), mean = as.numeric(m),
                    stringsAsFactors = FALSE)
  agg$sig <- as.logical(sg[as.character(agg$.x)]); agg$sig[is.na(agg$sig)] <- FALSE
  agg$grp <- ifelse(agg$sig, "Significant", "Not significant")
  agg <- agg[order(agg$.x), ]
  brks <- xlev[seq(1, length(xlev), by = LABEL_EVERY)]
  ggplot(agg, aes(.x, mean, fill = grp)) +
    geom_col(width = 0.82, colour = "grey35", linewidth = 0.18) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
    scale_fill_manual(values = c("Significant" = COL_SIG, "Not significant" = COL_NS),
                      name = "Enrichment") +
    scale_x_discrete(breaks = brks, drop = FALSE) +
    labs(title = title, x = xlab, y = expression(log[2] ~ "fold-change")) +
    polished +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = XTICK))
}

## ------------------------------- build it -----------------------------------
xlab_unsat  <- "Total unsaturation"
xlab_length <- "Total chain length (carbons)"
long_titles <- c(Q3 = "C1 \u00b7 alcohol (mutant bg)", Q4 = "C2 \u00b7 alcohol (WT bg)",
                 Q2 = "C3 \u00b7 mutation",            Q5 = "C4 \u00b7 mutation (alcohol bg)")
c_short     <- c(Q3 = "C1", Q4 = "C2", Q2 = "C3", Q5 = "C4")   # row-2 short labels

row1 <- lapply(CONTRASTS, function(q) box_panel(q, long_titles[[q]], xlab_unsat))
row2 <- lapply(CONTRASTS, function(q) bar_panel(q, c_short[[q]], xlab_length))

top <- Reduce(`|`, row1)
bot <- Reduce(`|`, row2)
fig <- (top / bot) + plot_layout(guides = "collect") &
  theme(text = element_text(family = FONT))

## -------------------------------- save --------------------------------------
ggsave(OUT_PNG, fig, width = FIG_W, height = FIG_H, dpi = 300, bg = "white")
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(OUT_SVG, fig, width = FIG_W, height = FIG_H, device = svglite::svglite, bg = "white")
cat("\nSaved:\n  ", OUT_PNG, "\n",
    if (requireNamespace("svglite", quietly = TRUE)) paste0("  ", OUT_SVG, "\n") else "")
