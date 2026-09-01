# =============================================================================
# Figure 2 (Banini-style) — per-species arrow grid
#   Mirrors Banini 2021 Fig 3D: each lipid SPECIES is one arrow plotted at
#   (acyl-chain double bonds  x ,  carbon number  y). The arrow points up/down
#   by logFC sign; RED = that species is significant, GREY = not.
#   NO enrichment here, by design -- this is the per-species view in Dr Sanyal's
#   own figure. Keep your mountains script as the enrichment companion.
#
# Reuses the SAME `de` object as fig2_mountains_2x3.R. Run AFTER
# lipid_LSEA_Q3_vs_Q4.R (so `de` exists) or it will source it for you.
# =============================================================================

## ------------------------------- CONFIG -------------------------------------
LSEA_SCRIPT <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 3/lipid_LSEA_Q3_vs_Q4.R")
OUT_DIR <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 3/outs")   # always set (avoid stale OUT_DIR from a prior script in the same session)

CONTRASTS  <- c("Q3", "Q4", "Q2", "Q5")   # display order C1,C2,C3,C4 (keys stay Q); drop last to go back to 3 panels
CLASS_KEEP <- "TG"                  # Banini match = TG only.  NULL = ALL classes.
SIG_USE    <- "adj"                 # "adj" = adj.P.Val (BH FDR) -> matches A/B FDR<0.10 |
                                   # "adj" = adj.P.Val (strict FDR) | "col" = de$significant
SIG_P      <- 0.10                  # FDR cutoff for the arrow colour (matches LSEA FDR<0.10 in A/B)
COL_SIG    <- "#D32F2F"; COL_NS <- "#BDBDBD"
JITTER_W   <- 0                     # 0 = sit on the gridlines (TG: 1 species/cell).
                                   # bump to ~0.12 only if you set CLASS_KEEP <- NULL
FONT       <- "Arial"; BASE <- 14   # set FONT <- "" if you get a font warning
FIG_W      <- 15.5; FIG_H <- 4.6    # inches (one row of panels)
Y_RANGE    <- NULL                  # NULL = your data's carbon range.
                                   # set c(31, 61) to match Banini's axis exactly

OUT_PNG <- file.path(OUT_DIR, "fig2_banini_arrows.png")
OUT_SVG <- file.path(OUT_DIR, "fig2_banini_arrows.svg")

suppressWarnings(suppressMessages({ library(ggplot2) }))

## --------------------- make sure `de` is available --------------------------
if (!exists("de")) {
  message("`de` not found -- sourcing LSEA script:\n  ", LSEA_SCRIPT)
  if (!file.exists(LSEA_SCRIPT)) stop("Set LSEA_SCRIPT to the real path.")
  source(LSEA_SCRIPT)
}
if (!exists("short"))
  short <- c("Q2 (mutation | no-alcohol)" = "Q2", "Q3 (alcohol | mutant bg)" = "Q3",
             "Q4 (alcohol | WT bg)" = "Q4", "Q5 (mutation | alcohol bg)" = "Q5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## --------------- find the columns robustly (lipidr naming varies) -----------
pick <- function(df, cands) { h <- cands[cands %in% names(df)]; if (length(h)) h[1] else NA }
name_col <- pick(de, c("Molecule", "clean_name", "INDIVIDUAL.NOTATION", "variable", "Lipid"))
lfc_col  <- pick(de, c("logFC", "log2FC", "log2FoldChange", "estimate"))
adj_col  <- pick(de, c("adj.P.Val", "padj", "qvalue", "q.value", "FDR"))
raw_col  <- pick(de, c("P.Value", "pvalue", "p.value", "PValue", "p"))
sig_col  <- pick(de, c("significant"))
cls_col  <- pick(de, c("Class", "class"))
if (is.na(name_col) || is.na(lfc_col))
  stop("Couldn't find a molecule-name or logFC column in `de`. names(de): ",
       paste(names(de), collapse = ", "))
if (!is.null(CLASS_KEEP) && is.na(cls_col))
  message("No Class column found -- CLASS_KEEP ignored, plotting ALL species.")

## -------- parse TOTAL carbons & TOTAL double bonds from the species name -----
# e.g. "TG(52:5)" -> 52 carbons, 5 db ;  "PC(16:0/18:1)" -> 34 carbons, 1 db
parse_cl_cs <- function(x) {
  hits <- regmatches(x, gregexpr("([0-9]+):([0-9]+)", x))
  cl <- cs <- rep(NA_real_, length(x))
  for (i in seq_along(hits)) {
    h <- hits[[i]]; if (!length(h)) next
    nums <- do.call(rbind, lapply(strsplit(h, ":"), as.numeric))
    cl[i] <- sum(nums[, 1]); cs[i] <- sum(nums[, 2])
  }
  data.frame(cl = cl, cs = cs)
}

## ------------------------- build the plotting frame -------------------------
build_df <- function(qcode) {
  cn  <- names(short)[match(qcode, short)]
  sub <- de[de$contrast == cn, , drop = FALSE]
  if (!nrow(sub)) { message("No rows for ", qcode); return(NULL) }
  if (!is.null(CLASS_KEEP) && !is.na(cls_col))
    sub <- sub[sub[[cls_col]] %in% CLASS_KEEP, , drop = FALSE]
  if (!nrow(sub)) { message("No kept species for ", qcode); return(NULL) }
  cc  <- parse_cl_cs(as.character(sub[[name_col]]))
  lfc <- suppressWarnings(as.numeric(sub[[lfc_col]]))
  pcol <- switch(SIG_USE,
                 raw = raw_col, adj = adj_col,
                 col = NA, stop("SIG_USE must be 'raw', 'adj', or 'col'"))
  if (SIG_USE == "col") {
    if (is.na(sig_col)) stop("SIG_USE='col' but no `significant` column in `de`.")
    sig <- as.logical(sub[[sig_col]])
  } else {
    if (is.na(pcol)) {                                  # asked-for p not found -> fall back
      alt <- if (SIG_USE == "raw") adj_col else raw_col
      if (is.na(alt)) stop("No p-value column found in `de`.")
      message("SIG_USE='", SIG_USE, "' column missing; using '", alt, "' instead.")
      pcol <- alt
    }
    sig <- suppressWarnings(as.numeric(sub[[pcol]])) < SIG_P
  }
  sig[is.na(sig)] <- FALSE
  out <- data.frame(contrast = qcode, cl = cc$cl, cs = cc$cs, lfc = lfc, sig = sig)
  out[!is.na(out$cl) & !is.na(out$cs) & !is.na(out$lfc), , drop = FALSE]
}

df <- do.call(rbind, lapply(CONTRASTS, build_df))
if (is.null(df) || !nrow(df)) stop("No species to plot after parsing/filtering.")
df$contrast <- factor(df$contrast, levels = CONTRASTS)
df$dir  <- ifelse(df$lfc > 0, "up", "down")
df$sigf <- ifelse(df$sig, "Significant", "Not significant")

## --------------------------- Banini-style arrow grid ------------------------
ylo <- if (is.null(Y_RANGE)) floor(min(df$cl) / 2) * 2 else Y_RANGE[1]
yhi <- if (is.null(Y_RANGE)) ceiling(max(df$cl) / 2) * 2 else Y_RANGE[2]
yb  <- seq(ylo, yhi, by = 2)
xb <- seq(0, max(df$cs))
p <- ggplot(df, aes(cs, cl)) +
  geom_jitter(aes(shape = dir, fill = sigf, colour = sigf),
              width = JITTER_W, height = 0, size = 2.6, stroke = 0.3) +
  scale_shape_manual(values = c(up = 24, down = 25), name = NULL,
                     labels = c(down = "Decreased", up = "Increased")) +
  scale_fill_manual(values = c("Significant" = COL_SIG, "Not significant" = COL_NS), name = NULL) +
  scale_colour_manual(values = c("Significant" = COL_SIG, "Not significant" = "#9E9E9E"),
                      guide = "none") +
  scale_x_continuous(breaks = xb) +
  scale_y_continuous(breaks = yb) +
  facet_wrap(~contrast, nrow = 1,
             labeller = as_labeller(c(Q3 = "C1 \u00b7 alcohol (mutant bg)",
                                      Q4 = "C2 \u00b7 alcohol (WT bg)",
                                      Q2 = "C3 \u00b7 mutation",
                                      Q5 = "C4 \u00b7 mutation (alcohol bg)"))) +
  labs(x = "Acyl-chain double bonds", y = "Carbon number") +
  theme_classic(base_size = BASE) +
  theme(text             = element_text(family = FONT),
        strip.text       = element_text(face = "bold", size = BASE + 1),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom") +
  guides(fill  = guide_legend(override.aes = list(shape = 24, size = 3.5)),
         shape = guide_legend(override.aes = list(fill = "#777777", size = 3.5)))

if (!is.null(Y_RANGE)) p <- p + coord_cartesian(ylim = Y_RANGE)

cat("\nSignificant species per contrast (SIG_USE=", SIG_USE, ", p<", SIG_P, "):\n", sep = "")
print(tapply(df$sig, df$contrast, sum))
cat("Total species per contrast:\n"); print(table(df$contrast))

ggsave(OUT_PNG, p, width = FIG_W, height = FIG_H, dpi = 300, bg = "white")
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(OUT_SVG, p, width = FIG_W, height = FIG_H, device = svglite::svglite, bg = "white")

cat("\nSaved:\n  ", OUT_PNG, "\n",
    if (requireNamespace("svglite", quietly = TRUE)) paste0("  ", OUT_SVG, "\n") else "")
