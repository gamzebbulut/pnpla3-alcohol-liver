# =============================================================================
# Figure 2 (A–C) — C1-significant metabolite heatmaps, HORIZONTAL, 3 class panels.
# RELABELED + REORDERED version of F1/Fig1_metabolite_panels_3groups.R:
#   contrasts now read top->bottom  C1, C2, C3, C4  (was Q2/Q3/Q4/Q5)
#     C1 = Q3 = alcohol on PNPLA3-I148M (mutant) background  (LEAD)
#     C2 = Q4 = alcohol on WT background
#     C3 = Q2 = mutation alone (no alcohol)
#     C4 = Q5 = mutation on alcohol background
# Only labels + row order change; the data and significance are identical
# (anchor set = metabolites significant in C1 / old Q3; reads the same FC_*.csv).
#
# Reads FC_Q2/Q3/Q4/Q5_metabolites.csv from IN_DIR. Run:
#   source("fig2_metabolite_heatmaps_C1C4.R")
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })

## ---- paths + settings (edit here) ------------------------------------------
IN_DIR   <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/F1"  # where FC_*.csv live
OUT_DIR  <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/F2"
P_CUTOFF <- 0.05
Q_CUTOFF <- 0.10
PER_W    <- 0.26          # inches per metabolite column (same tile size in all panels)
ORDER_BY <- "log2FC"      # "log2FC" (gradient) or "class"
FONT <- "Arial"; RED <- "#C62828"; BLUE <- "#1565C0"
RELABEL <- c(
  "Ursodeoxycholic + Hyodeoxycholic acid" = "UDCA + HDCA",
  "Flavin adenine dinucleotide"           = "FAD",
  "Adenosine 5'-diphosphate"              = "ADP",
  "N-Acetylglutamic acid"                 = "N-AcGlu",
  "SM(d18:1/24:1) + SM(d18:2/24:0)"       = "SM(42:2)"
)
## ---------------------------------------------------------------------------

IN_DIR  <- path.expand(IN_DIR)
OUT_DIR <- path.expand(OUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

rd <- function(q) {
  f <- file.path(IN_DIR, sprintf("FC_%s_metabolites.csv", q))
  if (!file.exists(f)) stop("Missing ", f, " - export it first.")
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}
fc <- list(Q2 = rd("Q2"), Q3 = rd("Q3"), Q4 = rd("Q4"), Q5 = rd("Q5"))
FC_C1 <- fc$Q3   # C1 == old Q3 (the lead, anchor set)

# clabel: names = CSV keys (Q*), values = NEW C display labels, ordered C1..C4
clabel <- c(Q3 = "C1 · alcohol|mut", Q4 = "C2 · alcohol|WT",
            Q2 = "C3 · mutation",    Q5 = "C4 · mut|alcohol")

sig_of <- function(d) { s <- d$p < P_CUTOFF & d$qvalue < Q_CUTOFF; s[is.na(s)] <- FALSE; s }

# C1-significant row set (== old Q3)
sigC1 <- sig_of(FC_C1) &
  !(is.na(FC_C1$INDIVIDUAL.NOTATION) |
      trimws(as.character(FC_C1$INDIVIDUAL.NOTATION)) %in% c("", "NA"))
cat("C1-significant metabolites:", sum(sigC1), "\n")   # should print 107

is_blank <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA")
lab <- as.character(FC_C1$INDIVIDUAL.NOTATION)
if ("COMMON.NAME" %in% names(FC_C1)) lab[is_blank(lab)] <- as.character(FC_C1$COMMON.NAME)[is_blank(lab)]
lab[is_blank(lab)] <- as.character(FC_C1$variable)[is_blank(lab)]
hit <- lab %in% names(RELABEL); lab[hit] <- RELABEL[lab[hit]]
lab_by_var <- setNames(lab, as.character(FC_C1$variable))

GPL <- "Glycerophospholipids"; SG <- c("Sphingolipids", "Glycerolipids")
groups <- list(
  list(file = "Fig2A_glycerophospholipids",  title = "Glycerophospholipids",            sel = FC_C1$CLASS %in% GPL),
  list(file = "Fig2B_sphingo_glycerolipids",  title = "Sphingolipids & glycerolipids",   sel = FC_C1$CLASS %in% SG),
  list(file = "Fig2C_aminoacids_other",       title = "Amino acids & other metabolites", sel = !(FC_C1$CLASS %in% c(GPL, SG)))
)

all_vars <- FC_C1$variable[sigC1]
all_fc   <- unlist(lapply(fc, function(d) d$log2_FC[match(all_vars, d$variable)]))
LIM      <- max(abs(all_fc), na.rm = TRUE)

make_panel <- function(g) {
  vars <- FC_C1$variable[sigC1 & g$sel]
  meta <- FC_C1[match(vars, FC_C1$variable), c("CLASS", "log2_FC")]
  vars <- vars[if (ORDER_BY == "class") order(meta$CLASS, meta$log2_FC) else order(meta$log2_FC)]

  long <- do.call(rbind, lapply(names(clabel), function(q) {
    d <- fc[[q]]; i <- match(vars, d$variable)
    data.frame(variable = vars,
               log2_FC  = d$log2_FC[i],
               star     = ifelse(sig_of(d)[i] %in% TRUE, "*", ""),
               contrast = clabel[[q]], stringsAsFactors = FALSE)
  }))
  long$variable <- factor(long$variable, levels = vars)
  # C1 on TOP, C4 on bottom: ggplot puts the FIRST level at the bottom,
  # so the levels run C4 -> C1 (rev of the C1..C4 display order).
  long$contrast <- factor(long$contrast, levels = rev(unname(clabel)))

  p <- ggplot(long, aes(variable, contrast, fill = log2_FC)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = star), vjust = 0.74, size = 5, fontface = "bold", color = "grey15") +
    scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 0,
                         limits = c(-LIM, LIM), na.value = "grey92", name = expression(log[2]~FC)) +
    scale_x_discrete(labels = lab_by_var) +
    labs(x = NULL, y = NULL,
         title = paste0(g$title, " - log2 fold-change across contrasts"),
         subtitle = "Metabolites significant in C1 (p<0.05, q<0.10)  ·  * = significant in that contrast") +
    theme_minimal(base_family = FONT, base_size = 13) +
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9, face = "bold", color = "black"),
          axis.text.y  = element_text(size = 14, face = "bold", color = "black"),
          panel.grid   = element_blank(),
          plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 10.5, color = "grey40"),
          legend.title  = element_text(size = 12, face = "bold"),
          legend.text   = element_text(size = 11, face = "bold"),
          legend.key.height = unit(1.1, "cm"), legend.key.width = unit(0.4, "cm"))

  w <- max(6, length(vars) * PER_W) + 1.3
  ggsave(file.path(OUT_DIR, paste0(g$file, ".png")), p, width = w, height = 4.8, dpi = 300)
  cat(sprintf("  %-32s %3d metabolites  ->  %.1f x 4.8 in\n", g$title, length(vars), w))
}

invisible(lapply(groups, make_panel))
cat("Done -> ", OUT_DIR, "\n", sep = "")
