# =============================================================================
#  Figure 2, panels a-c — metabolites significant under alcohol on the variant
#  background, shown across all four comparisons. Nature Communications style.
#
#  CHANGED FROM fig2_metabolite_heatmaps_C1C4.R:
#    - Paths fixed. The old IN_DIR/OUT_DIR pointed at "F1"/"F2", which do not
#      exist; the fold-change CSVs live in "Fig 1".
#    - Sized for the page. Panel a was 396 mm wide (2.2x over the 180 mm limit).
#      a  -> 180 mm, 60 metabolites, 3.0 mm per column
#      b  ->  88 mm, 21 metabolites
#      c  ->  88 mm, 21 metabolites          b and c sit side by side.
#    - C1-C4 replaced with plain language on the y axis (task A2).
#    - No panel titles or subtitles; that text belongs in the legend.
#    - No baked-in panel letters (NC_TAGS is FALSE) — add them in Illustrator.
#    - Writes Source Data: exact log2 fold-change, P and q for every metabolite
#      in every comparison (task J9). A heatmap cannot carry 240 exact P values,
#      so the marks stay and the numbers go to the file.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN_DIR  <- Sys.getenv("FC_DIR",     file.path(FIGROOT, "Fig 1"))
OUT     <- Sys.getenv("FIG2_OUT",   file.path(FIGROOT, "Fig 2/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

P_CUTOFF <- 0.05; Q_CUTOFF <- 0.10
RELABEL <- c("Ursodeoxycholic + Hyodeoxycholic acid" = "UDCA + HDCA",
             "Flavin adenine dinucleotide" = "FAD",
             "Adenosine 5'-diphosphate"    = "ADP",
             "N-Acetylglutamic acid"       = "N-AcGlu",
             "SM(d18:1/24:1) + SM(d18:2/24:0)" = "SM(42:2)")

rd <- function(q) {
  f <- file.path(IN_DIR, sprintf("FC_%s_metabolites.csv", q))
  if (!file.exists(f)) stop("Missing ", f)
  utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}
fc    <- list(Q2 = rd("Q2"), Q3 = rd("Q3"), Q4 = rd("Q4"), Q5 = rd("Q5"))
ANCH  <- fc$Q3                       # alcohol on the variant background = the lead

## Plain language, top to bottom, alcohol first (tasks A1 and A2).
## Sentence case, and the exact strings in nc_style.R's CONTRAST_LABEL, so a
## reader moving between figures sees one name per comparison.
clabel <- c(Q3 = "Alcohol on variant", Q4 = "Alcohol on wild type",
            Q2 = "Variant alone",      Q5 = "Variant on alcohol")

sig_of <- function(d) { s <- d$p < P_CUTOFF & d$qvalue < Q_CUTOFF; s[is.na(s)] <- FALSE; s }
blank  <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA")

keep <- sig_of(ANCH) & !blank(ANCH$INDIVIDUAL.NOTATION)
cat(sprintf("  annotated metabolites carried into the figure: %d\n", sum(keep)))

lab <- as.character(ANCH$INDIVIDUAL.NOTATION)
if ("COMMON.NAME" %in% names(ANCH)) lab[blank(lab)] <- as.character(ANCH$COMMON.NAME)[blank(lab)]
lab[blank(lab)] <- as.character(ANCH$variable)[blank(lab)]
hit <- lab %in% names(RELABEL); lab[hit] <- RELABEL[lab[hit]]
lab_by_var <- setNames(lab, as.character(ANCH$variable))

GPL <- "Glycerophospholipids"; SG <- c("Sphingolipids", "Glycerolipids")
GROUPS <- list(
  list(file = "Fig2a_glycerophospholipids", sel = ANCH$CLASS %in% GPL,             w = 180, h = 44, key = TRUE),
  list(file = "Fig2b_sphingo_glycerolipids", sel = ANCH$CLASS %in% SG,             w = 88,  h = 44, key = FALSE),
  list(file = "Fig2c_aminoacids_other",      sel = !(ANCH$CLASS %in% c(GPL, SG)),  w = 88,  h = 44, key = FALSE))

## One colour scale across all three panels.
all_vars <- ANCH$variable[keep]
LIM <- max(abs(unlist(lapply(fc, function(d) d$log2_FC[match(all_vars, d$variable)]))), na.rm = TRUE)

long_of <- function(vars) do.call(rbind, lapply(names(clabel), function(q) {
  d <- fc[[q]]; i <- match(vars, d$variable)
  data.frame(variable = vars, log2_FC = d$log2_FC[i], p = d$p[i], q = d$qvalue[i],
             star = ifelse(sig_of(d)[i] %in% TRUE, "*", ""),
             contrast = clabel[[q]], stringsAsFactors = FALSE) }))

make_panel <- function(g) {
  vars <- ANCH$variable[keep & g$sel]
  vars <- vars[order(ANCH$log2_FC[match(vars, ANCH$variable)])]
  long <- long_of(vars)
  long$variable <- factor(long$variable, levels = vars)
  ## ggplot draws the first level at the bottom, so reverse for top-down order
  long$contrast <- factor(long$contrast, levels = rev(unname(clabel)))

  p <- ggplot(long, aes(variable, contrast, fill = log2_FC)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    geom_text(aes(label = star), vjust = 0.72, size = NC_MIN / .pt, colour = "grey15") +
    scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                         limits = c(-LIM, LIM), na.value = "grey92",
                         breaks = c(-4, 0, 4), name = expression("Log"[2]~"FC")) +
    scale_x_discrete(labels = lab_by_var, expand = c(0, 0)) +
    ## Reader-facing contrast names from nc_style.R's display layer, parsed
    ## for the I148M superscript. `clabel` stays the internal key.
    scale_y_discrete(expand = c(0, 0), labels = nc_label_expr) +
    labs(x = NULL, y = NULL) +
    theme_nc() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
          axis.text.y = element_text(size = NC_BASE),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = if (isTRUE(g$key)) "right" else "none",
          legend.justification = "centre",
          legend.title = element_text(size = 5, hjust = 0,
                                      margin = ggplot2::margin(b = 2)),
          legend.text  = element_text(size = 5),
          legend.margin = ggplot2::margin(2, 2, 2, 2)) +
    guides(fill = guide_colourbar(title.position = "top",
                                  barheight = grid::unit(24, "pt"),
                                  barwidth  = grid::unit(3.5, "pt"),
                                  frame.colour = NA, ticks.colour = "white")) +
    nc_tag("")

  save_nc(p, file.path(OUT, g$file), width_mm = g$w, height_mm = g$h)
  long$panel <- g$file
  long
}

src <- do.call(rbind, lapply(GROUPS, make_panel))

## ---- Source Data (task J9) --------------------------------------------------
src$metabolite <- lab_by_var[as.character(src$variable)]
src <- src[, c("panel", "metabolite", "variable", "contrast", "log2_FC", "p", "q")]
names(src) <- c("panel", "metabolite", "owl_code", "comparison", "log2_fold_change",
                "P_value", "q_value")
utils::write.csv(src, file.path(OUT, "SourceData_Fig2a-c.csv"), row.names = FALSE)
cat(sprintf("\n  Source Data written: %d rows (%d metabolites x 4 comparisons)\n",
            nrow(src), nrow(src) / 4))
