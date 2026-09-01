# NOTE (30 Aug 2026): superseded by fig5_aminoacid_boxplots_nc.R.
# That version also FIXES A BUG here: the per-panel q appended the metabolite's
# own p to the metabolome-wide p vector before re-running qvalue, but that p was
# already in the vector, so every q was estimated with one value duplicated.
## ============================================================
##  Figure 5 — add-on strip: 6 SLIM boxplots, C1 amino-acid metabolites
##  (C1 = Q3 = alcohol on PNPLA3-I148M background)
##  Style matched to fig2_C1_slim_boxplots.R (theme_classic, Arial,
##  blue/red, p-value with red bar, q calibrated vs metabolome-wide C1).
##  Reads the SAME metabolite matrix the integration used (met_mat_contrast.csv),
##  so Figure 5B is consistent with Figure 5A.
##  HOW TO RUN: set DATA_DIR / OUT_DIR, then source().
## ============================================================

DATA_DIR <- "~/Desktop/Lipidomics project with Dr Smirnova/GB_metabolomics_20260628/Q3_StoufferPMAX_pooled"
OUT_DIR  <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 5/outs"

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
safe_q <- function(p) tryCatch(qvalue::qvalue(p)$qvalues, error = function(e) p.adjust(p, "BH"))

# ── Load metabolite matrix (cpd_id x sample) + sample groups ──────────────────
mat <- readr::read_csv(file.path(DATA_DIR, "met_mat_contrast.csv"), show_col_types = FALSE)
grp <- readr::read_csv(file.path(DATA_DIR, "sample_groups_met.csv"), show_col_types = FALSE)

REF <- "MUT_WSD_NAL"   # MUT noEtOH
TRT <- "MUT_WSD_PAL"   # MUT +EtOH
g <- setNames(grp$group, grp$sample)

long <- mat |>
  tidyr::pivot_longer(-cpd_id, names_to = "sample", values_to = "val") |>
  mutate(group = g[sample]) |>
  filter(group %in% c(REF, TRT))

# full-distribution p-vector for C1 (q calibration, all metabolites)
C1_pvec <- long |> group_by(cpd_id) |>
  summarise(p = tryCatch(t.test(val[group==TRT], val[group==REF], var.equal=FALSE)$p.value,
                         error=function(e) NA_real_), .groups="drop") |>
  pull(p)
C1_pvec <- C1_pvec[!is.na(C1_pvec)]

# ── 6 panels: KEGG compound id -> display name ───────────────────────────────
PANELS <- tibble::tribble(
  ~cpd,       ~name,
  "C00037",   "Glycine",
  "C00188",   "Threonine",
  "C00213",   "Sarcosine",
  "C00079",   "Phenylalanine",
  "C00183",   "Valine",
  "C00097",   "Cysteine")

FONT <- "Arial"
grp_cols <- c("MUT noEtOH" = "#1565C0", "MUT +EtOH" = "#C62828")

make_slim <- function(cpd, nm, show_y) {
  df <- long |> filter(cpd_id == cpd) |>
    transmute(grp = factor(ifelse(group==REF,"MUT noEtOH","MUT +EtOH"),
                           levels=c("MUT noEtOH","MUT +EtOH")), val)
  a <- df$val[df$grp=="MUT +EtOH"]; b <- df$val[df$grp=="MUT noEtOH"]
  p <- t.test(a, b, var.equal = FALSE)$p.value
  q <- tail(safe_q(c(C1_pvec, p)), 1)
  rng <- range(df$val, na.rm = TRUE); pad <- diff(rng) * 0.16; ytop <- max(df$val)
  lab <- if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
  ggplot(df, aes(grp, val, fill = grp)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.8, fatten = 1.6, colour = "black") +
    geom_jitter(width = 0.12, size = 2.2, shape = 21, colour = "black", stroke = 0.4) +
    annotate("segment", x = 1, xend = 2, y = ytop + pad*0.6, yend = ytop + pad*0.6,
             colour = "#C62828", linewidth = 0.9) +
    annotate("text", x = 1.5, y = ytop + pad*1.05, label = lab,
             size = 4.2, colour = "#C62828", fontface = "bold", family = FONT) +
    scale_fill_manual(values = grp_cols, name = NULL) +
    scale_x_discrete(labels = c("−EtOH","+EtOH")) +
    coord_cartesian(ylim = c(rng[1] - pad*0.5, ytop + pad*2)) +
    labs(title = nm, x = NULL, y = if (show_y) "relative abundance" else NULL) +
    theme_classic(base_size = 14, base_family = FONT) +
    theme(plot.title  = element_text(face="bold", size=14, hjust=0.5, colour="black"),
          axis.text.x = element_text(size=13, face="bold", colour="black"),
          axis.text.y = element_text(size=11, colour="black"),
          axis.title.y = element_text(size=13, face="bold"),
          axis.line  = element_line(linewidth=0.7, colour="black"),
          axis.ticks = element_line(linewidth=0.7, colour="black"),
          legend.position = "none", plot.margin = margin(3,5,3,5))
}

plots <- lapply(seq_len(nrow(PANELS)), function(i)
  make_slim(PANELS$cpd[i], PANELS$name[i], show_y = (i == 1)))

strip <- wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "C1 (alcohol on PNPLA3-I148M background) — amino-acid metabolites",
    caption = paste0("Two groups define C1: MUT noEtOH vs MUT +EtOH. p = Welch t; ",
                     "q calibrated vs the metabolome-wide C1 distribution."),
    theme = theme(plot.title = element_text(face="bold", size=16, family=FONT),
                  plot.caption = element_text(size=10, hjust=0.5, face="italic", family=FONT))) &
  theme(legend.position = "bottom",
        legend.text = element_text(size=13, face="bold", family=FONT))

ggsave(file.path(OUT_DIR, "Fig5_C1_aminoacid_slim_boxplots.png"), strip, width = 14, height = 3.8, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig5_C1_aminoacid_slim_boxplots.pdf"), strip, width = 14, height = 3.8, device = cairo_pdf)

cat("\nC1 amino-acid panel stats:\n")
for (i in seq_len(nrow(PANELS))) {
  df <- long |> filter(cpd_id == PANELS$cpd[i])
  a <- df$val[df$group==TRT]; b <- df$val[df$group==REF]
  p <- t.test(a,b,var.equal=FALSE)$p.value; q <- tail(safe_q(c(C1_pvec,p)),1)
  cat(sprintf("  %-14s p=%.4f  q=%.4f\n", PANELS$name[i], p, q))
}
