## ============================================================
##  Figure 2 — add-on strip: 7 SLIM boxplots, C1 contrast only
##  (C1 = Q3 = alcohol on PNPLA3-I148M background, the LEAD)
##
##  Each panel shows the two groups that define C1:
##      MUT noEtOH (ref)  vs  MUT +EtOH (test)
##  p-value shown (not stars); q calibrated against the
##  metabolome-wide p-distribution of C1 (same engine as the heatmap).
##
##  Panels (left -> right):
##    4 composite indices : Lyso-PC, Free cholesterol, Total PC, Cholesteryl-ester
##    3 individual species: LPI(22:6), LPG(22:6), LPC 20:5  [PC(0:0/20:5)]
##
##  Style: theme_bw + Arial, group colours matched to the heatmap deck.
##  Designed as a thin row to paste UNDER the existing heatmap in PPT.
##
##  HOW TO RUN: set DATA_FILE / OUT_DIR, then source().
## ============================================================

DATA_FILE <- "~/Desktop/Lipidomics project with Dr Smirnova/Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx"
OUT_DIR   <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/F2"

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(ggplot2); library(reshape2)
  library(rstatix); library(qvalue); library(patchwork)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
q_val_contr <- 0.10

# ── Load (identical to metabolite_heatmaps.R) ─────────────────────────────────
data <- read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 16)
metabol_class <- read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 0, n_max = 16)
names(data) <- names(metabol_class)
metabol_class.m <- metabol_class[, -c(2:4)] %>%
  { tmp <- t(.[-1]); colnames(tmp) <- .[["COD.OWL"]]; data.frame(tmp) }
metabol_names <- names(data)[-c(1:4)]
for (col in metabol_names)
  data[[col]][is.na(data[[col]])] <- min(data[[col]], na.rm = TRUE) / 5
data$Group <- dplyr::recode_factor(data$Sample,
  "Male 24w CD" = "CD", "AAV Luc" = "HF", "PNPLA3 MUT" = "HF + PNPLA3",
  "PNPLA3 MUT + EtOH" = "HF+PNPLA3+ALC", "PNPLA3 Luc-WD + EtOH" = "HF+ALC")

# C1 groups
REF <- "HF + PNPLA3"      # MUT noEtOH
TRT <- "HF+PNPLA3+ALC"    # MUT +EtOH

# full-distribution p-vector for C1 (for q calibration), same engine as heatmap
data.m <- reshape2::melt(data[, c("SAMPLE.ID","Group", metabol_names)],
                         id.vars = c("SAMPLE.ID","Group"))
safe_qvalue <- function(p) qvalue(p)$qvalues   # Storey q-value (qvalue package); BH fallback removed (never triggered)
C1_pvec <- data.m %>% filter(Group %in% c(REF, TRT)) %>%
  group_by(variable) %>% t_test(value ~ Group) %>% pull(p)

# ── Define the 7 panels ───────────────────────────────────────────────────────
simp <- setNames(metabol_class.m$SIMPLIFIED.NAME, rownames(metabol_class.m))
notn <- setNames(metabol_class.m$INDIVIDUAL.NOTATION, rownames(metabol_class.m))
members_simp <- function(codes) names(simp)[simp %in% codes]
members_notn <- function(x)     names(notn)[notn %in% x]

PANELS <- list(
  list(name = "Lyso-PC index",      cols = members_simp(c("LPC_1MAPC","LPC_2MAPC","LPC_1MEPC"))),
  list(name = "Free cholesterol",   cols = members_simp("Cho")),
  list(name = "Total PC index",     cols = members_simp("PC_DAPC")),
  list(name = "Cholesteryl-ester",  cols = members_simp("ChoE")),
  list(name = "LPI(22:6)",          cols = members_notn("LPI(22:6)")),
  list(name = "LPG(22:6)",          cols = members_notn("LPG(22:6)")),
  list(name = "LPC 20:5",           cols = members_notn("PC(0:0/20:5)"))
)

# Match Figure 1 exactly: COL_DOWN/-EtOH = #1565C0, COL_UP/+EtOH = #C62828, font Arial
FONT     <- "Arial"
grp_cols <- c("MUT noEtOH" = "#1565C0", "MUT +EtOH" = "#C62828")

# per-sample value (index = mean of members; single species = that column)
panel_df <- function(cols) {
  d <- data[data$Group %in% c(REF, TRT), c("Group", cols), drop = FALSE]
  val <- if (length(cols) == 1) d[[cols]] else rowMeans(d[, cols, drop = FALSE], na.rm = TRUE)
  data.frame(grp = ifelse(d$Group == REF, "MUT noEtOH", "MUT +EtOH"), val = val)
}
panel_stat <- function(df) {
  a <- df$val[df$grp == "MUT +EtOH"]; b <- df$val[df$grp == "MUT noEtOH"]
  p <- t.test(a, b, var.equal = FALSE)$p.value
  q <- tail(safe_qvalue(c(C1_pvec, p)), 1)
  list(p = p, q = q, ytop = max(df$val, na.rm = TRUE))
}

make_slim <- function(P, show_y) {
  df <- panel_df(P$cols); df$grp <- factor(df$grp, levels = c("MUT noEtOH","MUT +EtOH"))
  s  <- panel_stat(df); rng <- range(df$val, na.rm = TRUE); pad <- diff(rng) * 0.16
  lab <- if (s$p < 0.001) "p<0.001" else sprintf("p=%.3f", s$p)
  ggplot(df, aes(grp, val, fill = grp)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.8, fatten = 1.6,
                 colour = "black") +
    geom_jitter(width = 0.12, size = 2.2, shape = 21, colour = "black", stroke = 0.4) +
    annotate("segment", x = 1, xend = 2, y = s$ytop + pad*0.6, yend = s$ytop + pad*0.6,
             colour = "#C62828", linewidth = 0.9) +
    annotate("text", x = 1.5, y = s$ytop + pad*1.05, label = lab,
             size = 4.2, colour = "#C62828", fontface = "bold", family = FONT) +
    scale_fill_manual(values = grp_cols, name = NULL) +
    scale_x_discrete(labels = c("−EtOH","+EtOH")) +
    coord_cartesian(ylim = c(rng[1] - pad*0.5, s$ytop + pad*2)) +
    labs(title = P$name, x = NULL, y = if (show_y) "relative abundance" else NULL) +
    theme_classic(base_size = 14, base_family = FONT) +
    theme(plot.title  = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
          axis.text.x = element_text(size = 13, face = "bold", colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.title.y = element_text(size = 13, face = "bold"),
          axis.line  = element_line(linewidth = 0.7, colour = "black"),
          axis.ticks = element_line(linewidth = 0.7, colour = "black"),
          legend.position = "none",
          plot.margin = margin(3, 5, 3, 5))
}

plots <- lapply(seq_along(PANELS), function(i) make_slim(PANELS[[i]], show_y = (i == 1)))

strip <- wrap_plots(plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "C1 (alcohol on PNPLA3-I148M background) — lipotoxic readouts",
    caption = paste0("Two groups define C1: MUT noEtOH vs MUT +EtOH. p = Welch t. ",
                     "All panels q<0.10 (FDR, calibrated vs the metabolome-wide C1 distribution). ",
                     "Indices = mean of member species."),
    theme = theme(plot.title   = element_text(face = "bold", size = 16, family = FONT),
                  plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", family = FONT))) &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 13, face = "bold", family = FONT))

ggsave(file.path(OUT_DIR, "Figure2_C1_slim_boxplots.png"), strip, width = 14, height = 3.8, dpi = 300)
pdf(file.path(OUT_DIR, "Figure2_C1_slim_boxplots.pdf"), width = 14, height = 3.8); print(strip); invisible(dev.off())

cat("\nC1 panel stats:\n")
for (P in PANELS) { df <- panel_df(P$cols); s <- panel_stat(df)
  cat(sprintf("  %-18s p=%.4f  q=%.4f  n(ref/test)=%d/%d\n",
              P$name, s$p, s$q, sum(df$grp=="MUT noEtOH"), sum(df$grp=="MUT +EtOH"))) }
cat("\nSaved Figure2_C1_slim_boxplots.png/.pdf to", OUT_DIR, "\n")
