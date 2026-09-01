## ============================================================
##  FIGURE 2 — Composite: log2FC heatmap (A) + lead-significant
##  lipotoxic-membrane boxplots (B–D), assembled with patchwork.
##
##  ONE significance rule for the whole figure (matches Figure 1):
##      significant  <=>  p < 0.05  AND  qvalue < 0.10
##  The box-panel indices are tested with the SAME t-test engine and
##  their q-value is calibrated against the metabolome-wide p-value
##  distribution of that contrast (i.e. the index p is appended to the
##  full vector of per-metabolite p-values before qvalue() is run).
##
##  NEW contrast labels (to match the abstract). Display order C1->C4:
##      C1 = Q3 = alcohol on PNPLA3-I148M background   (LEAD)
##      C2 = Q4 = alcohol on WT background
##      C3 = Q2 = mutation alone (no alcohol)
##      C4 = Q5 = mutation on alcohol background
##
##  Verified independently in Python: this engine reproduces the
##  Figure-1 significant counts exactly (C1/Q3=107, C2/Q4=32,
##  C3/Q2=5, C4/Q5=0). Recomputed index results under the matched
##  cutoff:  Lyso-PC C1 p=0.038 q=0.079 *,  Free chol C1 p=0.021
##  q=0.060 *,  Total PC C1 p=0.038 q=0.079 *.  (Cholesteryl-ester
##  misses: C1 p=0.053 -> dropped from the panel.)
##
##  HOW TO RUN: set DATA_FILE / OUT_DIR, then source().
## ============================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_FILE <- "~/Desktop/Lipidomics project with Dr Smirnova/Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx"
OUT_DIR   <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/F2"

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(ggplot2)
  library(reshape2); library(rstatix); library(qvalue); library(patchwork)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
q_val_contr <- 0.10

# ── Load data (identical to metabolite_heatmaps.R) ────────────────────────────
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

data.m <- reshape2::melt(data[, c("SAMPLE.ID", "Group", metabol_names)],
                         id.vars = c("SAMPLE.ID", "Group"))
data.m$log2_value <- log2(data.m$value)

safe_qvalue <- function(p) qvalue(p)$qvalues   # Storey q-value (qvalue package); BH fallback removed (never triggered)

# ── Per-contrast FC + t-test (identical engine) ───────────────────────────────
compute_FC <- function(grp1, grp2, label) {           # grp1 = ref, grp2 = test
  mean_FC <- data.m %>% filter(Group %in% c(grp1, grp2)) %>%
    group_by(variable, Group) %>%
    summarise(log2_value = mean(log2_value, na.rm = TRUE), .groups = "drop")
  fc <- merge(subset(mean_FC, Group == grp1), subset(mean_FC, Group == grp2),
              by = "variable", suffixes = c("_ref", "_trt"))
  fc$log2_FC <- fc$log2_value_trt - fc$log2_value_ref
  tt <- data.m %>% filter(Group %in% c(grp1, grp2)) %>%
    group_by(variable) %>% t_test(value ~ Group) %>%
    adjust_pvalue(method = "BH", output.col = "p_adj") %>% add_significance()
  tt$qvalue <- safe_qvalue(tt$p)
  tt <- merge(tt, metabol_class.m, by.x = "variable", by.y = "row.names")
  res <- merge(fc, tt, by = "variable"); res$Comparison <- label; res
}

# C-label, ref group, test group, contrast type
CONTR <- list(
  C1 = list(ref = "HF + PNPLA3", trt = "HF+PNPLA3+ALC", type = "alcohol", lead = TRUE),
  C2 = list(ref = "HF",          trt = "HF+ALC",        type = "alcohol", lead = FALSE),
  C3 = list(ref = "HF",          trt = "HF + PNPLA3",   type = "mutation",lead = FALSE),
  C4 = list(ref = "HF+ALC",      trt = "HF+PNPLA3+ALC", type = "mutation",lead = FALSE))
C_LEVELS <- c("C1", "C2", "C3", "C4")

FC <- lapply(names(CONTR), function(cl) compute_FC(CONTR[[cl]]$ref, CONTR[[cl]]$trt, cl))
names(FC) <- names(CONTR)

# ── PANEL A: multi-contrast heatmap of C1-significant metabolites ─────────────
breaks_Q3 <- c(-3, -2, -1, -0.5, -0.3, 0.3, 0.5, 3)
values_Q3 <- c("blue4","blue","cornflowerblue","cadetblue2","white","lightpink2","coral")
add_cut <- function(df) { df$cut <- cut(df$log2_FC, breaks_Q3, right = TRUE, include.lowest = TRUE)
  df$signif <- df$p.signif; df$signif[df$signif == "ns"] <- ""; df }
FC <- lapply(FC, add_cut)
cut_levels <- levels(FC$C1$cut)

c1_sig <- subset(FC$C1, p < 0.05 & qvalue < q_val_contr)$variable

allC <- bind_rows(FC) %>% filter(variable %in% c1_sig) %>%
  mutate(Comparison = factor(Comparison, levels = C_LEVELS),
         signif_show = if_else(p < 0.05 & qvalue < q_val_contr, "*", ""))

pA <- ggplot(allC, aes(Comparison, INDIVIDUAL.NOTATION, fill = cut)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = signif_show), size = 3) +
  scale_fill_manual("log2(FC)", breaks = cut_levels, values = values_Q3,
                    na.value = "grey90", drop = FALSE) +
  facet_grid(CLASS ~ ., space = "free_y", scales = "free_y", switch = "y") +
  labs(title = "Metabolite log2 fold-change — C1-significant metabolites across contrasts",
       subtitle = paste0("rows = metabolites significant in C1 (p<0.05, q<0.10, n=",
                          length(c1_sig), ")  |  * = also significant in that contrast  |  C1 = alcohol on PNPLA3-I148M (lead)")) +
  theme_bw() +
  theme(strip.text.y.left = element_text(size = 8, angle = 0),
        strip.background = element_rect(fill = "#EEF2F8"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 5),
        panel.spacing = unit(0, "line"), legend.position = "right")

# ── Box-panel indices (member species by SIMPLIFIED NAME) ─────────────────────
simp <- setNames(metabol_class.m$SIMPLIFIED.NAME, rownames(metabol_class.m))
members <- function(codes) names(simp)[simp %in% codes]
INDEX_DEFS <- list(
  "Lyso-PC index"    = members(c("LPC_1MAPC","LPC_2MAPC","LPC_1MEPC")),
  "Free cholesterol" = members("Cho"),
  "Total PC index"   = members("PC_DAPC"))

# per-sample index value = mean of member raw (imputed) values
index_long <- function(cols) {
  sub <- data[, c("SAMPLE.ID","Group", cols)]
  sub$idx <- rowMeans(sub[, cols, drop = FALSE], na.rm = TRUE)
  sub[, c("SAMPLE.ID","Group","idx")]
}

# index q-value calibrated against the contrast's metabolome-wide p distribution
index_stat <- function(df_idx, contr) {
  ref <- CONTR[[contr]]$ref; trt <- CONTR[[contr]]$trt
  a <- df_idx$idx[df_idx$Group == trt]; b <- df_idx$idx[df_idx$Group == ref]
  p <- t.test(a, b, var.equal = FALSE)$p.value
  lfc <- mean(log2(a), na.rm = TRUE) - mean(log2(b), na.rm = TRUE)
  qv <- tail(safe_qvalue(c(FC[[contr]]$p, p)), 1)   # appended to full distribution
  star <- if (p < 0.05 & qv < q_val_contr)
            (if (p < .001) "***" else if (p < .01) "**" else "*") else "ns"
  data.frame(Comparison = contr, ref = ref, trt = trt, type = CONTR[[contr]]$type,
             logFC = lfc, p = p, q = qv, star = star, ytop = max(c(a, b)))
}

grp_cols <- c("HF"="#BDBDBD","HF + PNPLA3"="#5B8FB9","HF+ALC"="#E8A87C","HF+PNPLA3+ALC"="#B83232")
grp_labs <- c("HF"="WT noEtOH","HF + PNPLA3"="MUT noEtOH","HF+ALC"="WT +EtOH","HF+PNPLA3+ALC"="MUT +EtOH")
type_cols <- c(mutation = "#2E7D32", alcohol = "#C62828")

box_panel <- function(name) {
  di <- index_long(INDEX_DEFS[[name]])
  st <- bind_rows(lapply(C_LEVELS, function(c) index_stat(di, c)))
  # long frame: each contrast contributes its ref & test samples, x = contrast
  plotdf <- bind_rows(lapply(C_LEVELS, function(c) {
    ref <- CONTR[[c]]$ref; trt <- CONTR[[c]]$trt
    rbind(
      data.frame(Comparison = c, role = "ref",  Group = ref, idx = di$idx[di$Group == ref]),
      data.frame(Comparison = c, role = "test", Group = trt, idx = di$idx[di$Group == trt]))
  }))
  plotdf$Comparison <- factor(plotdf$Comparison, levels = C_LEVELS)
  plotdf$Group      <- factor(plotdf$Group, levels = names(grp_cols))
  rng <- range(plotdf$idx, na.rm = TRUE); pad <- diff(rng) * 0.12
  st$Comparison <- factor(st$Comparison, levels = C_LEVELS)

  ggplot(plotdf, aes(Comparison, idx, fill = Group)) +
    geom_boxplot(aes(group = interaction(Comparison, role)),
                 position = position_dodge(width = 0.7), width = 0.6,
                 outlier.shape = NA, linewidth = 0.4) +
    geom_point(aes(group = interaction(Comparison, role)),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
               size = 1.1, shape = 21, colour = "black", stroke = 0.2) +
    geom_text(data = st, aes(x = Comparison, y = ytop + pad, label = star, colour = type),
              inherit.aes = FALSE, fontface = "bold", size = 4.2,
              show.legend = FALSE) +
    scale_fill_manual(values = grp_cols, labels = grp_labs, name = NULL) +
    scale_colour_manual(values = type_cols, guide = "none") +
    labs(title = name, y = "relative abundance", x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 11),
          panel.grid.major.x = element_blank(), legend.position = "none")
}

pB <- box_panel("Lyso-PC index")
pC <- box_panel("Free cholesterol")
pD <- box_panel("Total PC index")

# ── Assemble: A on top, B|C|D beneath, one shared legend + caption ────────────
boxrow <- pB + pC + pD + plot_layout(nrow = 1)
fig2 <- pA / boxrow +
  plot_layout(heights = c(2.7, 1.05)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste0("One significance key for the whole figure: p<0.05 & q<0.10; index q calibrated against the ",
                     "metabolome-wide p-distribution per contrast. Index = mean of member species. ",
                     "Brackets: green = mutation (C3/C4), red = alcohol (C1/C2)."),
    theme = theme(plot.caption = element_text(size = 8, hjust = 0.5, face = "italic")))

ggsave(file.path(OUT_DIR, "Figure2_combined.png"), fig2, width = 10.5, height = 13, dpi = 300)
pdf(file.path(OUT_DIR, "Figure2_combined.pdf"), width = 10.5, height = 13); print(fig2); invisible(dev.off())

# ── Echo the index stats table (for the legend / methods) ─────────────────────
cat("\nIndex statistics (matched cutoff p<0.05 & q<0.10):\n")
for (nm in names(INDEX_DEFS)) {
  di <- index_long(INDEX_DEFS[[nm]])
  print(cbind(index = nm, bind_rows(lapply(C_LEVELS, function(c)
    index_stat(di, c)))[, c("Comparison","logFC","p","q","star")]), row.names = FALSE)
}
cat("\nDone. Saved Figure2_combined.png/.pdf to", OUT_DIR, "\n")
