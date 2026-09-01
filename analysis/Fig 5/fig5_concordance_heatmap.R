# NOTE (30 Aug 2026): the HEATMAP section below is superseded by
# fig5_heatmap_nc.R (Nature Communications style). Content is identical --
# same significance rule, same rows, same order. The GRAPHICAL-ABSTRACT
# section further down is NOT superseded and is still the source for that.
# =============================================================================
# Figure 5 — 4-contrast pathway heatmap, CONCORDANCE rule, in Isabel's original style.
#   Same look as heatmap.Rmd (white→red scale, blue/red ▲▼ glyphs, C-axis labels,
#   ggtext subtitle) — ONLY the significance rule and the input run are changed.
#
#   Significance per contrast:  gene FDR < 0.10  AND  metabolite nominal p < 0.05
#   Cell shade: -log10(FDR_gene) on the original white→red scale (NOT brown).
#   Data: canonical 0628 run (same as Figure 4), all four contrasts.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2)})
if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext"); library(ggtext)

IN_BASE  <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/GB_metabolomics_20260628")
DIR      <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 5/outs")
GENE_FDR <- 0.10
MET_P    <- 0.05
qmap <- c(C1 = "Q3", C2 = "Q4", C3 = "Q2", C4 = "Q5")
xlab <- c(C1 = "C1\nAlcohol\non MUT", C2 = "C2\nAlcohol\non WT",
          C3 = "C3\nMutation\n(no alc.)", C4 = "C4\nMutation\n(+alcohol)")

long <- bind_rows(lapply(names(qmap), function(cc){
  f <- file.path(IN_BASE, paste0(qmap[[cc]], "_StoufferPMAX_pooled"),
                 "modalities_pmax_gene_plus_metDec12.csv")
  read_csv(f, show_col_types = FALSE) |>
    transmute(KEGG_id, pw_label = dplyr::coalesce(KEGG_title, KEGG_id),
              contrast = cc, Z_gene, FDR_gene, Z_met, P_met, concordance)
})) |>
  mutate(pw_label = gsub(" metabolism", " metab.", pw_label, fixed = TRUE),
         pw_label = gsub(" biosynthesis", " biosynth.", pw_label, fixed = TRUE),
         sig = !is.na(FDR_gene) & !is.na(P_met) & FDR_gene < GENE_FDR & P_met < MET_P,
         contrast = factor(contrast, levels = c("C1", "C2", "C3", "C4")))

cat("Per-contrast joint hits (gene FDR<0.10 & metabolite p<0.05):\n")
print(long |> filter(sig) |> count(contrast))

ids_show <- long |> group_by(KEGG_id) |> summarise(any = any(sig), .groups = "drop") |>
  filter(any) |> pull(KEGG_id)

hd <- long |> filter(KEGG_id %in% ids_show) |>
  mutate(
    neg_log_fdr = ifelse(sig, pmin(-log10(P_met + 1e-10), 4), 0),   # shade by metabolite-arm p (the limiting arm); non-sig cells white
    gene_mark = dplyr::case_when(!sig ~ "", grepl("↑", concordance) ~ "▲",
                                 grepl("↓", concordance) ~ "▼",
                                 Z_gene > 0 ~ "▲", TRUE ~ "▼"),
    met_mark  = dplyr::case_when(!sig ~ "", grepl("↑", concordance) ~ "▲",
                                 grepl("↓", concordance) ~ "▼",
                                 Z_met  > 0 ~ "▲", TRUE ~ "▼"))

# order rows: most-significant (lowest gene FDR among sig cells) at top
ord <- hd |> filter(sig) |> group_by(KEGG_id, pw_label) |>
  summarise(m = min(FDR_gene), .groups = "drop") |> arrange(m)
hd$pw_label <- factor(hd$pw_label, levels = rev(ord$pw_label))

sub_label <- paste0(
  "<span style='color:#0D47A1; font-size:16pt; font-weight:bold;'>▲▼</span>",
  " <span style='color:#0D47A1;'>gene</span> &nbsp;&nbsp; ",
  "<span style='color:#B71C1C; font-size:16pt; font-weight:bold;'>▲▼</span>",
  " <span style='color:#B71C1C;'>metabolite</span>")

p_heat <- ggplot(hd, aes(x = contrast, y = pw_label, fill = neg_log_fdr)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_text(data = dplyr::filter(hd, gene_mark != ""), aes(label = gene_mark),
            size = 8, fontface = "bold", color = "#0D47A1", nudge_x = -0.18) +
  geom_text(data = dplyr::filter(hd, met_mark != ""), aes(label = met_mark),
            size = 8, fontface = "bold", color = "#B71C1C", nudge_x = +0.18) +
  scale_fill_gradient(low = "white", high = "#D32F2F",
                      name = expression(-log[10]~italic(p)[metabolite]),
                      na.value = "grey85", limits = c(0, 4),
                      breaks = c(0,1,2,3,4), labels = c("0","1","2","3","≥4")) +
  scale_x_discrete(labels = xlab) +
  labs(x = NULL, y = NULL, title = "Integrated pathways (concordance)", subtitle = sub_label) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = ggtext::element_markdown(size = 13, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, lineheight = 1.05,
                                   face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "right", legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.key.height = unit(1.5, "cm"), legend.key.width = unit(0.6, "cm"),
        panel.grid = element_blank(), plot.margin = margin(5, 5, 5, 5))

h <- min(48, max(6, 0.30 * length(ids_show) + 3))
ggsave(file.path(DIR, "heatmap_concordance_4contrasts.png"), p_heat, width = 11, height = h, dpi = 300, limitsize = FALSE)
ggsave(file.path(DIR, "heatmap_concordance_4contrasts.pdf"), p_heat, width = 11, height = h, device = cairo_pdf, limitsize = FALSE)
message("wrote heatmap_concordance_4contrasts.{png,pdf}")

# =============================================================================
# GRAPHICAL-ABSTRACT VERSION  (accurate, data-driven — no hand-typed arrows)
#   Reuses `long` above. Same significance rule (gene FDR<0.10 & met p<0.05).
#   Every glyph is derived from the real sign of Z_gene / Z_met, so it cannot
#   drift from the data. Curated to the story pathways, grouped into two hubs,
#   clean/large styling for placement in the graphical abstract.
# -----------------------------------------------------------------------------
# Curated story pathways: KEGG id -> short display label -> block.
# Blocks mirror the Results paragraph: C1-only (11 pw), shared (3 pw), C2-only (6 pw).
# Only CONCORDANT (both arms same direction) pathways are shown, so every arrow pair agrees.
ga_meta <- tibble::tribble(
  ~KEGG_id,    ~ga_label,                       ~hub,
  # -- significant ONLY on the mutant background (alcohol on MUT) --
  "mmu00590",  "Arachidonic acid metab.",       "MUTANT-specific",
  "mmu00591",  "Linoleic acid metab.",          "MUTANT-specific",
  "mmu00260",  "Glycine, serine & threonine",   "MUTANT-specific",
  "mmu00400",  "Phe, Tyr & Trp biosynth.",      "MUTANT-specific",
  # -- significant on BOTH backgrounds (genotype-independent alcohol response) --
  "mmu00270",  "Cysteine & methionine metab.",  "Shared",
  # -- significant ONLY on the WT background (alcohol on WT) --
  "mmu04923",  "Lipolysis regulation",          "WT-specific",
  "mmu00600",  "Sphingolipid metab.",           "WT-specific"
)

# Which contrasts to show. Default = the hypersensitivity comparison (C1 vs C2).
# Set GA_CONTRASTS <- c("C1","C2","C3","C4") to show all four.
GA_CONTRASTS <- c("C1", "C2")
ga_xlab <- c(C1 = "Alcohol\non MUTANT", C2 = "Alcohol\non WT",
             C3 = "Mutation\n(no alcohol)", C4 = "Mutation\n(+alcohol)")

ga <- long |>
  dplyr::filter(KEGG_id %in% ga_meta$KEGG_id, contrast %in% GA_CONTRASTS) |>
  dplyr::left_join(ga_meta, by = "KEGG_id") |>
  dplyr::mutate(
    shade     = ifelse(sig, pmin(-log10(P_met + 1e-10), 4), 0),  # white when not sig
    gene_mark = dplyr::case_when(!sig ~ "", Z_gene > 0 ~ "▲", TRUE ~ "▼"),
    met_mark  = dplyr::case_when(!sig ~ "", Z_met  > 0 ~ "▲", TRUE ~ "▼"),
    contrast  = factor(contrast, levels = GA_CONTRASTS),
    ga_label  = factor(ga_label, levels = rev(ga_meta$ga_label)),
    hub       = factor(hub, levels = unique(ga_meta$hub)))

cat("\nGA cells that are significant (gene FDR<0.10 & met p<0.05):\n")
print(ga |> dplyr::filter(sig) |>
        dplyr::transmute(hub, ga_label, contrast,
                         gene = gene_mark, met = met_mark,
                         P_met = signif(P_met, 2), FDR_gene = signif(FDR_gene, 2)))

ga_sub <- paste0(
  "<span style='color:#0D47A1; font-size:20pt; font-weight:bold;'>▲▼</span>",
  " <span style='color:#0D47A1;'>gene</span> &nbsp;&nbsp; ",
  "<span style='color:#B71C1C; font-size:20pt; font-weight:bold;'>▲▼</span>",
  " <span style='color:#B71C1C;'>metabolite</span>")

p_ga <- ggplot(ga, aes(x = contrast, y = ga_label, fill = shade)) +
  geom_tile(color = "black", linewidth = 0.6) +
  geom_text(data = dplyr::filter(ga, gene_mark != ""), aes(label = gene_mark),
            size = 11, fontface = "bold", color = "#0D47A1", nudge_x = -0.19) +
  geom_text(data = dplyr::filter(ga, met_mark != ""), aes(label = met_mark),
            size = 11, fontface = "bold", color = "#B71C1C", nudge_x = +0.19) +
  scale_fill_gradient(low = "white", high = "#D32F2F",
                      name = expression(-log[10]~italic(p)),
                      na.value = "grey85", limits = c(0, 4),
                      breaks = c(0, 2, 4), labels = c("0", "2", "≥4")) +
  scale_x_discrete(labels = ga_xlab, position = "top") +
  facet_grid(rows = vars(hub), scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 15) +
  theme(plot.subtitle = ggtext::element_markdown(size = 15, color = "black", hjust = 0.5),
        axis.text.x.top = element_text(size = 16, face = "bold", color = "black", lineheight = 1.05),
        axis.text.y = element_text(size = 15, color = "black"),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 12,
                                         color = "#333333"),
        strip.placement = "outside",
        panel.spacing = unit(6, "pt"),
        legend.position = "right", legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.key.height = unit(1.1, "cm"), legend.key.width = unit(0.5, "cm"),
        panel.grid = element_blank(), plot.margin = margin(8, 12, 8, 8)) +
  labs(subtitle = ga_sub)

ggsave(file.path(DIR, "GA_integrated_pathways.png"), p_ga,
       width = 8.5, height = 6.8, dpi = 400, bg = "white")
ggsave(file.path(DIR, "GA_integrated_pathways.pdf"), p_ga,
       width = 8.5, height = 6.8, device = cairo_pdf)          # vector, opens in Illustrator
if (requireNamespace("svglite", quietly = TRUE))               # editable SVG for Illustrator
  ggsave(file.path(DIR, "GA_integrated_pathways.svg"), p_ga, width = 8.5, height = 6.8)
message("wrote GA_integrated_pathways.{png,pdf,svg}")
