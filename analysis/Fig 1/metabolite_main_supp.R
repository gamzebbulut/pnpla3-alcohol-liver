## ============================================================
##  Metabolite figures — Main + Supplement
##
##  MAIN FIGURE (one page, ~9x10 in):
##    Panel A: Volcano plot — all metabolites, colored by class
##    Panel B: Class summary bar — n up / n down per class
##
##  SUPPLEMENT (tall bar charts, broken into 3 parts):
##    Supp A: Amino acids & small polar classes
##    Supp B: Glycerolipids + Sphingolipids + Sterols
##    Supp C: Glycerophospholipids
##    Each ordered by |log2FC| descending (longest bars first)
##
##  BONUS: Q2/Q3/Q4/Q5 class-level comparison bar
##    Shows n significant metabolites per class per contrast
##    Justifies focusing on Q3
##
##  HOW TO RUN: set DATA_FILE and OUT_DIR, hit Source
## ============================================================

DATA_FILE <- "~/Desktop/Lipidomics project with Dr Smirnova/Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx"
OUT_DIR   <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/F1"

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl)
  library(ggplot2); library(reshape2); library(rstatix)
  library(qvalue); library(stringr); library(forcats)
  library(ggrepel); library(patchwork)
})

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
q_val_contr <- 0.10

# ── Load & prep ───────────────────────────────────────────────────────────────
cat("Loading data...\n")
data <- read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 16)
metabol_class <- read_excel(DATA_FILE, sheet = "Data per metabolite",
                             skip = 0, n_max = 16)
names(data) <- names(metabol_class)

metabol_class.m <- metabol_class[, -c(2:4)] %>%
  { tmp <- t(.[-1]); colnames(tmp) <- .[["COD.OWL"]]; data.frame(tmp) }

metabol_names <- names(data)[-c(1:4)]
data[, metabol_names] <- data[, metabol_names] %>%
  dplyr::mutate_if(is.numeric, ~replace(., is.na(.), min(., na.rm = TRUE) / 5))

data$Group <- dplyr::recode_factor(data$Sample,
  "Male 24w CD"           = "CD",
  "AAV Luc"               = "HF",
  "PNPLA3 MUT"            = "HF + PNPLA3",
  "PNPLA3 MUT + EtOH"     = "HF+PNPLA3+ALC",
  "PNPLA3 Luc-WD + EtOH"  = "HF+ALC"
)

data.m <- reshape2::melt(
  data[, c("SAMPLE.ID", "Group", metabol_names)],
  id.vars = c("SAMPLE.ID", "Group")
)
data.m$log2_value <- log2(data.m$value)
data.m <- merge(data.m, metabol_class.m, by.x = "variable", by.y = "row.names")

safe_qvalue <- function(p) {
  qvalue(p)$qvalues   # Storey q-value (qvalue package); BH fallback removed (never triggered)
}

# ── Compute FC + t-test for one contrast ─────────────────────────────────────
compute_FC <- function(grp1, grp2, label) {
  mfc <- data.m %>%
    dplyr::filter(Group %in% c(grp1, grp2)) %>%
    dplyr::group_by(variable, Group) %>%
    dplyr::summarise(log2_value = mean(log2_value, na.rm = TRUE), .groups = "drop")
  fc <- merge(subset(mfc, Group==grp1), subset(mfc, Group==grp2),
              by="variable", suffixes=c("_ref","_trt"))
  fc$log2_FC <- fc$log2_value_trt - fc$log2_value_ref
  tt <- data.m %>%
    dplyr::filter(Group %in% c(grp1, grp2)) %>%
    dplyr::group_by(variable) %>%
    t_test(value ~ Group) %>%
    adjust_pvalue(method="BH", output.col="p_adj") %>%
    add_significance()
  tt$qvalue <- safe_qvalue(tt$p)
  tt <- merge(tt, metabol_class.m, by.x="variable", by.y="row.names")
  result <- merge(fc, tt, by="variable") %>%
    dplyr::mutate(Q = label,
                  sig = p < 0.05 & qvalue < q_val_contr,
                  neg_log10_p = -log10(p),
                  direction = dplyr::if_else(log2_FC >= 0, "Up", "Down"),
                  signif_star = dplyr::case_when(
                    p < 0.0001 ~ "****", p < 0.001 ~ "***",
                    p < 0.01   ~ "**",  p < 0.05  ~ "*", TRUE ~ ""))
  result
}

cat("Running t-tests for Q2, Q3, Q4, Q5...\n")
FC_Q2 <- compute_FC("HF",          "HF + PNPLA3",   "Q2")
FC_Q3 <- compute_FC("HF + PNPLA3", "HF+PNPLA3+ALC", "Q3")
FC_Q4 <- compute_FC("HF",          "HF+ALC",         "Q4")
FC_Q5 <- compute_FC("HF+ALC",      "HF+PNPLA3+ALC", "Q5")
# Note: Q5 has 0 significant metabolites — excluded from comparison figure

for (q in list(FC_Q2,FC_Q3,FC_Q4,FC_Q5))
  cat(sprintf("  %s: %d significant\n", q$Q[1], sum(q$sig)))

# ── Colours ───────────────────────────────────────────────────────────────────
COL_UP   <- "#C62828"
COL_DOWN <- "#1565C0"
COL_NS   <- "#BBBBBB"

# Class color palette (consistent across all figures)
class_pal <- c(
  "Amino acid"                     = "#E53935",
  "Amino acid derivatives"         = "#E53935",
  "Amino acids"                    = "#E53935",
  "Amino Acids and Derivatives"    = "#E53935",
  "Amino acids & derivatives"      = "#E53935",
  "Bile acids"                     = "#F57F17",
  "Disaccharides"                  = "#FF8F00",
  "Fatty acids"                    = "#558B2F",
  "Fatty esters"                   = "#33691E",
  "Flavin Nucleotides"             = "#00695C",
  "Glycerolipids"                  = "#42A5F5",
  "Glycerophospholipids"           = "#1565C0",
  "Purine Nucleotides"             = "#6A1B9A",
  "Sphingolipids"                  = "#9C27B0",
  "Sterols"                        = "#880E4F",
  "Sulfinic Acids and Derivatives" = "#BF360C"
)

# Collapse amino acid sub-classes for display
collapse_class <- function(x) case_when(
  x %in% c("Amino acid","Amino acid derivatives",
            "Amino acids","Amino Acids and Derivatives") ~ "Amino acids & derivatives",
  TRUE ~ x
)

theme_fig <- function(base=9) {
  theme_classic(base_size=base) +
    theme(plot.title    = element_text(face="bold", size=base+1),
          plot.subtitle = element_text(size=base-1, colour="#444444"),
          plot.caption  = element_text(size=6.5, colour="#666666",
                                       hjust=0, margin=margin(t=4)),
          panel.grid.major = element_line(colour="#EEEEEE", linewidth=0.35),
          panel.grid.minor = element_blank(),
          plot.margin = margin(6,8,4,6))
}

# ════════════════════════════════════════════════════════════════════════════════
#  MAIN FIGURE — Panel A: Volcano
# ════════════════════════════════════════════════════════════════════════════════
cat("\nBuilding main figure...\n")

# Report 107 significant, but DISPLAY only the 102 with a species-level name.
# The 5 without an INDIVIDUAL.NOTATION (MAPE09, MAPE36, ChoE_Col_02/03, FFA42)
# are counted as significant but not drawn, so every panel shows the same 102.
FC_Q3$displayable <- !is.na(FC_Q3$INDIVIDUAL.NOTATION) & FC_Q3$INDIVIDUAL.NOTATION != "NA"
# Collapse the 4 identically-colored amino-acid sub-classes into one legend entry
FC_Q3$CLASS_disp <- collapse_class(FC_Q3$CLASS)
n_sig    <- sum(FC_Q3$sig)
n_disp   <- sum(FC_Q3$sig & FC_Q3$displayable)
n_hidden <- n_sig - n_disp
cat(sprintf("  Q3 significant: %d  |  shown (annotated): %d  |  hidden (NA notation): %d\n",
            n_sig, n_disp, n_hidden))

p_thresh <- FC_Q3 %>% filter(sig) %>%
  summarise(p_max = max(p, na.rm=TRUE)) %>% pull(p_max)

label_df <- FC_Q3 %>%
  filter(sig,
         # Raise FC threshold slightly to avoid clutter near the threshold line
         # Hypotaurine (FC~1.7, p~0.05) excluded — label manually if needed
         abs(log2_FC) > 1.8 | p < 0.001,
         !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION != "NA") %>%
  mutate(display_name = str_trunc(INDIVIDUAL.NOTATION, 22, ellipsis="…"))

pA_vol <- ggplot(FC_Q3, aes(x=log2_FC, y=neg_log10_p)) +
  geom_point(data=filter(FC_Q3,!sig),
             colour=COL_NS, size=1.7, alpha=0.38, shape=16) +
  geom_hline(yintercept=-log10(p_thresh),
             linetype="dashed", colour="#888888", linewidth=0.5) +
  annotate("text", x=max(FC_Q3$log2_FC,na.rm=TRUE)*0.85,
           y=-log10(p_thresh)+0.09,
           label=paste0("q=0.10 (p=",round(p_thresh,3),")"),
           colour="#888888", size=3.8, hjust=1) +
  geom_vline(xintercept=0, linewidth=0.35, colour="#CCCCCC") +
  geom_point(data=filter(FC_Q3, sig & displayable),
             aes(colour=CLASS_disp, size=neg_log10_p),
             alpha=0.88, shape=16) +
  geom_text_repel(data=label_df,
                  aes(label=display_name, colour=CLASS_disp),
                  size=4.4, max.overlaps=20, seed=42,
                  box.padding=0.4, point.padding=0.22,
                  min.segment.length=0.15,
                  segment.size=0.3, segment.alpha=0.5,
                  show.legend=FALSE) +
  scale_colour_manual(values=class_pal, na.value=COL_NS,
                      name="Metabolite class") +
  scale_size_continuous(range=c(2.4,9), name="-log\u2081\u2080(p)",
                        breaks=c(1.3,2,3,4),
                        labels=c("0.05","0.01","0.001","0.0001")) +
  labs(title="Volcano: Q3 metabolites",
       subtitle=paste0(n_sig," significant (q<0.10) of ",nrow(FC_Q3),
                       " tested; ",n_disp," shown, ",n_hidden,
                       " unannotated not displayed  |  colored by class"),
       x="log\u2082 fold-change  (MUT+Alc vs MUT)",
       y="-log\u2081\u2080(p)") +
  theme_fig(base=15) +
  theme(legend.position="bottom",
        legend.box="vertical",
        legend.justification="left",
        legend.box.just="left",
        legend.margin=margin(t=2,b=2),
        legend.spacing.y=unit(1,"mm"),
        legend.key.size=unit(0.5,"cm"),
        legend.key.spacing.x=unit(1.2,"cm"),   # spread the 3 class columns toward full box width
        legend.text=element_text(size=13),
        legend.title=element_text(size=14,face="bold"),
        plot.subtitle=element_text(size=10),
        # size legend dropped INSIDE the panel, bottom-right corner
        legend.position.inside=c(0.99, 0.02),
        legend.justification.inside=c(1, 0),
        legend.background=element_rect(fill=alpha("white",0.85), colour=NA)) +
  guides(colour=guide_legend(title="Metabolite class",
                             title.position="top",
                             override.aes=list(size=5.5),
                             nrow=4, byrow=TRUE),
         size=guide_legend(title="-log\u2081\u2080(p)",
                           title.position="top",
                           position="inside",
                           ncol=1))

# ════════════════════════════════════════════════════════════════════════════════
#  MAIN FIGURE — Panel B: Class summary (n up / n down per class)
# ════════════════════════════════════════════════════════════════════════════════
class_sum <- FC_Q3 %>%
  filter(sig,
         !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION != "NA",
         !is.na(CLASS), CLASS != "NA") %>%
  mutate(CLASS_clean = collapse_class(CLASS)) %>%
  group_by(CLASS_clean, direction) %>%
  summarise(n = n(), .groups="drop") %>%
  mutate(n_signed = if_else(direction=="Down", -n, n))

# Order classes by total n significant
class_order_sum <- class_sum %>%
  group_by(CLASS_clean) %>%
  summarise(total=sum(n)) %>%
  arrange(total) %>%
  pull(CLASS_clean)

class_sum <- class_sum %>%
  mutate(CLASS_clean = factor(CLASS_clean, levels=class_order_sum))

pB_sum <- ggplot(class_sum, aes(x=CLASS_clean, y=n_signed, fill=direction)) +
  geom_col(width=0.70) +
  geom_hline(yintercept=0, linewidth=0.5, colour="#555555") +
  # n labels: sit just outside bar tip, black text, always readable
  geom_text(aes(label=n,
                y=n_signed + if_else(direction=="Up", 0.8, -0.8),
                hjust=if_else(direction=="Up", 0, 1)),
            size=4.2, colour="#333333", fontface="bold") +
  scale_fill_manual(values=c(Up=COL_UP, Down=COL_DOWN),
                    labels=c(Up="Increased", Down="Decreased")) +
  scale_y_continuous(labels=abs) +
  coord_flip() +
  labs(title="Significant metabolites by class",
       subtitle=paste0("Q3: ",n_sig," significant (q<0.10); ",n_disp,
                       " annotated shown  |  red = increased, blue = decreased"),
       x=NULL, y="Number of significant metabolites", fill=NULL) +
  theme_fig(base=15) +
  theme(legend.position="bottom",
        legend.key.size=unit(0.5,"cm"),
        legend.text=element_text(size=11.5),
        plot.subtitle=element_text(size=11),
        axis.text.y=element_text(size=12.5),
        panel.grid.major.x=element_line(colour="#EEEEEE",linewidth=0.35),
        panel.grid.major.y=element_blank())

# ── Standalone panel exports (drop into PowerPoint separately) ────────────────
# NOTE: at the enlarged label font, a few metabolite labels in the dense central
# cloud will overlap — ggrepel can't fully separate them at this size. The .svg
# is exported so final label positions can be nudged by hand in Illustrator;
# the .png is the auto-placed quick-look version.
ggsave(file.path(OUT_DIR,"volcano_Q3.png"),       pA_vol, width=8.4, height=8.4, dpi=300)
ggsave(file.path(OUT_DIR,"volcano_Q3.svg"),       pA_vol, width=8.4, height=8.4, dpi=300)
ggsave(file.path(OUT_DIR,"barplot_class_Q3.png"), pB_sum, width=6.2, height=5.5, dpi=300)
cat("  Saved standalone: volcano_Q3.png + barplot_class_Q3.png\n")

# ── Assemble main figure ──────────────────────────────────────────────────────
main_fig <- pA_vol + pB_sum +
  plot_layout(widths=c(1.6, 1)) +
  plot_annotation(
    title   = "Q3: Metabolite-level analysis \u2014 Alcohol on PNPLA3-I148M background",
    subtitle= "HF+PNPLA3+ALC vs HF+PNPLA3 (G3\u2212G2)  |  t-test per metabolite  |  BH + q-value FDR correction (q<0.10)",
    caption = paste0(
      "log\u2082FC = mean(log\u2082[MUT+Alc]) \u2212 mean(log\u2082[MUT]). ",
      "Dashed line = p-value threshold at q=0.10. ",
      "Point size \u221d \u2212log\u2081\u2080(p). ",
      "Labels on metabolites with |log\u2082FC|>1.5 or p<0.001."
    ),
    theme=theme(plot.title    = element_text(face="bold", size=11),
                plot.subtitle = element_text(size=8.5, colour="#444444"),
                plot.caption  = element_text(size=7, colour="#666666",
                                             hjust=0, margin=margin(t=4)))
  )

pdf(file.path(OUT_DIR,"main_fig_volcano_classsummary.pdf"), width=11, height=6)
print(main_fig); invisible(dev.off())
ggsave(file.path(OUT_DIR,"main_fig_volcano_classsummary.png"),
       main_fig, width=11, height=6, dpi=300)
cat("  Saved: main_fig_volcano_classsummary\n")

# ════════════════════════════════════════════════════════════════════════════════
#  BONUS — Q2/Q3/Q4/Q5 class-level comparison
#  Justifies focusing on Q3
# ════════════════════════════════════════════════════════════════════════════════
cat("Building Q2/Q3/Q4/Q5 class comparison...\n")

# Q5 has 0 significant metabolites — only show Q2, Q3, Q4
all_q <- bind_rows(FC_Q2, FC_Q3, FC_Q4) %>%
  filter(sig,
         !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION != "NA",
         !is.na(CLASS), CLASS != "NA") %>%
  mutate(
    CLASS_clean = collapse_class(CLASS),
    Q = factor(Q, levels=c("Q2","Q3","Q4"),
               labels=c("Q2\nMut vs WT",
                        "Q3 [LEAD]\nAlc on Mut",
                        "Q4\nAlc on WT"))
  ) %>%
  group_by(Q, CLASS_clean, direction) %>%
  summarise(n=n(), .groups="drop") %>%
  mutate(n_signed = if_else(direction=="Down", -n, n))

# Total per contrast for annotation
q_totals <- bind_rows(FC_Q2,FC_Q3,FC_Q4) %>%
  filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION != "NA") %>%
  mutate(Q=factor(Q,levels=c("Q2","Q3","Q4"),
                  labels=c("Q2\nMut vs WT","Q3 [LEAD]\nAlc on Mut",
                           "Q4\nAlc on WT"))) %>%
  count(Q, name="total") %>%
  mutate(label=paste0("n=",total))

# Class order by Q3 total
class_ord_q <- all_q %>%
  filter(grepl("Q3",as.character(Q))) %>%
  group_by(CLASS_clean) %>%
  summarise(t=sum(n)) %>%
  arrange(t) %>% pull(CLASS_clean)

all_q <- all_q %>%
  mutate(CLASS_clean=factor(CLASS_clean, levels=class_ord_q))

# dynamic per-contrast summary (annotated metabolites only) for the caption
.disp <- function(d) d[d$sig & !is.na(d$INDIVIDUAL.NOTATION) & d$INDIVIDUAL.NOTATION != "NA", ]
.fmt  <- function(d, nm) { dd <- .disp(d); up <- sum(dd$direction=="Up"); dn <- sum(dd$direction=="Down")
                           sprintf("%s n=%d (%d up, %d down)", nm, up+dn, up, dn) }
summ_str <- paste(.fmt(FC_Q2,"Q2"), .fmt(FC_Q3,"Q3"), .fmt(FC_Q4,"Q4"),
                  sprintf("Q5 n=%d", nrow(.disp(FC_Q5))), sep="  |  ")

p_compare <- ggplot(all_q, aes(x=CLASS_clean, y=n_signed, fill=direction)) +
  geom_col(width=0.72) +
  geom_hline(yintercept=0, linewidth=0.5, colour="#555555") +
  # Total n label at top of each facet
  geom_label(data=q_totals,
             aes(x=Inf, y=Inf, label=label),
             inherit.aes=FALSE,
             hjust=1.05, vjust=1.2,
             size=3, fontface="bold", fill="white",
             label.size=0.3, colour="#333333") +
  facet_wrap(~Q, nrow=1) +
  scale_fill_manual(values=c(Up=COL_UP, Down=COL_DOWN),
                    labels=c(Up="Increased", Down="Decreased")) +
  scale_y_continuous(labels=abs) +
  coord_flip() +
  labs(
    title   = "Significant metabolites per class across all contrasts (Q2\u2013Q5)",
    subtitle= paste0(
      "Q3 (alcohol on PNPLA3-I148M background) drives the broadest response across metabolite classes\n",
      "p<0.05 & q<0.10 | red=increased, blue=decreased | n shown per contrast"
    ),
    x=NULL, y="Number of significant metabolites", fill=NULL,
    caption = paste0(
      "Q2 (G2\u2212G1): mutation alone vs WT.  ",
      "Q3 (G3\u2212G2): alcohol on mutant background.  ",
      "Q4 (G4\u2212G1): alcohol on WT background.  ",
      "Q5 (G3\u2212G4): 0 significant metabolites, not shown.\n",
      "Summary (annotated only): ", summ_str
    )
  ) +
  theme_fig(base=8.5) +
  theme(
    strip.text       = element_text(face="bold", size=9),
    strip.background = element_rect(fill="#EEF2F8", colour=NA),
    legend.position  = "bottom",
    legend.key.size  = unit(0.35,"cm"),
    axis.text.y      = element_text(size=7.5),
    panel.grid.major.x = element_line(colour="#EEEEEE",linewidth=0.35),
    panel.grid.major.y = element_blank(),
    panel.spacing    = unit(6,"mm")
  ) +
  # Bold Q3 strip
  theme(strip.text=element_text(face="bold", size=9))

pdf(file.path(OUT_DIR,"comparison_Q2Q3Q4Q5_byclass.pdf"), width=11, height=6)
print(p_compare); invisible(dev.off())
ggsave(file.path(OUT_DIR,"comparison_Q2Q3Q4Q5_byclass.png"),
       p_compare, width=11, height=6, dpi=300)
cat("  Saved: comparison_Q2Q3Q4Q5_byclass\n")

# ════════════════════════════════════════════════════════════════════════════════
#  SUPPLEMENT — Horizontal bar charts in 3 parts
#  Ordered by |log2FC| descending within each class (longest bars first)
# ════════════════════════════════════════════════════════════════════════════════
cat("Building supplementary bar charts...\n")

# Order: by |log2FC| descending within class so longest bars come first
sig_ordered <- FC_Q3 %>%
  filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION != "NA") %>%
  mutate(CLASS_clean = collapse_class(CLASS),
         abs_fc = abs(log2_FC)) %>%
  arrange(CLASS_clean,
          # Down bars first (negative FC), then up bars; within each group by |FC| desc
          desc(direction == "Down"), desc(abs_fc)) %>%
  mutate(display_name = str_trunc(INDIVIDUAL.NOTATION, 30, ellipsis="…"),
         display_name = fct_inorder(display_name))

make_bar_supp <- function(df, panel_title, subtitle) {
  n <- nrow(df)
  ggplot(df, aes(x=log2_FC, y=display_name, fill=direction)) +
    geom_vline(xintercept=0, linewidth=0.5, colour="#555555") +
    geom_col(width=0.75, colour=NA) +
    geom_text(
      aes(x = log2_FC + if_else(log2_FC>=0, 0.04, -0.04),
          label = signif_star,
          hjust = if_else(log2_FC>=0, 0, 1),
          colour = direction),
      size=2.8, vjust=0.5, show.legend=FALSE
    ) +
    facet_grid(CLASS_clean~., scales="free_y", space="free_y", switch="y") +
    scale_fill_manual(values=c(Up=COL_UP,Down=COL_DOWN),
                      labels=c(Up=paste0("Increased (n=",sum(df$direction=="Up"),")"),
                               Down=paste0("Decreased (n=",sum(df$direction=="Down"),")"))) +
    scale_colour_manual(values=c(Up=COL_UP,Down=COL_DOWN), guide="none") +
    scale_x_continuous(expand=expansion(mult=c(0.05,0.12)),
                       breaks=seq(-3,3,by=0.5)) +
    labs(title=panel_title, subtitle=subtitle,
         x="log\u2082 fold-change", y=NULL, fill=NULL,
         caption=paste0("Q3 (G3\u2212G2): t-test, BH+q-value FDR (q<0.10). ",
                        "Ordered by |log\u2082FC| descending within class. ",
                        "* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001")) +
    theme_fig(base=8) +
    theme(strip.text.y.left  = element_text(size=7.5,angle=0,face="bold"),
          strip.background   = element_rect(fill="#EEF2F8",colour=NA),
          strip.placement    = "outside",
          axis.text.y        = element_text(size=7.5),
          panel.grid.major.x = element_line(colour="#EEEEEE",linewidth=0.35),
          panel.grid.major.y = element_blank(),
          panel.spacing      = unit(1.5,"mm"),
          legend.position    = "bottom",
          legend.key.size    = unit(0.38,"cm"))
}

# Part 1: Amino acids & small polars
supp_classes_A <- c("Amino acids & derivatives","Bile acids","Disaccharides",
                    "Fatty acids","Fatty esters","Flavin Nucleotides",
                    "Purine Nucleotides","Sulfinic Acids and Derivatives")
dfA <- filter(sig_ordered, CLASS_clean %in% supp_classes_A)
pSA <- make_bar_supp(dfA,
  "Supplementary: Q3 metabolite log\u2082FC — amino acids & small polar classes",
  paste0("n=",nrow(dfA)," metabolites  |  HF+PNPLA3+ALC vs HF+PNPLA3  |  ordered by |log\u2082FC| within class"))

# Part 2: Glycerolipids + Sphingolipids + Sterols
supp_classes_B <- c("Glycerolipids","Sphingolipids","Sterols")
dfB <- filter(sig_ordered, CLASS_clean %in% supp_classes_B)
pSB <- make_bar_supp(dfB,
  "Supplementary: Q3 metabolite log\u2082FC — glycerolipids, sphingolipids & sterols",
  paste0("n=",nrow(dfB)," metabolites  |  ordered by |log\u2082FC| within class"))

# Part 3: Glycerophospholipids (largest class)
dfC <- filter(sig_ordered, CLASS_clean == "Glycerophospholipids")
pSC <- make_bar_supp(dfC,
  "Supplementary: Q3 metabolite log\u2082FC — glycerophospholipids",
  paste0("n=",nrow(dfC)," metabolites  |  ordered by |log\u2082FC| descending"))

# Save supplements
save_supp <- function(p, stem, df) {
  n <- nrow(df)
  h <- n*0.20 + 2.5
  w <- 8.5
  cat(sprintf("  Saving %s: %.1f x %.1f in (%d metabolites)\n", stem, w, h, n))
  pdf(file.path(OUT_DIR, paste0(stem,".pdf")), width=w, height=h)
  print(p); invisible(dev.off())
  ggsave(file.path(OUT_DIR, paste0(stem,".png")), p, width=w, height=h, dpi=300)
  cat("  Saved:", stem, "\n")
}

save_supp(pSA, "supp_bar_aminoacids_smallpolar", dfA)
save_supp(pSB, "supp_bar_glycerolipids_sphingolipids", dfB)
save_supp(pSC, "supp_bar_glycerophospholipids", dfC)

cat("\nAll done! Files saved to:\n ", OUT_DIR, "\n\n")
cat("Summary of significant metabolites by contrast:\n")
print(as.data.frame(
  bind_rows(FC_Q2,FC_Q3,FC_Q4,FC_Q5) %>%
    group_by(Q) %>%
    summarise(n_sig=sum(sig), n_up=sum(sig&direction=="Up"),
              n_down=sum(sig&direction=="Down"))
), row.names=FALSE)

# ════════════════════════════════════════════════════════════════════════════════
#  CANDIDATE FIGURE — Q3 vs Q4 side-by-side (alcohol effect by genotype background)
#  Self-contained add-on; does NOT touch the locked Q3 standalone panels above.
#  Q3 = alcohol on mutant (HF+PNPLA3 -> +ALC);  Q4 = alcohol on WT (HF -> +ALC).
#  All four panels share identical x/y limits, point-size scale, class palette and
#  class ordering, with ONE collected legend, so the contrasts are comparable.
#  Emits: candidate_Fig1_Q3vsQ4.png / .svg   (show to Katia for feedback)
# ════════════════════════════════════════════════════════════════════════════════
cat("\nBuilding candidate Q3-vs-Q4 side-by-side figure...\n")

# ---- Shared volcano scales across BOTH contrasts ------------------------------
# x-limit driven by SIGNIFICANT points only (not the sparse grey NS cloud), so the
# margins stay tight; still shared across Q3/Q4 and never clips a significant hit.
.cmp_xm       <- max(abs(c(FC_Q3$log2_FC[FC_Q3$sig], FC_Q4$log2_FC[FC_Q4$sig])), na.rm=TRUE)
.cmp_xlim     <- c(-.cmp_xm, .cmp_xm) * 1.08
.cmp_ymax     <- max(c(FC_Q3$neg_log10_p, FC_Q4$neg_log10_p), na.rm=TRUE)
.cmp_ylim     <- c(-0.15, .cmp_ymax * 1.03)
.cmp_size_lim <- range(c(FC_Q3$neg_log10_p[FC_Q3$sig],
                         FC_Q4$neg_log10_p[FC_Q4$sig]), na.rm=TRUE)

# Union of (collapsed) significant classes -> both volcano legends match exactly
.cmp_vol_levels <- bind_rows(
    FC_Q3 %>% filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA"),
    FC_Q4 %>% filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA")
  ) %>%
  mutate(CC = collapse_class(CLASS)) %>%
  filter(!is.na(CC), CC!="NA") %>%
  distinct(CC) %>% arrange(CC) %>% pull(CC)

make_volcano_cmp <- function(FC, main_title, x_lab, extra_labels=character(0)) {
  FC <- FC %>% mutate(
    displayable = !is.na(INDIVIDUAL.NOTATION) & INDIVIDUAL.NOTATION!="NA",
    CLASS_disp  = collapse_class(CLASS))
  pth <- FC %>% filter(sig) %>% summarise(m=max(p, na.rm=TRUE)) %>% pull(m)
  lab <- FC %>% filter(sig,
                       abs(log2_FC) > 1.8 | p < 0.001 | INDIVIDUAL.NOTATION %in% extra_labels,
                       !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA") %>%
    mutate(display_name = str_trunc(INDIVIDUAL.NOTATION, 22, ellipsis="\u2026"))
  ns <- sum(FC$sig); nd <- sum(FC$sig & FC$displayable)
  ggplot(FC, aes(log2_FC, neg_log10_p)) +
    geom_point(data=filter(FC, !sig), colour=COL_NS, size=1.4, alpha=0.35, shape=16) +
    geom_hline(yintercept=-log10(pth), linetype="dashed", colour="#888888", linewidth=0.5) +
    geom_vline(xintercept=0, linewidth=0.35, colour="#CCCCCC") +
    geom_point(data=filter(FC, sig & displayable),
               aes(colour=CLASS_disp, size=neg_log10_p), alpha=0.88, shape=16) +
    geom_text_repel(data=lab, aes(label=display_name, colour=CLASS_disp),
                    size=5.0, max.overlaps=18, seed=42, box.padding=0.4,
                    point.padding=0.22, min.segment.length=0.15,
                    segment.size=0.3, segment.alpha=0.5, show.legend=FALSE) +
    scale_colour_manual(values=class_pal, na.value=COL_NS,
                        name="Metabolite class", limits=.cmp_vol_levels) +
    scale_size_continuous(range=c(2.0,8), name="-log\u2081\u2080(p)",
                          breaks=c(1.3,2,3,4), labels=c("0.05","0.01","0.001","0.0001"),
                          limits=.cmp_size_lim) +
    coord_cartesian(xlim=.cmp_xlim, ylim=.cmp_ylim) +
    labs(title=main_title,
         subtitle=paste0(ns, " significant (q<0.10); ", nd, " shown"),
         x=x_lab, y="-log\u2081\u2080(p)") +
    theme_fig(base=18) +
    theme(plot.subtitle=element_text(size=12),
          legend.text=element_text(size=13),
          legend.title=element_text(size=14, face="bold"))
}

# ---- Shared class order + count axis for the bars -----------------------------
.cmp_bar_levels <- bind_rows(FC_Q3, FC_Q4) %>%
  filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA",
         !is.na(CLASS), CLASS!="NA") %>%
  mutate(CC = collapse_class(CLASS)) %>%
  count(CC) %>% arrange(n) %>% pull(CC)

.cmp_dir_counts <- bind_rows(FC_Q3 %>% mutate(Qx="Q3"), FC_Q4 %>% mutate(Qx="Q4")) %>%
  filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA",
         !is.na(CLASS), CLASS!="NA") %>%
  mutate(CC = collapse_class(CLASS)) %>%
  count(Qx, CC, direction)
.cmp_up_max <- max(.cmp_dir_counts$n[.cmp_dir_counts$direction=="Up"])
.cmp_dn_max <- max(.cmp_dir_counts$n[.cmp_dir_counts$direction=="Down"])
# Asymmetric shared axis: trims the dead far-left (decreases never reach the up-max),
# while keeping Q3 and Q4 on ONE axis so the magnitude difference stays visible.
.cmp_bar_lo <- -(.cmp_dn_max + max(3, ceiling(.cmp_dn_max * 0.10)))
.cmp_bar_hi <-  (.cmp_up_max + max(3, ceiling(.cmp_up_max * 0.06)))

make_bar_cmp <- function(FC, contrast_tag) {
  cs <- FC %>%
    filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA",
           !is.na(CLASS), CLASS!="NA") %>%
    mutate(CC = collapse_class(CLASS)) %>%
    group_by(CC, direction) %>% summarise(n=n(), .groups="drop") %>%
    tidyr::complete(CC=.cmp_bar_levels, direction=c("Up","Down"), fill=list(n=0)) %>%
    mutate(n_signed = if_else(direction=="Down", -n, n),
           CC = factor(CC, levels=.cmp_bar_levels))
  ns <- sum(FC$sig)
  ggplot(cs, aes(CC, n_signed, fill=direction)) +
    geom_col(width=0.70) +
    geom_hline(yintercept=0, linewidth=0.5, colour="#555555") +
    geom_text(data=filter(cs, n>0),
              aes(label=n, y=n_signed + if_else(direction=="Up", 0.8, -0.8),
                  hjust=if_else(direction=="Up", 0, 1)),
              size=4.2, colour="#333333", fontface="bold") +
    scale_fill_manual(values=c(Up=COL_UP, Down=COL_DOWN),
                      labels=c(Up="Increased", Down="Decreased")) +
    scale_x_discrete(labels=function(x) str_wrap(x, width=14)) +
    scale_y_continuous(labels=abs, limits=c(.cmp_bar_lo, .cmp_bar_hi)) +
    coord_flip() +
    labs(title=paste0(contrast_tag, " \u2014 significant by class"),
         subtitle=paste0(ns, " significant (q<0.10)"),
         x=NULL, y="n significant metabolites", fill=NULL) +
    theme_fig(base=16) +
    theme(axis.text.y=element_text(size=12),
          plot.subtitle=element_text(size=11),
          panel.grid.major.x=element_line(colour="#EEEEEE", linewidth=0.35),
          panel.grid.major.y=element_blank())
}

# ---- Build panels + assemble (columns = contrast, rows = volcano/bar) ----------
vol_Q3 <- make_volcano_cmp(FC_Q3, "Q3 \u2014 alcohol on mutant (PNPLA3-I148M)",
                           "log\u2082 fold-change  (MUT+Alc vs MUT)")
vol_Q4 <- make_volcano_cmp(FC_Q4, "Q4 \u2014 alcohol on WT",
                           "log\u2082 fold-change  (WT+Alc vs WT)")
bar_Q3 <- make_bar_cmp(FC_Q3, "Q3")
bar_Q4 <- make_bar_cmp(FC_Q4, "Q4")

# ---- Save each panel legend-less, as PNG + SVG (assemble by hand in PPTX) ------
save_both <- function(p, stem, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stem,".png")), p, width=w, height=h, dpi=300)
  ggsave(file.path(OUT_DIR, paste0(stem,".svg")), p, width=w, height=h, dpi=300)
}
no_legend <- theme(legend.position="none")

# Panel export sizes match a stacked layout on a ~7.5in content width:
#   volcanoes = 2/3 width (stacked), bars = 1/3 width (stacked). Placing each near
#   1:1 (don't shrink) keeps the text as large as it renders here.
save_both(vol_Q3 + no_legend, "candidate_volcano_Q3", 4.8, 4.6)
save_both(vol_Q4 + no_legend, "candidate_volcano_Q4", 4.8, 4.6)
save_both(bar_Q3 + no_legend, "candidate_bar_Q3",     2.7, 4.6)
save_both(bar_Q4 + no_legend, "candidate_bar_Q4",     2.7, 4.6)

# ---- Shared legends as SEPARATE pieces (place each by hand) -------------------
# Invisible points generate the scales; override.aes makes the keys visible.
.vol_leg_df <- bind_rows(
    FC_Q3 %>% filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA"),
    FC_Q4 %>% filter(sig, !is.na(INDIVIDUAL.NOTATION), INDIVIDUAL.NOTATION!="NA")
  ) %>% mutate(CLASS_disp = collapse_class(CLASS))

# (a) CLASS legend — wide horizontal strip, spread across the figure width
vol_class_legend <- ggplot(.vol_leg_df,
                     aes(log2_FC, neg_log10_p, colour=CLASS_disp)) +
  geom_point(alpha=0) +
  scale_colour_manual(values=class_pal, na.value=COL_NS,
                      name="Metabolite class", limits=.cmp_vol_levels) +
  guides(colour=guide_legend(title="Metabolite class", title.position="top",
                             override.aes=list(size=6.5, alpha=1), nrow=2, byrow=TRUE)) +
  theme_void(base_size=14) +
  theme(legend.position="bottom", legend.justification="center",
        legend.key.spacing.x=unit(1.0,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=15, face="bold"),
        legend.margin=margin(4,4,4,4))
save_both(vol_class_legend, "candidate_volcano_class_legend", 14, 2.0)

# (b) SIZE legend (-log10 p) — standalone VERTICAL block to drop into a corner
size_legend <- ggplot(.vol_leg_df, aes(log2_FC, neg_log10_p, size=neg_log10_p)) +
  geom_point(alpha=0) +
  scale_size_continuous(range=c(2,8), name="-log\u2081\u2080(p)",
                        breaks=c(1.3,2,3,4), labels=c("0.05","0.01","0.001","0.0001"),
                        limits=.cmp_size_lim) +
  guides(size=guide_legend(title="-log\u2081\u2080(p)", title.position="top",
                           override.aes=list(alpha=1), ncol=1)) +
  theme_void(base_size=14) +
  theme(legend.position="right", legend.justification="center",
        legend.text=element_text(size=14),
        legend.title=element_text(size=15, face="bold"),
        legend.margin=margin(4,4,4,4))
save_both(size_legend, "candidate_size_legend", 2.3, 3.3)

# ---- Increased / Decreased legend — standalone VERTICAL block -----------------
bar_legend <- ggplot(data.frame(direction=factor(c("Up","Down"), levels=c("Up","Down")), n=1),
                     aes(direction, n, fill=direction)) +
  geom_col(alpha=0) +
  scale_fill_manual(values=c(Up=COL_UP, Down=COL_DOWN),
                    labels=c(Up="Increased", Down="Decreased"), name=NULL) +
  guides(fill=guide_legend(override.aes=list(alpha=1), ncol=1)) +
  theme_void(base_size=14) +
  theme(legend.position="right", legend.text=element_text(size=15),
        legend.margin=margin(4,4,4,4))
save_both(bar_legend, "candidate_bar_legend", 2.5, 1.7)

cat("  Saved candidate panels: volcano_Q3/Q4, bar_Q3/Q4 (PNG+SVG)\n")
cat("  Saved legends: volcano_class (wide), size (vertical), bar (vertical) (PNG+SVG)\n")
