# =============================================================================
#  Figure 5 panel b — the amino-acid response, grouped by pathway
#
#  Replaces the six-boxplot strip. Six boxplots showed six metabolites chosen to
#  span panel a's pathways; a reader could not tell whether those six were the
#  story or six picks from a longer list. This panel shows EVERY amino-acid and
#  sulfur metabolite measured in the contrast, grouped by the pathway it belongs
#  to, with significance marked. Nothing is selected on the result.
#
#  Grouping is by textbook pathway membership, fixed before looking at any
#  direction or magnitude, and it deliberately crosses the platform's own class
#  labels: hypotaurine is filed under sulfinic acids, S-adenosylhomocysteine
#  under nucleoside analogues and glutathione under peptides, yet all three sit
#  on the transsulfuration axis and belong beside cysteine. Grouping by the
#  vendor's classes would have hidden the strongest signal in the figure.
#
#  WHAT THE PANEL IS FOR. Almost every amino acid rises. Against that
#  background the informative observation is not that things go up -- it is the
#  ONE that goes down. Cysteine falls while the metabolite immediately upstream
#  (cystathionine) and both products downstream (hypotaurine, taurine) rise.
#  That is a flux signature, and it does not rest on any single P value.
#
#  P and q are read from the Figure 1 table so the two figures cannot disagree.
#
#  ORIENTATION. Horizontal: metabolites along x, effect on y, pathway groups as
#  facet strips across the top. The vertical version was 116 x 132 mm and could
#  not sit under the 180 x 148 mm cross-contrast heatmap inside the 225 mm page.
#  Names are abbreviated only where the platform's string is long enough to eat
#  the panel; the full strings are in the Source Data.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
FC      <- Sys.getenv("FIG5_FC",    file.path(FIGROOT, "Fig 1/FC_Q3_metabolites.csv"))
OUT     <- Sys.getenv("FIG5_OUT",   file.path(FIGROOT, "Fig 5/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

P_CUT <- 0.05; Q_CUT <- 0.10

## ---- pathway groups, by biochemistry, fixed in advance ----------------------
## Order sets the row order in the figure: the transsulfuration axis first,
## because that is the group that carries a direction, then the rest.
GROUPS <- list(
  "Transsulfuration and taurine" = c(
    "Methionine", "S-Adenosylhomocysteine", "Homocysteine", "Cystathionine",
    "Cysteine", "Cystine-double charge", "Hypotaurine", "Taurine",
    "Glutathione (reduced)", "Glutathione (oxidized)"),
  "One-carbon donors" = c("Serine", "Glycine", "Sarcosine", "Threonine", "L-Homoserine"),
  "Urea cycle and nitrogen" = c(
    "N-Acetylglutamic acid", "L-Citrulline", "Aspartic Acid", "Asparagine",
    "Glutamic Acid", "Glutamine", "N-Acetyl-glutamine", "Pyroglutamic acid"),
  "BCAA" = c("Valine", "Leucine", "Isoleucine"),   # short label: the facet is only 3 wide
  "Aromatic" = c("Phenylalanine", "Tyrosine", "Tryptophan"),
  "Other amino acids" = c(
    "Alanine", "Proline", "Lysine", "Histidine", "1-Methyl-L-histidine",
    "Gamma-Aminobutyric acid", "DL-2-Aminoadipic acid", "Hippuric acid"))

z <- utils::read.csv(FC, check.names = FALSE, stringsAsFactors = FALSE)
z$nm <- trimws(as.character(z$INDIVIDUAL.NOTATION))

d <- do.call(rbind, lapply(names(GROUPS), function(g) {
  m <- GROUPS[[g]]
  found <- m[m %in% z$nm]
  if (length(setdiff(m, found)))
    message("  not on the platform, skipped: ", paste(setdiff(m, found), collapse = ", "))
  r <- z[match(found, z$nm), ]
  data.frame(group = g, name = found, logFC = r$log2_FC, p = r$p, q = r$qvalue,
             stringsAsFactors = FALSE)
}))
d$sig <- !is.na(d$q) & !is.na(d$p) & d$p < P_CUT & d$q < Q_CUT
d$group <- factor(d$group, levels = names(GROUPS))
## Within each group, order by effect so the direction reads down the column.
d <- d[order(d$group, d$logFC), ]
d$name <- factor(d$name, levels = unique(d$name))

## Shorter display names. Rotated labels are the height budget in a horizontal
## panel, so anything over ~14 characters is abbreviated. Full names in the
## Source Data file.
SHORT <- c("S-Adenosylhomocysteine" = "SAH",
           "Glutathione (reduced)"  = "GSH",
           "Glutathione (oxidized)" = "GSSG",
           "Cystine-double charge"  = "Cystine",
           "N-Acetylglutamic acid"  = "N-Ac-glutamate",
           "N-Acetyl-glutamine"     = "N-Ac-glutamine",
           "L-Citrulline"           = "Citrulline",
           "L-Homoserine"           = "Homoserine",
           "Aspartic Acid"          = "Aspartate",
           "Glutamic Acid"          = "Glutamate",
           "Pyroglutamic acid"      = "Pyroglutamate",
           "Hippuric acid"          = "Hippurate",
           "Gamma-Aminobutyric acid" = "GABA",
           "DL-2-Aminoadipic acid"  = "2-Aminoadipate",
           "1-Methyl-L-histidine"   = "1-Me-histidine")
for (k in names(SHORT)) levels(d$name)[levels(d$name) == k] <- SHORT[[k]]

## Within each group the strongest effect goes LEFT, so each block reads
## downhill and the one negative bar in the transsulfuration block sits at its
## right-hand end where the eye stops.
d <- d[order(d$group, -d$logFC), ]
d$name <- factor(as.character(d$name), levels = unique(as.character(d$name)))

## Group labels wrapped, or they overrun their facet strip.
levels(d$group) <- vapply(levels(d$group), function(x)
  paste(strwrap(x, width = 16), collapse = "\n"), character(1))

cat(sprintf("\n  %d metabolites in %d pathway groups; %d significant (P < %.2f and q < %.2f)\n",
            nrow(d), nlevels(d$group), sum(d$sig), P_CUT, Q_CUT))
cat(sprintf("  direction, all measured: %d up, %d down (sign test P = %s)\n",
            sum(d$logFC > 0), sum(d$logFC < 0),
            fmt_p(stats::binom.test(sum(d$logFC > 0), nrow(d), 0.5)$p.value)))
cat("  significant and DECREASED: ",
    paste(as.character(d$name)[d$sig & d$logFC < 0], collapse = ", "), "\n")

p <- ggplot(d, aes(name, logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NC_INK) +
  geom_segment(aes(x = name, xend = name, y = 0, yend = logFC, colour = logFC > 0),
               linewidth = 0.35) +
  geom_point(aes(fill = logFC > 0, size = sig, alpha = sig),
             shape = 21, colour = NC_INK, stroke = 0.2) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c("TRUE" = NC_RED_TINT, "FALSE" = NC_BLUE_TINT), guide = "none") +
  scale_fill_manual(values = c("TRUE" = NC_RED, "FALSE" = NC_BLUE), guide = "none") +
  scale_size_manual(values = c("TRUE" = 2.1, "FALSE" = 1.2),
                    labels = c("TRUE" = "Yes", "FALSE" = "No"), name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
  ## The caption states the BACKGROUND, not just the exposure. "no alcohol vs
  ## alcohol" alone would let a reader think the groups are pooled across
  ## genotypes; both groups here carry PNPLA3 I148M. Named in plain language
  ## only -- the C1..C4 shorthand is not to appear anywhere a reader can see
  ## (revision task A2, Sanyal).
  ## The wording is the approved one from nc_style.R's display layer, so the
  ## caption and every axis in the paper say the same thing; parsed, so the
  ## superscript matches the figures.
  labs(y = expression("Log"[2]~"fold-change"), x = NULL,
       caption = expression("Alcohol on the "*PNPLA3^"I148M"*" background")) +
  guides(size = guide_legend(override.aes = list(fill = NC_GREY, alpha = c(0.45, 1)),
                             title = "Significant")) +
  theme_nc(grid_y = TRUE) +
  theme(strip.background = element_blank(),
        ## 8 pt, not 6. The co-authors' complaint about this panel was that
        ## nobody could read the metabolite names, and 6 pt was the smallest
        ## type in the paper. The labels are ROTATED, so their size is limited
        ## by panel HEIGHT, not by the 4.8 mm between columns -- at 8 pt the 37
        ## names still clear each other; the panel just grows 72 -> 84 mm.
        strip.text.x = element_text(size = NC_BASE, colour = NC_INK2, lineheight = 0.9,
                                    margin = margin(b = 1)),
        panel.spacing.x = grid::unit(3, "pt"),
        axis.text.x = element_text(size = NC_BASE, angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top", legend.key.size = grid::unit(7, "pt"),
        legend.title = element_text(size = NC_MIN, colour = NC_INK2),
        plot.caption = element_text(size = NC_MIN, colour = NC_INK2, hjust = 0)) +
  nc_tag("b")

## 84 mm, up from 72, to carry the 8 pt names. Figure 5 then totals
## 135 (panel a) + 84 = 219 mm against the 225 mm page maximum, leaving 6 mm
## for the gutter and the panel letters. If that is too tight, the lever is
## panel a's row height in fig5_heatmap_nc.R (H <- 26 + 5.2 * n): 5.0 mm per
## row gives back 4 mm and 21 rows will still take it.
save_nc(p, file.path(OUT, "Fig5b_aminoacid_groups"), width_mm = NC_W2, height_mm = 84)

utils::write.csv(d[, c("group","name","logFC","p","q","sig")],
                 file.path(OUT, "SourceData_Fig5b_groups.csv"), row.names = FALSE)

cat("\n-- the transsulfuration axis, in pathway order --\n")
AX <- c("Methionine","S-Adenosylhomocysteine","Homocysteine","Cystathionine",
        "Cysteine","Hypotaurine","Taurine","Glutathione (reduced)")
for (m in AX) { r <- z[z$nm == m, ]
  if (nrow(r)) cat(sprintf("  %-24s %+6.2f   P = %-9s q = %-7s %s\n", m, r$log2_FC[1],
      fmt_p(r$p[1]), fmt_p(r$qvalue[1]),
      if (!is.na(r$qvalue[1]) && r$p[1] < P_CUT && r$qvalue[1] < Q_CUT) "*" else "")) }
