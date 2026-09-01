# =============================================================================
#  Figure 4 panel e — the Cyp strip
#
#  ONE CLAIM, THREE ROWS:
#    row 1  how far each Cyp moved with alcohol, in the population panel d maps
#    rows 2-4  which of the three hit pathways each Cyp belongs to
#
#  The point of the panel is that arachidonic acid metabolism, linoleic acid
#  metabolism and PPAR signaling look like three findings and are largely ONE:
#  their KEGG gene sets overlap in the cytochrome P450s, and it is the Cyps that
#  are down. Filled circles under a blue tile say that in one look.
#
#  ORIENTATION. Genes run along x because panel e is a wide, short strip under a
#  full-width pathway map. VERTICAL <- TRUE flips it to genes down the y axis if
#  the layout changes.
#
#  WHICH NUMBERS -- READ THIS, THEY ARE NOT THE SAME THING
#    SOURCE "periportal"  edgeR glmQLFTest on the periportal population, the
#                         route Steve's README specifies, one column per mouse.
#                         This is the DEFAULT because it is the same test that
#                         colours the map in panel d and stars the boxplots. One
#                         population, one model, one set of numbers, so nothing
#                         in Figure 4 can disagree with anything else in it.
#    SOURCE "animal"      moderated limma on per-mouse means across the six
#                         hepatocyte populations. Also correct, also mouse-level
#                         -- but it is a different population set from panel d,
#                         so choosing it means panel e and panel d are computed
#                         on different data and the legend has to say so.
#  Both are "the approved way" in the sense that matters -- the replication unit
#  is the mouse in both. They differ in which cells were measured, not in how.
#
#  SHOWING BOTH ARMS (TWO_ARMS <- TRUE)
#  A second row for alcohol on the wild-type background is easy to draw and
#  hard to caption. Two rows side by side invite the reader to compare them and
#  conclude the variant is what makes alcohol suppress these Cyps. It does not:
#  the difference contrast reaches nothing, and for the Cyp4a genes the wild-type
#  arm is the LARGER effect (Cyp4a14 -3.03 against -4.44). So the two-row version
#  tells a true story only if the Results says that out loud -- alcohol suppresses
#  this programme on both backgrounds, somewhat less so with the variant present.
#  If the plan is still to keep the wild-type arm in the supplement, leave this
#  FALSE. If it goes in the main figure, it needs a Results sentence and Katia's
#  sign-off, because it re-opens the cross-background comparison.
#
#  Significance is an OUTLINE on the tile, not a saturation step: at n = 5 vs 4
#  in one population few genes clear a transcriptome-wide FDR, and tinting the
#  rest would read as missing data. Same grammar as the population grid.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
SOURCE     <- "periportal"      # or "animal" -- see above
POP        <- "Periportal-hepatocyte"
POP_FULL   <- "Periportal hepatocyte"    # as the RDS names it
## TWO_ARMS draws alcohol on the wild-type background as a SECOND effect row.
## Before switching it on, read the note under "SHOWING BOTH ARMS" below --
## this is not a neutral presentation choice.
TWO_ARMS   <- FALSE
ARM_LAB    <- c(variant = "Alcohol on variant", wild_type = "Alcohol on wild type")
QUESTION   <- "Q3"
## 0.10, not 0.05 (changed 31 Aug), and it is the manuscript-wide threshold:
## the gene arm enters the integration at FDR < 0.10, LSEA is called at 0.10,
## and Figure 1 reports q < 0.10. At 0.05 exactly one of the 33 drawn Cyps was
## outlined, which reads as a null panel; at 0.10 it is two, Cyp4a10 and
## Cyp4a14 -- both PPAR-alpha omega-hydroxylases, which is the same programme
## panel d shows suppressed. The panel's claim is the coordinated shift
## (23 of 33 decreased), not per-gene significance, either way.
FDR_CUT    <- 0.10
LIMIT      <- NULL              # NULL = symmetric to the data; or set e.g. 3
VERTICAL   <- FALSE
## 170 x 45, set 30 Aug to fit the space under the pathway map. Type stays at
## NC_BIG: the panel was originally drawn at NC_MIN (6 pt), below theme_nc's own
## 7 pt default, and read as an afterthought beside the rest of the figure. The
## rotated gene labels are what set the floor on height -- if they start to
## collide with the tiles, take height back before you take type size back.
W_MM       <- 170
H_MM       <- 45
SHOW_FDR   <- FALSE             # print the exact FDR under each significant gene
PATHWAYS   <- c(AA = "Arachidonic acid metabolism",
                LA = "Linoleic acid metabolism",
                PPAR = "PPAR signaling")
## "alcohol effect" said nothing about WHICH background, which is the same
## omission task A2 removes with C1-C4. Name it.
EFFECT_ROW <- unname(ARM_LAB["variant"])
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
OUT     <- Sys.getenv("FIG4_OUT",  file.path(FIGROOT, "Fig 4/nc_panels"))
CYP     <- Sys.getenv("FIG4_CYP",  file.path(FIGROOT, "Fig 4/shared_Cyps_AA_LA_PPAR_C1.csv"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

cyp <- utils::read.csv(CYP, stringsAsFactors = FALSE)
for (k in c("AA", "LA", "PPAR")) cyp[[k]] <- as.logical(cyp[[k]])

SPEC <- NULL
if (TWO_ARMS) {
  cand <- c(file.path(dirname(OUTS), "fig4_outs_ALLCELLS", "cyp_specificity_qlf_long.csv"),
            file.path(OUTS, "cyp_specificity_qlf_long.csv"))
  f2 <- cand[file.exists(cand)][1]
  if (is.na(f2)) stop("TWO_ARMS needs cyp_specificity_qlf_long.csv -- run fig4_cyp_specificity_qlf.R")
  SPEC <- utils::read.csv(f2, stringsAsFactors = FALSE)
  SPEC <- SPEC[SPEC$cell_type == POP_FULL & SPEC$detected %in% c(TRUE, "TRUE") &
                 SPEC$contrast %in% c("variant", "wild_type"), ]
  cat(sprintf("  two arms, %s, from %s\n", POP_FULL, basename(f2)))
}

if (identical(SOURCE, "periportal")) {
  f <- file.path(OUTS, sprintf("singlepop_%s_%s_genes.csv", POP, QUESTION))
  if (!file.exists(f)) stop("Run fig4_single_population.R first -- missing ", f)
  de <- utils::read.csv(f, stringsAsFactors = FALSE)
  de <- data.frame(SYMBOL = de$SYMBOL, logFC = de$logFC, FDR = de$FDR)
  SRC_LAB <- "periportal hepatocytes, edgeR quasi-likelihood F test"
} else {
  f <- file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION))
  if (!file.exists(f)) stop("Run fig4_animal_level_limma.R first -- missing ", f)
  de <- utils::read.csv(f, stringsAsFactors = FALSE)
  de <- data.frame(SYMBOL = de$SYMBOL, logFC = de$logFC_animal, FDR = de$FDR_animal)
  SRC_LAB <- "six hepatocyte populations, per-mouse means, moderated limma"
}
cat(sprintf("  effect source: %s\n", SRC_LAB))

## Genes with no measured effect carry nothing to draw. They are counted in the
## legend rather than drawn as empty columns.
i <- match(cyp$symbol, de$SYMBOL)
cyp$logFC_use <- de$logFC[i]; cyp$FDR_use <- de$FDR[i]
n_all <- nrow(cyp)
d <- cyp[!is.na(cyp$logFC_use), ]
d <- d[order(d$logFC_use), ]
LEV <- d$symbol
cat(sprintf("  %d of %d Cyps detected and drawn; %d at FDR < %.2f\n",
            nrow(d), n_all, sum(d$FDR_use < FDR_CUT, na.rm = TRUE), FDR_CUT))

LIM <- if (!is.null(LIMIT)) LIMIT else
  max(abs(c(d$logFC_use, if (TWO_ARMS) SPEC$logFC)), na.rm = TRUE)

## Long form. Row 1 is the effect tile; rows 2-4 are membership circles. They
## share a y scale so the strip is ONE ggplot -- a montage would drift out of
## column alignment on any re-render, which is the same reason panel e's
## membership matrix was drawn inside its own plot rather than beside it.
if (TWO_ARMS) {
  ## Rows stay ordered by the variant arm so the strip still reads left to right
  ## as the reported contrast does; the wild-type row is then read against it.
  EFF_ROWS <- unname(ARM_LAB[c("variant", "wild_type")])
  tile <- do.call(rbind, lapply(c("variant", "wild_type"), function(a) {
    x <- SPEC[SPEC$contrast == a, ]
    i <- match(d$symbol, x$symbol)
    data.frame(gene = factor(d$symbol, levels = LEV), row = unname(ARM_LAB[a]),
               logFC = x$logFC[i], sig = !is.na(x$FDR[i]) & x$FDR[i] < FDR_CUT,
               stringsAsFactors = FALSE) }))
  tile <- tile[!is.na(tile$logFC), ]
} else {
  EFF_ROWS <- EFFECT_ROW
  tile <- data.frame(gene = factor(d$symbol, levels = LEV), row = EFFECT_ROW,
                     logFC = d$logFC_use, sig = !is.na(d$FDR_use) & d$FDR_use < FDR_CUT,
                     stringsAsFactors = FALSE)
}
ROWS <- c(EFF_ROWS, unname(PATHWAYS))
mem <- do.call(rbind, lapply(names(PATHWAYS), function(k) data.frame(
  gene = factor(d$symbol, levels = LEV), row = unname(PATHWAYS[k]),
  member = d[[k]], stringsAsFactors = FALSE)))
tile$row <- factor(tile$row, levels = rev(ROWS))
mem$row  <- factor(mem$row,  levels = rev(ROWS))

p <- ggplot(mapping = aes(gene, row)) +
  geom_tile(data = tile, aes(fill = logFC), colour = "white", linewidth = 0.4) +
  geom_tile(data = tile[tile$sig, ], fill = NA, colour = NC_INK, linewidth = 0.55) +
  ## Membership is drawn as two point layers rather than a mapped fill, because
  ## a second fill scale cannot coexist with the effect scale in one ggplot.
  geom_point(data = mem[mem$member, ],  shape = 16, size = 1.8, colour = NC_INK) +
  geom_point(data = mem[!mem$member, ], shape = 1,  size = 1.8, colour = "#C4C4C4",
             stroke = 0.35) +
  scale_fill_gradient2(low = NC_BLUE, mid = "white", high = NC_RED, midpoint = 0,
                       limits = c(-LIM, LIM), na.value = "#F2F2F2",
                       name = expression("Log"[2]~"fold-change"),
                       breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(expand = expansion(add = 0.04)) +
  ## a little y padding, or the circles on the top row are sliced in half by
  ## the panel edge
  scale_y_discrete(expand = expansion(add = 0.14)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = grid::unit(34, "pt"),
                                barheight = grid::unit(5.5, "pt"),
                                title.position = "left", ticks.colour = NC_INK)) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = NC_BIG, angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic"),
        axis.text.y = element_text(size = NC_BIG),
        legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = NC_BIG),
        legend.text = element_text(size = NC_BIG),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "pt")) +
  nc_tag("e")

if (SHOW_FDR) {
  lb <- tile[tile$sig, ]
  p <- p + geom_text(data = lb, aes(gene, row, label = fmt_p(d$FDR_use[match(lb$gene, d$symbol)])),
                     size = NC_MIN / .pt * 0.8, colour = NC_INK, vjust = -1.6,
                     inherit.aes = FALSE)
}
if (VERTICAL) p <- p + coord_flip() +
  theme(axis.text.x = element_text(size = NC_MIN, angle = 0, hjust = 0.5, face = "plain"),
        axis.text.y = element_text(size = NC_MIN, face = "italic"))

## One more row needs one more row's height; the gene labels take the rest.
H_USE <- H_MM + if (TWO_ARMS) 8 else 0
save_nc(p, file.path(OUT, sprintf("Fig4e_cyp_strip_%s%s", SOURCE,
                                  if (TWO_ARMS) "_both_arms" else "")),
        if (VERTICAL) H_USE else W_MM, if (VERTICAL) W_MM else H_USE)

## ---- what the panel says ---------------------------------------------------
cat(sprintf("\n  %d of %d drawn Cyps decreased\n", sum(d$logFC_use < 0), nrow(d)))
for (k in names(PATHWAYS))
  cat(sprintf("    %-28s %2d genes, mean log2FC %+.2f\n", PATHWAYS[k],
              sum(d[[k]]), mean(d$logFC_use[d[[k]]])))
sh <- d$symbol[rowSums(d[, c("AA","LA","PPAR")]) >= 2]
cat(sprintf("    in two or more of the three : %d  (%s)\n", length(sh),
            paste(sh, collapse = ", ")))
if (TWO_ARMS) {
  v <- tile$logFC[tile$row == ARM_LAB["variant"]]
  w <- tile$logFC[tile$row == ARM_LAB["wild_type"]]
  cat(sprintf("\n  both arms in %s: mean %+.2f on variant, %+.2f on wild type\n",
              POP_FULL, mean(v, na.rm = TRUE), mean(w, na.rm = TRUE)))
  cat("  DO NOT caption this as specificity. The difference contrast reaches\n")
  cat("  nothing and the Cyp4a genes move further on the wild-type background.\n")
}
cat("\n  LEGEND MATERIAL: the three pathways share their cytochrome P450 members,\n")
cat("  so they are largely one finding seen three times. Say shared membership,\n")
cat("  not three independent pathways.\n")
