# =============================================================================
#  Figure 4 panel c  -- gene-arm-only top KEGG pathways (slim, right of the row)
#  Supplement       -- metabolite-arm-only top KEGG pathways (full width)
#
#  WHY A GENE-ONLY PANEL AT ALL, when panel b already shows both arms?
#  Because the joint integration is metabolite-limited by construction: a
#  pathway can only be a joint hit if it carries measured metabolites, and the
#  metabolite arm spans only +/-2.8 against the gene arm's +/-6.7. The joint
#  list is therefore dominated by what the metabolome could see -- the amino
#  acid cluster. The strongest transcriptional signals in the experiment,
#  PPAR signalling foremost, never enter it. This panel is the gene arm on its
#  own terms, and it is where PPARa (Z = -6.67, the most decreased pathway in
#  the whole gene arm) becomes discussable.
#
#  Read left to right, the row then tells one story: pipeline -> what the two
#  arms agree on -> what the transcriptome says by itself. Arachidonic acid and
#  linoleic acid metabolism sit in panel b, concordantly down; PPAR signalling
#  sits at the top of panel c, most decreased. All three are the same Cyp
#  biology, which is what the Cyp panel then resolves gene by gene.
#
#  Colour carries CATEGORY, not direction: which side of zero a bar sits on
#  already says increased or decreased, so spending colour on that would repeat
#  the axis. Categories come from pathway_categories.csv (editable, with a note
#  column recording every judgement call).
#
#  NOTE the ranking is by |Z| with NO filtering. Vascular smooth muscle
#  contraction, autoimmune thyroid disease and allograft rejection do reach the
#  top 12; they are left in, in grey, rather than quietly dropped. Their KEGG
#  gene sets are largely MHC and cytoskeletal genes that are not liver-specific,
#  and grey is doing exactly the job it should here -- de-emphasising them
#  without hiding them.
#
#  NOTE the two arms are NOT on a shared x axis and are no longer shown side by
#  side. Both are Z scores, but they are separate measurements of separate
#  modalities, not two halves of one comparison.
# =============================================================================

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("FIG4_IN",    file.path(ROOT, "GB_metabolomics_20260628/Q3_StoufferPMAX_pooled"))
OUT     <- Sys.getenv("FIG4_OUT",   file.path(FIGROOT, "Fig 4/nc_panels"))
CAT     <- Sys.getenv("FIG4_CATS",  file.path(FIGROOT, "Fig 4/pathway_categories.csv"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

TOP_MAIN <- 12   # slim panel c -- 12 is what fits at 5 pt in 54 mm
TOP_SUPP <- 15   # supplement has the width for 15

combo  <- utils::read.csv(file.path(IN, "modalities_pmax_gene_plus_metDec12.csv"),
                          stringsAsFactors = FALSE)
titles <- stats::setNames(combo$KEGG_title, combo$KEGG_id)
cats   <- utils::read.csv(CAT, stringsAsFactors = FALSE)
catof  <- stats::setNames(cats$category, cats$KEGG_id)

gene <- utils::read.csv(file.path(IN, "genes_liver_ct_stoufferDec12.csv"), stringsAsFactors = FALSE)
gene <- data.frame(id = gene$KEGG_id, title = gene$KEGG_title, Z = gene$Z_gene_liver,
                   stringsAsFactors = FALSE)

metab <- utils::read.csv(file.path(IN, "metabolite_limma_ZDec12.csv"), stringsAsFactors = FALSE)
metab$id <- sub(":.*$", "", metab$pathway)
metab <- data.frame(id = metab$id, title = unname(titles[metab$id]), Z = metab$Z_met,
                    stringsAsFactors = FALSE)
metab$title[is.na(metab$title) | metab$title == ""] <- metab$id[is.na(metab$title) | metab$title == ""]

top_n <- function(d, n) {
  d <- d[!is.na(d$Z) & !is.na(d$title) & d$title != "", ]
  d <- d[order(-abs(d$Z)), ][seq_len(min(n, nrow(d))), ]
  d[order(d$Z), ]                      # most decreased at the bottom
}

## KEGG titles run to 35 characters; unwrapped they leave a slim panel no room
## and clip the axis title. Wrap rather than truncate -- a shortened pathway
## name is a different pathway name.
wrap_lab <- function(x, width) vapply(x, function(s)
  paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)

arm_panel <- function(d, xlab, tag, n, wrap, lab_size, show_legend,
                      legend_nrow = 2) {
  d <- top_n(d, n)
  d$title <- factor(wrap_lab(d$title, wrap), levels = wrap_lab(d$title, wrap))
  d$cat <- unname(catof[d$id]); d$cat[is.na(d$cat)] <- "Other"
  d$cat <- factor(d$cat, levels = names(NC_CAT))
  lim <- max(abs(d$Z)) * 1.12
  ggplot(d, aes(Z, title, fill = cat)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = NC_INK) +
    geom_col(width = 0.7, colour = NC_INK, linewidth = 0.1) +
    scale_fill_manual(values = NC_CAT, name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(-lim, lim)) +
    labs(x = xlab, y = NULL) +
    theme_nc(grid_y = FALSE) +
    guides(fill = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
    ## A legend at "top" is laid out INSIDE the fixed device height, so the bars
    ## compress to make room and the panel stays 86 mm -- which is what we want
    ## here. It also spans the full panel width rather than the plot area, so
    ## the long category names fit without being shortened.
    theme(legend.position = if (show_legend) "top" else "none",
          legend.justification = "left",
          legend.key.size = grid::unit(7, "pt"),
          legend.margin = ggplot2::margin(0, 0, 3, 0),
          axis.text.y = element_text(size = lab_size, lineheight = 0.82)) +
    nc_tag(tag)
}

## Panel c DOES NOT carry the category key (reverted 31 Aug). The key went in
## here on 30 Aug because the standalone strip kept not getting placed, and the
## assembled figure had five bar colours and nothing to read them by. The cost
## was visible in the assembly: a five-row vertical key ate about 16 mm off an
## 86 mm panel and the bars were noticeably squeezed.
##
## The key is a strip again -- `Fig4_category_key.svg`, written by
## fig4_concordance_scatter_nc.R at 132 x 12 mm, which is panel b + gutter +
## panel c. That is the honest home for it: panels b and c encode the SAME five
## categories in the same colours, so one key serves both, and a strip spanning
## exactly those two panels says which panels it governs without a caption.
## Placing it is now a numbered step in LAYOUT_Fig4_NatComms.md rather than
## something to remember.
##
## Short axis title: the full "Z (Stouffer, snRNA-seq)" overruns 54 mm at the
## 8 pt base. The method is named in the figure legend.
g <- arm_panel(gene, "Z  (Stouffer)", "c", TOP_MAIN, 18, NC_MIN, FALSE)
save_nc(g, file.path(OUT, "Fig4c_gene_arm_slim"), width_mm = 54, height_mm = 86)

## Supplement: the metabolite arm on its own, full single-column width. It is
## no longer in the main figure -- panel b carries the metabolite arm now. This
## one KEEPS its own key: it is placed alone in the supplement, with no shared
## strip beside it to read the colours by.
m <- arm_panel(metab, "Z  (ssGSEA, metabolome)", "a", TOP_SUPP, 26, NC_BASE, TRUE)
save_nc(m, file.path(OUT, "SuppFig_metabolite_arm"), width_mm = 88, height_mm = 86)

src <- rbind(cbind(arm = "gene",       top_n(gene,  TOP_MAIN)),
             cbind(arm = "metabolite", top_n(metab, TOP_SUPP)))
src$category <- unname(catof[src$id]); src$category[is.na(src$category)] <- "Other"
utils::write.csv(src, file.path(OUT, "SourceData_Fig4c_and_supp.csv"), row.names = FALSE)

cat("\n-- top pathways plotted --\n")
for (a in c("gene", "metabolite")) {
  d <- src[src$arm == a, ]; d <- d[order(-abs(d$Z)), ]
  cat(sprintf("  %s arm:\n", a))
  for (i in seq_len(nrow(d))) cat(sprintf("     %+6.2f  %-38s %s\n", d$Z[i], d$title[i], d$category[i]))
}
