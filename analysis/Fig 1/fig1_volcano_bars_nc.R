# =============================================================================
#  Figure 1, panels a-c — Nature Communications style
#    a  volcano, alcohol on the variant background
#    b  volcano, alcohol on the wild-type background   (axes identical to a)
#    c  significant metabolites per class, BOTH comparisons dodged, full width
#
#  Replaces volcano_standalone.R / volcano_Q4_standalone.R / barplot_*.R, and
#  supersedes the four-panel version kept in
#  "archived versions/fig1_volcano_bars_nc_2panel_20260830.R".
#
#  WHY THE BARS WERE MERGED (30 Aug 2026)
#    The old panels b and d shared a class order and a count axis so the two
#    comparisons could be read against each other. That was right, and it is
#    kept here — but it forced the wild-type panel to carry a row for every
#    class the variant had. Five of twelve rows were empty, and the axis ran
#    -60 to +60 for a largest bar of 17: the panel was about 75% white, in the
#    opener. Dodging the two comparisons inside ONE full-width panel keeps the
#    shared order and the shared axis, drops the empty rows, and puts 107
#    against 32 in a single reading instead of asking the eye to hold one panel
#    while it looks at another.
#
#  CHANGED, from the four-panel version:
#    - Panel letters: the wild-type volcano is now b, the bars are c.
#      Old nc_panels/Fig1b_classcounts_* and Fig1d_classcounts_* are STALE.
#    - Classes carrying fewer than MIN_CLASS_TOTAL significant metabolites
#      across both comparisons pool into one row. Deliberately NOT called
#      "Other": the volcano key already has an "Other" with a different
#      membership (it pools sterols and bile acids too), and two different
#      "Other"s in one figure is a trap. Membership is printed on every run.
#    - x tick labels are abs(): the old panels printed negative counts.
#    - A class present in one comparison and absent in the other prints an
#      explicit 0 rather than nothing. Alcohol on the wild type has no
#      significant sphingolipid, and that is a result, not a gap.
#
#  CARRIED OVER, unchanged:
#    - Q3/Q4 and C1-C4 gone; comparisons named in plain language (task A2).
#    - No titles, no baked-in panel letters (Illustrator does those).
#    - Point size is FIXED; it used to double-encode -log10(P).
#    - Glycerophospholipids stay muted grey in the volcanoes (Isabel's
#      treatment) so the smaller classes are visible. Panel c counts them full.
#    - The dashed line is drawn at the P that actually corresponds to q = 0.10
#      in each comparison, computed from the data, not a hard-coded 0.05.
#    - Volcano class colours avoid red and blue, which mean increased and
#      decreased. Worst adjacent pair dE 12.1 (protanopia), 29.1 (normal).
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING — the only block you should need to touch
## ---------------------------------------------------------------------------
MERGED_BARS     <- TRUE          # FALSE restores the old separate panels b / d
COLOUR_BY       <- c("contrast")   # decided 30 Aug -- see the note below
MIN_CLASS_TOTAL <- 2             # pool classes below this (both comparisons)
MINOR_LABEL     <- "Minor classes"
SHOW_ZERO       <- TRUE          # print 0 where a class is absent in one arm
PAIR_ORDER      <- c("Alcohol on wild type", "Alcohol on variant")  # top->bottom
## Class rows are ordered by the VARIANT arm even though the wild-type arm is
## drawn first. Ordering by the arm with 27 annotated species would sort the
## panel on the sparser of the two; the richer arm gives a stable order that
## both rows are then read against.
ORDER_BY        <- "Alcohol on variant"
## Panel c is half width, sharing its row with the class heatmap (panel d).
## The three things that broke at 88 mm the first time are fixed rather than
## avoided: the key moved from inside the panel to the top, the longest class
## name is shortened by CLASS_AXIS below, and the left pad grew so the count on
## the widest decreased bar cannot clip.
BAR_W_MM        <- 88
## Matched to fig1d_lsea_classes_nc.R's H_MM so the two bottom half-width panels
## are the same height in the assembled figure. Change both or neither.
BAR_H_MM        <- 50
## Key inside the panel, anchored to its bottom-right corner (npc units).
## Set LEGEND_POS <- "top" to put it back above the panel.
LEGEND_POS      <- c(0.995, 0.11)
LEGEND_JUST     <- c(1, 0)
## Axis-only shortening. Panel d spells the classes out, so nothing is lost.
CLASS_AXIS      <- c("Amino acids & derivatives" = "Amino acids")
DODGE_W         <- 0.72

##  COLOUR_BY, the one open design decision in this panel
##    "contrast"  — bar colour is the comparison (variant dark, wild type
##                  light, both in the blue family because both are alcohol
##                  contrasts, per CONTRAST_COL in nc_style.R). Direction is
##                  carried by the side of the zero line alone. This is the
##                  same reasoning as Fig 4b/4c: the axis already says
##                  direction, so spending colour on it repeats the axis.
##    "direction" — bar colour is red increased / blue decreased, as the old
##                  panels b and d had it, with the comparison named in small
##                  text at the left of each bar. Keeps Figure 1 as the place
##                  where the red/blue convention is established.
##  DECIDED 30 Aug: "contrast". The rule the figures now follow is that where a
##  diverging axis already carries direction, colour is free to carry category;
##  where there is no axis to carry it -- Fig 4e's tiles, which have no sign --
##  colour carries direction. Fig 1c has the axis, so colour is the comparison.
##  Put "direction" back in COLOUR_BY to re-render the alternative for comparison.
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN_DIR  <- Sys.getenv("FC_DIR",     file.path(FIGROOT, "Fig 1"))
OUT     <- Sys.getenv("FIG1_OUT",   file.path(FIGROOT, "Fig 1/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

P_CUT <- 0.05; Q_CUT <- 0.10
N_LABEL <- 12          # metabolites labeled per volcano, strongest by P
                       # 12, not 14: at the larger type below, 14 repelled
                       # labels start pushing each other off the panel.

## ---- volcano size and type ------------------------------------------------
## The volcanoes were 88 x 88 with 8 pt type and read as small text in a large
## empty square: most of the cloud sits in the middle third and the corners are
## bare. Shorter panel, larger type, and axis padding cut to almost nothing.
VOLC_W   <- 88
VOLC_H   <- 76
VOLC_TXT <- 9    # axis titles, tick labels and the class key
VOLC_LAB <- 8    # the repelled metabolite names

## Reading order is alcohol on the wild type FIRST, then on the variant, which
## is CONTRAST_ORDER in nc_style.R and revision task A1: alcohol, then variant,
## then the combination. The volcano panels follow it left to right and panel c
## follows it top to bottom within each class.
CONTRASTS <- list(
  list(key = "Q4", panel = "a", name = "Alcohol on wild type"),
  list(key = "Q3", panel = "b", name = "Alcohol on variant"))

## ---- classes ---------------------------------------------------------------
collapse_class <- function(x) {
  x <- as.character(x)
  x[x %in% c("Amino acid", "Amino acid derivatives", "Amino acids",
             "Amino Acids and Derivatives")] <- "Amino acids & derivatives"
  x
}
MUTED <- "Glycerophospholipids"
KEYCL <- c("Amino acids & derivatives", "Glycerolipids", "Sphingolipids")
CLASS_COL <- c("Amino acids & derivatives" = "#6A3FA0",
               "Glycerolipids"             = "#E65100",
               "Sphingolipids"             = "#00897B",
               "Other class"               = "#4A4A4A",
               "Glycerophospholipids"      = "#9E9E9E")
CLASS_ORDER <- c("Amino acids & derivatives", "Glycerolipids", "Sphingolipids",
                 "Other class", "Glycerophospholipids")
## shorter strings for the volcano legend only; panel c spells them out
CLASS_SHORT <- c("Amino acids & derivatives" = "Amino acids",
                 "Glycerolipids" = "Glycerolipids", "Sphingolipids" = "Sphingolipids",
                 "Other class" = "Other", "Glycerophospholipids" = "Glycerophospholipids")

rd <- function(k) {
  f <- file.path(IN_DIR, sprintf("FC_%s_metabolites.csv", k))
  if (!file.exists(f)) stop("Missing ", f)
  z <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  z$CLASS_disp <- collapse_class(z$CLASS)
  z$sig <- !is.na(z$p) & !is.na(z$qvalue) & z$p < P_CUT & z$qvalue < Q_CUT
  notn <- trimws(as.character(z$INDIVIDUAL.NOTATION))
  z$displayable <- !(is.na(notn) | notn %in% c("", "NA"))
  z$grp <- ifelse(z$CLASS_disp == MUTED, MUTED,
           ifelse(z$CLASS_disp %in% KEYCL, z$CLASS_disp, "Other class"))
  z$grp <- factor(z$grp, levels = CLASS_ORDER)
  z
}
FC <- lapply(CONTRASTS, function(cn) rd(cn$key))
names(FC) <- vapply(CONTRASTS, `[[`, "", "key")

## Shared axes across the two volcanoes.
XR <- 1.08 * max(vapply(FC, function(z)
        max(abs(z$log2_FC[z$sig & z$displayable]), na.rm = TRUE), numeric(1)))
YR <- 1.05 * max(vapply(FC, function(z)
        max(-log10(z$p), na.rm = TRUE), numeric(1)))

## ---- panels a / b : volcano ------------------------------------------------
volcano <- function(z, tag) {
  z$y <- -log10(z$p)
  ## the P at which q first reaches 0.10, from this comparison's own data
  qline <- suppressWarnings(max(z$p[z$qvalue < Q_CUT], na.rm = TRUE))
  ns    <- z[!z$sig, ]
  sg    <- z[z$sig & z$displayable, ]
  ## The key lists only classes this panel actually plots. drop = FALSE used to
  ## keep all five levels in both volcanoes, so panel b carried a Sphingolipids
  ## swatch with no sphingolipid point behind it -- alcohol on the wild type has
  ## no significant sphingolipid. Restricting the levels here fixes it in code,
  ## so it does not have to be deleted by hand in Illustrator on every render.
  ## Colours are keyed by name in CLASS_COL, so a class that IS in both panels
  ## keeps the same colour and the two keys stay consistent where they overlap.
  keylev <- CLASS_ORDER[CLASS_ORDER %in% as.character(sg$grp)]
  sg$grp  <- factor(as.character(sg$grp),  levels = keylev)
  cat(sprintf("     panel %s key: %s\n", tag, paste(keylev, collapse = ", ")))
  lab   <- sg[as.character(sg$grp) != MUTED, ]
  lab   <- lab[order(lab$p), ][seq_len(min(N_LABEL, nrow(lab))), ]

  ggplot() +
    geom_hline(yintercept = -log10(qline), linetype = "dashed",
               linewidth = 0.25, colour = "#888888") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "#CCCCCC") +
    geom_point(data = ns, aes(log2_FC, y), colour = "#E4E4E4", size = 0.45, shape = 16) +
    geom_point(data = sg, aes(log2_FC, y, colour = grp), size = 1.5, shape = 16,
               alpha = 0.9) +
    ggrepel::geom_text_repel(data = lab, aes(log2_FC, y, label = INDIVIDUAL.NOTATION,
                                             colour = grp),
                             size = VOLC_LAB / .pt, seed = 1, max.overlaps = Inf,
                             min.segment.length = 0, segment.size = 0.15,
                             segment.colour = "#AAAAAA", box.padding = 0.28,
                             point.padding = 0.12, force = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = CLASS_COL, labels = CLASS_SHORT,
                        drop = FALSE, name = NULL) +
    ## The default 5% padding on each side is most of the white margin inside
    ## the panel, and XR/YR already carry their own headroom (8% and 5%).
    scale_x_continuous(expand = expansion(mult = 0.005)) +
    scale_y_continuous(expand = expansion(mult = 0.005)) +
    coord_cartesian(xlim = c(-XR, XR), ylim = c(0, YR)) +
    labs(x = expression("Log"[2]~"fold-change"),
         y = expression(-log[10]~italic(P))) +
    guides(colour = guide_legend(ncol = 1,
                                 override.aes = list(size = 1.8, label = "", alpha = 1))) +
    theme_nc() +
    theme(legend.position = c(0.002, 0.002), legend.justification = c(0, 0),
          legend.key.size = grid::unit(7, "pt"),
          legend.text = element_text(size = VOLC_TXT),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.margin = ggplot2::margin(1, 3, 1, 1),
          axis.title = element_text(size = VOLC_TXT),
          axis.text  = element_text(size = VOLC_TXT),
          plot.margin = ggplot2::margin(1, 1, 1, 1, "pt")) +
    nc_tag(tag)
}

## ---- panel c : significant metabolites per class, both comparisons ---------
## Counts are built on the ANNOTATED significant set (sig & displayable), the
## same set the volcanoes plot, so the figure is internally consistent: 102 and
## 27. The unannotated 5 each are counted in the text as 107 and 32 and are
## reported in the legend, not drawn.
class_counts <- function() {
  rows <- do.call(rbind, lapply(CONTRASTS, function(cn) {
    s <- FC[[cn$key]]; s <- s[s$sig & s$displayable, ]
    data.frame(class      = collapse_class(s$CLASS),
               dir        = ifelse(s$log2_FC > 0, "increased", "decreased"),
               comparison = cn$name, stringsAsFactors = FALSE)
  }))
  tot   <- table(rows$class)
  minor <- names(tot)[tot < MIN_CLASS_TOTAL]
  if (length(minor)) {
    cat(sprintf("  pooled into '%s' (fewer than %d across both comparisons):\n     %s\n",
                MINOR_LABEL, MIN_CLASS_TOTAL, paste(sort(minor), collapse = ", ")))
    rows$class[rows$class %in% minor] <- MINOR_LABEL
  }
  ## table() over three factors fills the grid, so a class absent from one
  ## comparison gets an explicit zero and both arms keep the same bar height.
  cnt <- as.data.frame(table(class = rows$class, dir = rows$dir,
                             comparison = rows$comparison), stringsAsFactors = FALSE)
  cnt$n <- cnt$Freq
  cnt$Freq[cnt$dir == "decreased"] <- -cnt$Freq[cnt$dir == "decreased"]
  cnt
}

CNT <- class_counts()
## Class order: by total in the lead comparison, largest at the top. The pooled
## row is forced to the bottom whatever its count, so it never sits between two
## real classes.
lead <- CNT[CNT$comparison == ORDER_BY, ]
ordv <- tapply(abs(lead$Freq), lead$class, sum)
lev  <- names(sort(ordv, decreasing = FALSE))
lev  <- c(MINOR_LABEL, setdiff(lev, MINOR_LABEL))
CNT$class <- factor(CNT$class, levels = lev)
## ggplot puts the FIRST factor level at the bottom of a dodged band, so the
## levels are the reverse of the reading order in PAIR_ORDER.
CNT$comparison <- factor(CNT$comparison, levels = rev(PAIR_ORDER))

XPOS <- 1.14 * max(CNT$Freq)
XNEG <- 1.14 * min(CNT$Freq)

barpanel_merged <- function(mode) {
  dodge <- position_dodge(width = DODGE_W)
  lb <- CNT[CNT$Freq != 0, ]
  lb$hj <- ifelse(lb$Freq > 0, -0.35, 1.35)
  zt <- stats::aggregate(n ~ class + comparison, data = CNT, FUN = sum)
  zt <- zt[zt$n == 0, ]

  if (identical(mode, "direction")) {
    fill_aes <- "dir"
    fill_val <- c(increased = NC_RED, decreased = NC_BLUE)
    fill_brk <- c("increased", "decreased")
    fill_lab <- ggplot2::waiver()          # "increased"/"decreased", not a contrast
    ## room at the left for the per-bar comparison names
    xlo <- XNEG - 0.46 * (XPOS - XNEG)
  } else {
    fill_aes <- "comparison"
    fill_val <- CONTRAST_COL[PAIR_ORDER]
    fill_brk <- PAIR_ORDER
    ## Reader-facing names come from nc_style.R's display layer, never from the
    ## internal keys the data carries. Parsed, for the I148M superscript.
    fill_lab <- function(x) nc_label_expr(x)
    xlo <- XNEG - 0.12 * (XPOS - XNEG)   # room for the widest decreased label
  }
  ## Breaks come from the DATA range, not from the limits: in "direction" mode
  ## the limits reach far left to make a gutter for the comparison names, and
  ## pretty() over that gutter puts tick marks where no bar can ever be.
  brks <- pretty(c(XNEG, XPOS), n = 6)
  brks <- brks[brks >= XNEG & brks <= XPOS]

  p <- ggplot(CNT, aes(Freq, class, fill = .data[[fill_aes]], group = comparison)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = NC_INK) +
    ## group = comparison is explicit on every layer. An implicit grouping from
    ## the fill aesthetic is what put Fig 2e's P values over the wrong bars.
    geom_col(width = DODGE_W * 0.86, position = dodge,
             colour = NC_INK, linewidth = 0.12) +
    geom_text(data = lb, aes(label = n, hjust = hj, group = comparison),
              position = dodge, size = NC_BIG / .pt, colour = NC_INK,
              inherit.aes = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = fill_val, breaks = fill_brk, labels = fill_lab,
                      name = NULL) +
    ## Two rows, not one: at 88 mm the two contrast names side by side are wider
    ## than the panel and the second one clips.
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(labels = abs, limits = c(xlo, XPOS), breaks = brks,
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(drop = FALSE,
                     labels = function(x) ifelse(x %in% names(CLASS_AXIS),
                                                 unname(CLASS_AXIS[x]), x)) +
    labs(x = "Significant metabolites", y = NULL) +
    theme_nc() +
    ## The key sits INSIDE the panel, bottom right. Only glycerophospholipids
    ## reaches past x = 17, so everything right of about x = 25 and below the top
    ## band is empty plotting area -- putting the key there costs nothing, where
    ## a key above the panel costs ~7 mm of height off a 58 mm panel and leaves
    ## the bars visibly shorter than panel d's heatmap beside them.
    ## If the counts ever change so that a second class runs long, move it back
    ## to "top" rather than letting it overlap a bar.
    theme(legend.position = LEGEND_POS, legend.justification = LEGEND_JUST,
          legend.key.size = grid::unit(6, "pt"),
          ## hjust = 0 — parsed labels differ in width and ggplot sizes every
          ## legend cell to the widest, so without it the short label drifts
          ## off its own swatch.
          legend.text = element_text(size = NC_BIG, hjust = 0),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.margin = ggplot2::margin(1, 2, 1, 1),
          axis.title  = element_text(size = NC_BIG),
          axis.text   = element_text(size = NC_BIG),
          plot.margin = ggplot2::margin(2, 4, 2, 2, "pt")) +
    nc_tag("c")

  if (SHOW_ZERO && nrow(zt))
    p <- p + geom_text(data = zt, aes(x = 0, y = class, label = "0",
                                      group = comparison),
                       position = dodge, hjust = -0.45, size = NC_BIG / .pt,
                       colour = NC_INK2, inherit.aes = FALSE, show.legend = FALSE)

  if (identical(mode, "direction")) {
    ## The colour channel is spent on direction here, so the comparison has to
    ## be named. Small grey text at the left of each bar, not a second legend:
    ## a legend would make the reader map a colour to a comparison and then a
    ## position to a direction, which is one indirection too many.
    pl <- unique(CNT[, c("class", "comparison")])
    ## parse = TRUE: the plotmath string carries the I148M superscript.
    pl$lab <- unname(CONTRAST_PLOTMATH[as.character(pl$comparison)])
    p <- p + geom_text(data = pl, aes(x = xlo, y = class, label = lab,
                                      group = comparison), parse = TRUE,
                       position = dodge, hjust = 0, size = NC_MIN / .pt,
                       colour = NC_INK2, inherit.aes = FALSE, show.legend = FALSE)
  }
  p
}

## ---- the old two-panel bars, kept for MERGED_BARS <- FALSE -----------------
barpanel_split <- function(tb, tag, levels_, cmax) {
  tb$class <- factor(tb$class, levels = levels_)
  ggplot(tb, aes(Freq, class, fill = dir)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = NC_INK) +
    geom_col(width = 0.68, colour = NC_INK, linewidth = 0.12) +
    geom_text(aes(label = abs(Freq), hjust = ifelse(Freq > 0, -0.35, 1.35)),
              size = NC_BIG / .pt, colour = NC_INK) +
    scale_fill_manual(values = c("increased" = NC_RED, "decreased" = NC_BLUE),
                      name = NULL, breaks = c("increased", "decreased")) +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(limits = c(-cmax, cmax), labels = abs) +
    labs(x = "Significant metabolites", y = NULL) +
    theme_nc() +
    theme(legend.position = "top", legend.key.size = grid::unit(6, "pt"),
          legend.text = element_text(size = NC_BIG),
          axis.title  = element_text(size = NC_BIG),
          axis.text   = element_text(size = NC_BIG)) +
    nc_tag(tag)
}

## ---- build + report --------------------------------------------------------
for (cn in CONTRASTS) {
  z <- FC[[cn$key]]
  save_nc(volcano(z, cn$panel),
          file.path(OUT, sprintf("Fig1%s_volcano_%s", cn$panel,
                                 tolower(gsub(" ", "_", cn$name)))), VOLC_W, VOLC_H)
  s <- z[z$sig, ]
  cat(sprintf("  %-22s %3d significant, %3d annotated  (up %d / down %d)\n",
              cn$name, nrow(s), sum(s$displayable),
              sum(s$displayable & s$log2_FC > 0), sum(s$displayable & s$log2_FC < 0)))
}

if (MERGED_BARS) {
  for (mode in COLOUR_BY)
    save_nc(barpanel_merged(mode),
            file.path(OUT, sprintf("Fig1c_classcounts_both_%s", mode)),
            BAR_W_MM, BAR_H_MM)
  cat("\n  panel c, counts as drawn\n")
  wide <- stats::reshape(CNT[, c("class", "dir", "comparison", "n")],
                         idvar = c("class", "comparison"), timevar = "dir",
                         direction = "wide")
  names(wide) <- sub("^n\\.", "", names(wide))
  wide <- wide[order(wide$class, wide$comparison), ]
  print(wide, row.names = FALSE)
} else {
  bd <- split(CNT, CNT$comparison)
  cm <- 1.28 * max(abs(CNT$Freq))
  for (cn in CONTRASTS)
    save_nc(barpanel_split(bd[[cn$name]], cn$panel, levels(CNT$class), cm),
            file.path(OUT, sprintf("Fig1_classcounts_%s", tolower(gsub(" ", "_", cn$name)))),
            88, 88)
}

src <- do.call(rbind, lapply(CONTRASTS, function(cn) {
  z <- FC[[cn$key]]
  data.frame(comparison = cn$name, metabolite = as.character(z$INDIVIDUAL.NOTATION),
             class = as.character(z$CLASS_disp), log2_fold_change = z$log2_FC,
             P_value = z$p, q_value = z$qvalue, significant = z$sig,
             stringsAsFactors = FALSE) }))
utils::write.csv(src, file.path(OUT, "SourceData_Fig1.csv"), row.names = FALSE)
cat(sprintf("\n  Source Data written: %d rows\n", nrow(src)))
