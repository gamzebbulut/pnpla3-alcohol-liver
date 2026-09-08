# =============================================================================
#  nc_style.R  —  shared Nature Communications figure style
#
#  source() this at the top of every figure script:
#      source(file.path(FIGROOT, "nc_style.R"))
#
#  Specs enforced here (Nature Portfolio artwork guide, research content):
#    width      88 mm single column · 180 mm double column, hard maximum
#    height     225 mm absolute maximum, less with a long caption
#    type       sans-serif (Arial/Helvetica), 5-7 pt; panel letters 8 pt bold
#               lower-case
#    format     vector PDF with editable text, plus a 450 dpi raster preview
#    colour     RGB
#
#  Two house rules that are not Nature's but are ours:
#    1. No in-panel titles. Nature panels carry a letter, not a title; the
#       message lives in the legend. fmt_p() and the legend do that work.
#    2. No "C1..C4" anywhere a reader can see (revision task A2). Internal
#       keys are CONTRAST_LABEL / CONTRAST_SHORT; what a READER sees comes
#       from the display layer below — nc_label(), nc_label_expr(),
#       nc_labeller() — approved 8 September 2026. Never draw a raw key.
# =============================================================================

## ---- canvas -----------------------------------------------------------------
NC_W1   <- 88    # mm, single column
NC_W2   <- 180   # mm, double column
NC_HMAX <- 225   # mm, absolute page maximum

mm2in <- function(mm) mm / 25.4

## ---- type -------------------------------------------------------------------
## Arial on macOS, Helvetica next, whatever the device calls "sans" last.
## Never returns "" — svglite writes the family name straight into the SVG, and
## an empty font-family makes Illustrator fall back to its own default.
nc_font <- function() {
  fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
  if (length(fams)) {
    for (f in c("Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans")) if (f %in% fams) return(f)
    return("Helvetica")
  }
  FONT_DETECTED <<- FALSE
  "Arial"        # systemfonts absent: assume the Mac/Windows default sans
}
FONT_DETECTED <- TRUE
FONT <- nc_font()

NC_BASE <- 8   # pt — body text in panels
NC_MIN  <- 6   # pt — never go below this
## RETIRED 31 August 2026 — NOTHING USES THIS. Kept only so the archived
## scripts under Fig 4/archived/ still source cleanly.
##
## It existed for Fig 4b, whose 14 repelled pathway labels collided at 6 pt in a
## 74 mm scatter. The diagnosis was wrong: the labels were not too big, they had
## nowhere to go. Seven amino-acid pathways sit in one corner, and the y limit
## stopped just above them. Opening the y headroom (1.26 -> 1.55) and wrapping
## the labels narrower (17 -> 13 characters) placed all 14 at the 6 pt floor with
## no overlaps and no loss of width. Set at the top of that script.
##
## The lesson generalises: reach for the layout before the type size. Nature
## allows 5-7 pt and the set runs at 8/6, deliberately above the floor, because
## legibility was the repeated complaint.
NC_DENSE <- 5  # pt — retired; no current script reads this
NC_TAG  <- 8   # pt — panel letter
NC_P    <- 7   # pt — P values and significance marks. Change here and every
               #      panel follows; Nature's range is 5-7 pt, so 7 is the
               #      largest legal size.
## Figure 1 deliberately runs one step above Nature's stated 5-7 pt, because at
## 88 mm the volcanoes have empty space and legibility matters more than the
## letter of the guide. Production may rescale; if they object, set this to
## NC_BASE and re-run — nothing else changes.
NC_BIG  <- 8   # pt — axis titles, tick labels and point labels in Figure 1

## ---- colour -----------------------------------------------------------------
## Manuscript scheme, unchanged. Validated colourblind-safe: worst adjacent pair
## dE 14.8 (deuteranopia), 14.1 (tritanopia), 15.8 (normal vision).
NC_RED   <- "#C62828"; NC_RED_L  <- "#E57373"
NC_BLUE  <- "#1565C0"; NC_BLUE_L <- "#64B5F6"
NC_GREY  <- "#BBBBBB"; NC_INK    <- "#1A1A1A"; NC_INK2 <- "#555555"

## Tints of the same two hues, for fills that carry text on top — KEGG pathview
## nodes above all. Derived from NC_RED/NC_BLUE rather than hardcoded, so the
## hue family can never drift from Figures 1-3.
## Needed because Cyp4a10 and Cyp4a14 sit at 97% of the pathview colour limit:
## at full NC_BLUE their black labels are unreadable, which is the whole reason
## the original script used light fills.
nc_tint <- function(hex, frac = 0.45) {
  rgb <- grDevices::col2rgb(hex)[, 1]
  grDevices::rgb(t(round(frac * rgb + (1 - frac) * 255)), maxColorValue = 255)
}
NC_BLUE_TINT <- nc_tint(NC_BLUE)
NC_RED_TINT  <- nc_tint(NC_RED)

## KEGG pathway categories, for Figure 4a/b. Four hues plus grey — seven
## categories failed the colour-vision check (worst adjacent pair dE 1.3 under
## protanopia), so everything outside the four that carry this paper's story
## folds into "Other". Validated: worst adjacent dE 12.1 protan, 18.0 normal.
NC_CAT <- c("Amino acid metabolism"  = "#6A3FA0",
            "Lipid metabolism"       = "#E65100",
            "Cofactors and vitamins" = "#00897B",
            "Global and overview"    = "#A68B00",
            "Other"                  = "#9E9E9E")

## Plain-language contrast names. Never print C1-C4 (revision task A2).
CONTRAST_LABEL <- c(
  "C1 alcohol (mut)" = "Alcohol on variant",
  "C2 alcohol (WT)"  = "Alcohol on wild type",
  "C3 mutation"      = "Variant alone",
  "C4 mut (alc bg)"  = "Variant on alcohol")
CONTRAST_SHORT <- c(
  "C1 alcohol (mut)" = "Alcohol | variant",
  "C2 alcohol (WT)"  = "Alcohol | wild type",
  "C3 mutation"      = "Variant alone",
  "C4 mut (alc bg)"  = "variant | alcohol")
## Colour rule, approved 1 September 2026: PALE is the factor on its own,
## SATURATED is that factor on the other's background. The blue pair already
## read that way; the two reds were the wrong way round and were exchanged.
## A reader can now decode the key without reading it.
CONTRAST_COL <- c(
  "Alcohol on wild type" = NC_BLUE_L,
  "Alcohol on variant"   = NC_BLUE,
  "Variant alone"        = NC_RED_L,
  "Variant on alcohol"   = NC_RED)
## Reading order: alcohol first, then variant (revision task A1 — alcohol,
## then variant, then the combination).
CONTRAST_ORDER <- c("Alcohol on wild type", "Alcohol on variant",
                    "Variant alone", "Variant on alcohol")

relabel_contrast <- function(x, short = FALSE) {
  m <- if (short) CONTRAST_SHORT else CONTRAST_LABEL
  out <- unname(m[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

## ---- reader-facing contrast labels (approved 8 September 2026) --------------
## Everything above is INTERNAL. Scripts index by those strings and the
## published SourceData_*.csv files carry them in their `contrast` / `cmp`
## columns, so they must not move. What follows is a DISPLAY layer: the only
## thing a reader sees. A future wording change is one edit here plus a
## re-render, not twelve scripts again.
##
##   Alcohol on wild type -> Alcohol alone
##   Alcohol on variant   -> Alcohol (PNPLA3^I148M bg)
##   Variant alone        -> PNPLA3^I148M alone
##   Variant on alcohol   -> PNPLA3^I148M (alcohol bg)
##
## "Alone" is exact rather than shorthand: alcohol on the wild-type background
## IS alcohol alone, and the variant on the no-alcohol background IS the variant
## alone. The set reads as two pairs — one insult on two backgrounds, then one
## variant on two backgrounds.

## Flat form. For captions, console output and anything that cannot parse.
CONTRAST_DISPLAY <- c(
  "Alcohol on wild type" = "Alcohol alone",
  "Alcohol on variant"   = "Alcohol (PNPLA3-I148M bg)",
  "Variant alone"        = "PNPLA3-I148M alone",
  "Variant on alcohol"   = "PNPLA3-I148M (alcohol bg)")

## Plotmath form — the one that goes in figures. The superscript is
## load-bearing, not cosmetic: written flat, these labels clip off the right
## edge of an 88 mm panel (rendered evidence, notation_preview/, form D).
## plotmath is the house route for superscripts; label_parsed is already used in
## fig4_boxplots_nc.R and fig5_aminoacid_boxplots_nc.R, and ggtext was
## deliberately removed from fig5_heatmap_nc.R and is not coming back.
##
## SPACING: use *" (" — a quoted leading space. Measured at real 6 pt in Arial
## metrics: quoted space = 0.589 mm, one normal word space; `~` = 1.100 mm,
## 1.9x too wide. The `~` is what produced the wide gap in the 31 Aug preview,
## so notation_preview/render.R is NOT the form to copy.
CONTRAST_PLOTMATH <- c(
  "Alcohol on wild type" = '"Alcohol alone"',
  "Alcohol on variant"   = '"Alcohol"*" ("*PNPLA3^"I148M"*" bg)"',
  "Variant alone"        = 'PNPLA3^"I148M"*" alone"',
  "Variant on alcohol"   = 'PNPLA3^"I148M"*" (alcohol bg)"')

## Two-line form, for an axis too narrow to take the label on one line —
## Figure 5a's four columns are the only current case.
CONTRAST_PLOTMATH_2L <- c(
  "Alcohol on wild type" = 'atop("Alcohol", "alone")',
  "Alcohol on variant"   = 'atop("Alcohol", "("*PNPLA3^"I148M"*" bg)")',
  "Variant alone"        = 'atop(PNPLA3^"I148M", "alone")',
  "Variant on alcohol"   = 'atop(PNPLA3^"I148M", "(alcohol bg)")')

## Flat labels. Unknown strings pass through unchanged.
nc_label <- function(x) {
  out <- unname(CONTRAST_DISPLAY[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

## Parsed labels, one expression per input — for `labels =` in any scale, for
## annotate(parse = TRUE), and anywhere else a plotmath slot takes a vector.
## An unknown string is quoted rather than dropped, so a typo shows up as plain
## text in the panel instead of throwing.
nc_label_expr <- function(x, two_line = FALSE) {
  tab <- if (two_line) CONTRAST_PLOTMATH_2L else CONTRAST_PLOTMATH
  s <- unname(tab[as.character(x)])
  s[is.na(s)] <- sprintf('"%s"', as.character(x)[is.na(s)])
  parse(text = s)
}

## Facet labeller:  facet_wrap(~ cmp, labeller = nc_labeller())
## Built on call rather than at load, so sourcing this file never needs ggplot2.
nc_labeller <- function(two_line = FALSE) {
  tab <- if (two_line) CONTRAST_PLOTMATH_2L else CONTRAST_PLOTMATH
  ggplot2::as_labeller(tab, default = ggplot2::label_parsed)
}

## ---- group labels, as distinct from contrast labels ------------------------
## A CONTRAST is a comparison (four of them, above). A GROUP is one arm of the
## study. They were drifting apart: Fig 2d said "Variant / Variant + alcohol",
## the Fig 4d boxplots said "No alcohol / Alcohol", and Fig 2f said
## "WD / WD+var / WD+alc / WD+var+alc" -- three vocabularies for the same mice.
## Everything a reader sees now comes from here.

## Two-group boxplots (Fig 2d, the Fig 4d boxplots). BOTH groups are on the
## variant background, so naming it on every tick prints the same gene name
## twelve times across six facets, and at 8 pt in 180 mm the two names collide
## (measured: 15 mm per label, the name needs ~15 mm on its own). The
## background is therefore named ONCE, in the axis title, and the ticks carry
## only what differs.
##
## The tick marks are the ASCII hyphen and plus, NEVER U+2212. A true minus is
## the same class of character as the Unicode >= that rendered as an empty box
## in Illustrator, and as the triangles that came out as dots in Fig 5.
ALC_TICKS <- c("-", "+")

## x-axis title for those panels. background = FALSE where the surrounding
## figure has already established it.
alc_axis <- function(background = TRUE) {
  if (background) expression("Alcohol (all mice "*PNPLA3^"I148M"*")")
  else            expression("Alcohol")
}

## The four study arms (Fig 2f, and Supp Fig 1 when it is rewritten). Every arm
## is on Western diet, so "WD" carries no information inside a panel and is
## dropped from the label and stated in the legend instead: spelled out in full
## the codes truncate mid-superscript at 60 mm even at the 6 pt floor.
ARM_PLOTMATH <- c(
  "WD"         = 'atop("Control", "")',
  "WD+var"     = 'atop(PNPLA3^"I148M", "")',
  "WD+alc"     = 'atop("Alcohol", "")',
  "WD+var+alc" = 'atop(PNPLA3^"I148M", "+ alcohol")')

nc_arm_expr <- function(x) {
  s <- unname(ARM_PLOTMATH[as.character(x)])
  s[is.na(s)] <- sprintf('"%s"', as.character(x)[is.na(s)])
  parse(text = s)
}

## ---- exact P values (revision task J8) --------------------------------------
## Nature wants the number, not a threshold. Italic capital P is added by the
## caller where the grammar allows it; these return the numeral only.
##   fmt_p(0.00591)      -> "0.0059"
##   fmt_p(3.166e-05)    -> "3.2 x 10^-5"   (plotmath-ready via fmt_p_expr)
fmt_p <- function(p, digits = 2) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x >= 1e-4) return(formatC(signif(x, digits), format = "fg"))
    e <- floor(log10(x)); m <- signif(x / 10^e, digits)
    sprintf("%s x 10^%d", formatC(m, format = "fg"), e)
  }, character(1))
}
## plotmath form of a bare P value, for geom_text(parse = TRUE):
##   0.0059  ->  "0.0059"        3.166e-5 -> "3.2 %*% 10^-5"
fmt_p_math <- function(p, digits = 2) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x >= 1e-4) return(fmt_p(x, digits))
    e <- floor(log10(x)); m <- signif(x / 10^e, digits)
    sprintf("%s %%*%% 10^%d", formatC(m, format = "fg"), e)
  }, character(1))
}
## plotmath form including the italic P, for annotate(parse = TRUE)
fmt_p_expr <- function(p, digits = 2) {
  if (is.na(p)) return("italic(P)~'='~'NA'")
  if (p >= 1e-4) return(sprintf("italic(P)~'='~'%s'", fmt_p(p, digits)))
  e <- floor(log10(p)); m <- signif(p / 10^e, digits)
  sprintf("italic(P)~'='~%s %%*%% 10^%d", formatC(m, format = "fg"), e)
}
## Full test report for a legend: test type, tails, t, d.f., exact P.
fmt_test <- function(t, df, p, test = "two-tailed Welch t-test") {
  sprintf("%s, t = %.2f, d.f. = %.1f, P = %s", test, t, df, fmt_p(p))
}

## ---- theme ------------------------------------------------------------------
theme_nc <- function(base = NC_BASE, grid_y = FALSE, tick_size = NC_MIN + 1) {
  th <- ggplot2::theme_classic(base_size = base, base_family = FONT) +
    ggplot2::theme(
      text             = ggplot2::element_text(colour = NC_INK),
      axis.text        = ggplot2::element_text(size = tick_size, colour = NC_INK),
      axis.title       = ggplot2::element_text(size = base, colour = NC_INK),
      axis.line        = ggplot2::element_line(linewidth = 0.25, colour = NC_INK),
      axis.ticks       = ggplot2::element_line(linewidth = 0.25, colour = NC_INK),
      axis.ticks.length = grid::unit(1.2, "pt"),
      legend.title     = ggplot2::element_blank(),
      ## hjust = 0 is not cosmetic. ggplot sizes every legend cell to the
      ## widest label, so with parsed labels of unequal width a short one
      ## drifts away from its own key. Pinning it left keeps each label on its
      ## swatch. Any script that overrides legend.text must carry this too.
      legend.text      = ggplot2::element_text(size = base, colour = NC_INK, hjust = 0),
      legend.key.size  = grid::unit(6, "pt"),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.background = ggplot2::element_blank(),
      plot.title       = ggplot2::element_blank(),    # house rule 1: no titles
      plot.subtitle    = ggplot2::element_blank(),
      plot.caption     = ggplot2::element_text(size = NC_MIN, colour = NC_INK2, hjust = 0),
      plot.tag         = ggplot2::element_text(size = NC_TAG, face = "bold", family = FONT,
                                               colour = NC_INK, hjust = 0, vjust = 1),
      plot.tag.position = c(0, 1),
      plot.margin      = ggplot2::margin(6, 2, 2, 2, "pt"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = base, colour = NC_INK),
      panel.grid       = ggplot2::element_blank())
  if (grid_y) th <- th + ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.2, colour = "grey92"))
  th
}

## Panel letter: lower case, bold, 8 pt (Nature house style).
## NC_TAGS = FALSE (the default) exports panels WITHOUT a letter, which is what
## you want when the figure is assembled in Illustrator: letters placed there on
## one text layer align across panels, baked-in letters do not. Set NC_TAGS to
## TRUE before sourcing a figure script to get them back.
if (!exists("NC_TAGS")) NC_TAGS <- FALSE
nc_tag <- function(letter) if (isTRUE(NC_TAGS)) ggplot2::labs(tag = letter) else NULL

## ---- export -----------------------------------------------------------------
## Writes <file>.pdf (vector, text editable, fonts embedded) and <file>.png
## (450 dpi preview for slides and for eyeballing). Refuses to exceed the spec.
## SVG is the interchange format, not PDF.
##
## cairo_pdf embeds fonts as Type0/CIDFontType2 subsets with Identity encoding.
## Illustrator cannot map that encoding, so a placed cairo PDF shows the artwork
## with the text missing — which is exactly what happened on 29 Aug. This is a
## long-standing cairo/Illustrator incompatibility, not a setup problem.
##
## svglite writes real <text> elements naming the font family, so Illustrator
## renders them with your installed Arial and every label is live, editable text.
## Fallback is the base pdf() device, whose Type1/WinAnsi text Illustrator also
## reads fine — but it can only use its built-in families, so Arial becomes
## Helvetica. Nature names both, so that is a safe landing, just not identical.
save_nc <- function(plot, file, width_mm, height_mm, dpi = 450, quiet = FALSE) {
  stopifnot(width_mm  <= NC_W2 + 1e-9)
  stopifnot(height_mm <= NC_HMAX + 1e-9)
  w <- mm2in(width_mm); h <- mm2in(height_mm)
  png_file <- paste0(file, ".png")

  draw_into <- function(path, open_dev) {
    if (file.exists(path)) unlink(path)
    opened <- FALSE
    tryCatch({ open_dev(); opened <- TRUE; print(plot) },
             error = function(e) invisible(NULL),
             finally = if (opened) try(grDevices::dev.off(), silent = TRUE))
    file.exists(path) && file.info(path)$size > 500
  }

  vec_file <- paste0(file, ".svg"); how <- "svg (svglite, editable text)"
  ok <- requireNamespace("svglite", quietly = TRUE) &&
        draw_into(vec_file, function()
          svglite::svglite(vec_file, width = w, height = h, bg = "white"))
  ## svglite pins the width of every label with textLength +
  ## lengthAdjust='spacingAndGlyphs'. Illustrator honours that and re-spaces the
  ## glyphs using its own Arial metrics, which pulls words apart: "ADP" renders
  ## as "A D", "log2 FC" as "log2F C". Strip both attributes and Illustrator
  ## sets the text normally. R's own layout is unaffected; only the hint goes.
  if (ok) {
    txt <- readLines(vec_file, warn = FALSE)
    txt <- gsub(" textLength='[^']*'", "", txt)
    txt <- gsub(" lengthAdjust='[^']*'", "", txt)
    ## plotmath splits a label into several <text> runs, and some begin with a
    ## space (" = 0.038", " fold-change"). XML collapses leading whitespace, so
    ## Illustrator drops it and the runs collide: "P= 0.038". Mark the elements
    ## space-preserving AND make the leading/trailing space non-breaking, which
    ## no renderer is allowed to collapse.
    txt <- gsub("<text ", "<text xml:space='preserve' ", txt, fixed = TRUE)
    ## plotmath emits spacer runs that hold nothing but a space. svglite gives
    ## every run its own absolute x, so those spacers carry no position and are
    ## purely decorative — but Illustrator draws them as missing-glyph boxes
    ## ("log2[]fold-change"). Delete them; nothing moves.
    txt <- gsub("<text[^>]*>[[:space:]]*</text>", "", txt, useBytes = TRUE)
    writeLines(txt, vec_file, useBytes = TRUE)
  }
  if (!ok) {
    if (file.exists(vec_file)) unlink(vec_file)
    vec_file <- paste0(file, ".pdf"); how <- "pdf (base device, Helvetica not Arial)"
    ok <- draw_into(vec_file, function()
            grDevices::pdf(vec_file, width = w, height = h, useDingbats = FALSE))
  }
  if (!ok) stop("no vector file written for ", basename(file),
                " — nothing to place in Illustrator")

  ggplot2::ggsave(png_file, plot, width = w, height = h, units = "in", dpi = dpi)

  if (!quiet) {
    cat(sprintf("  %-34s %3.0f x %3.0f mm   %s\n",
                basename(vec_file), width_mm, height_mm, how))
    if (!grepl("^svg", how))
      cat("      install.packages('svglite') and re-run to keep Arial\n")
  }
  invisible(vec_file)
}

## ---- self-check -------------------------------------------------------------
nc_style_report <- function() {
  cat("nc_style.R loaded\n")
  cat(sprintf("  font        : %s%s\n", FONT,
              if (isTRUE(FONT_DETECTED)) " (found on this machine)"
              else " (assumed - install.packages('systemfonts') to verify)"))
  cat(sprintf("  widths      : %d mm single, %d mm double, %d mm max height\n", NC_W1, NC_W2, NC_HMAX))
  cat(sprintf("  type sizes  : base %d pt, floor %d pt, panel letter %d pt bold\n", NC_BASE, NC_MIN, NC_TAG))
  cat(sprintf("  P format    : %s | %s | %s\n", fmt_p(0.00591), fmt_p(3.166e-5), fmt_p(0.704)))
}
