# =============================================================================
#  Figure 2, panels e and f — phosphatidylcholine saturation
#  Nature Communications style. Answers Dr Sanyal's C13/C14 (revision task B1).
#
#  WHAT IS BEING COMPARED, AND WHAT THE P VALUE MEANS
#  --------------------------------------------------
#  The claim is that the alcohol effect DEPENDS ON SATURATION. So the unit of
#  analysis is the mouse, and the quantity is one number per mouse:
#
#      index = log2( total PC with 0-1 double bonds / total PC with >=4 )
#
#  The P value is the group difference in that index. It is grouped on
#  saturation and compared across contrasts, which is the claim itself.
#
#  This replaces the earlier species-level tests, which had two faults:
#    1. the unit was the species, not the mouse. Species within a class share an
#       acyl-CoA pool and a compositional constraint, so treating 14 of them as
#       14 independent observations is pseudoreplication and the P values were
#       far too small.
#    2. "the 0-1 bin is significant and the >=4 bin is not" is the
#       difference-in-significance fallacy — it does not test that the two
#       differ. The index does.
#
#  Species with 2-3 double bonds are PLOTTED but EXCLUDED from the index: at
#  sum-composition resolution such a species may carry two monounsaturated
#  chains or one 18:2, so it is genuinely ambiguous. One threshold, used
#  everywhere: <=1 versus >=4.
#
#  e  mean log2 fold-change of each saturation bin, grouped by contrast.
#     DESCRIPTIVE ONLY — no P value is printed here, deliberately. The index
#     test uses the 0-1 and >=4 bins and ignores the 2-3 bin, so a P centred
#     over four bars would imply it covers all four when it covers two. It is
#     also the same number panel f already carries. All four exact P values go
#     in the legend. Animal level; error bars are s.e. across mice.
#  f  the index itself, one point per mouse, with the group comparisons.
#
#  Exported separately at final size for assembly: e 120 mm, f 60 mm.
# =============================================================================

ROOT      <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT   <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
DATA_FILE <- Sys.getenv("OWL_XLSX",   file.path(ROOT, "Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx"))
OUT       <- Sys.getenv("FIG2_OUT",   file.path(FIGROOT, "Fig 2/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

## ---- load (OWL_CACHE is an optional CSV copy of the two sheets) -------------
CACHE <- Sys.getenv("OWL_CACHE", "")
if (nzchar(CACHE)) {
  d  <- utils::read.csv(file.path(CACHE, "owl_data.csv"),  check.names = FALSE, stringsAsFactors = FALSE)
  cl <- utils::read.csv(file.path(CACHE, "owl_class.csv"), check.names = FALSE, stringsAsFactors = FALSE)
} else {
  suppressPackageStartupMessages(library(readxl))
  d  <- readxl::read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 16)
  cl <- readxl::read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 0, n_max = 16)
}
names(d) <- names(cl)
cm <- local({ x <- cl[, -c(2:4)]; t2 <- t(x[-1]); colnames(t2) <- x[["COD.OWL"]]; data.frame(t2) })
mets <- names(d)[-c(1:4)]
for (m in mets) d[[m]][is.na(d[[m]])] <- min(d[[m]], na.rm = TRUE) / 5

## G is the KEY used by CON below — keep it free of line breaks.
## G_LAB is the two-line version used only for the x axis of panel f.
GRP_MAP <- c("AAV Luc" = "WD", "PNPLA3 MUT" = "WD+var",
             "PNPLA3 Luc-WD + EtOH" = "WD+alc", "PNPLA3 MUT + EtOH" = "WD+var+alc")
G_WRAP  <- c("WD" = "WD", "WD+var" = "WD\n+var",
             "WD+alc" = "WD\n+alc", "WD+var+alc" = "WD+var\n+alc")
d$G <- unname(GRP_MAP[as.character(d$Sample)])
d   <- d[!is.na(d$G), ]
d$G <- factor(d$G, levels = unname(GRP_MAP))

## ---- phosphatidylcholine species, binned by total double bonds -------------
pc <- intersect(rownames(cm)[cm$SIMPLIFIED.NAME == "PC_DAPC"], mets)
db <- suppressWarnings(as.integer(sub(".*\\((\\d+):(\\d+)\\).*", "\\2", cm[pc, "INDIVIDUAL.NOTATION"])))
pc <- pc[!is.na(db)]; db <- db[!is.na(db)]
BINS <- c("0-1", "2-3", "4-5", "6+")
bin  <- cut(db, c(-1, 1, 3, 5, 99), labels = BINS)
M    <- as.matrix(d[, pc]); storage.mode(M) <- "double"
cat(sprintf("  %d PC species; per bin: %s\n", length(pc),
            paste(sprintf("%s n=%d", BINS, as.integer(table(bin))), collapse = ", ")))

tot   <- sapply(BINS, function(b) rowSums(M[, bin == b, drop = FALSE]))
index <- log2(rowSums(M[, db <= 1, drop = FALSE]) / rowSums(M[, db >= 4, drop = FALSE]))

## These names are the panel e x-axis labels, so they are sentence case and
## match nc_style.R's CONTRAST_LABEL exactly. The two lookups below use them as
## keys -- change a name here and both must change with it.
CON <- list("Alcohol on wild type" = c("WD",     "WD+alc"),
            "Alcohol on variant"   = c("WD+var", "WD+var+alc"),
            "Variant alone"        = c("WD",     "WD+var"),
            "Variant on alcohol"   = c("WD+alc", "WD+var+alc"))

## ---- panel e: bin fold-changes, grouped by contrast ------------------------
bars <- do.call(rbind, lapply(names(CON), function(nm) {
  g <- CON[[nm]]
  do.call(rbind, lapply(BINS, function(b) {
    a <- log2(tot[d$G == g[2], b]); r <- log2(tot[d$G == g[1], b])
    data.frame(contrast = nm, ubin = b, mean = mean(a) - mean(r),
               se = sqrt(var(a)/length(a) + var(r)/length(r)),
               stringsAsFactors = FALSE) })) }))
bars$contrast <- factor(bars$contrast, levels = names(CON))
bars$ubin     <- factor(bars$ubin, levels = BINS)

## the index test, one per contrast — this is the P shown above each group
idx <- do.call(rbind, lapply(names(CON), function(nm) {
  g  <- CON[[nm]]
  a  <- index[d$G == g[2]]; r <- index[d$G == g[1]]
  tt <- t.test(a, r, var.equal = FALSE)
  data.frame(contrast = nm, diff = mean(a) - mean(r), t = unname(tt$statistic),
             df = unname(tt$parameter), p = tt$p.value,
             n_ref = length(r), n_trt = length(a), stringsAsFactors = FALSE) }))
idx$contrast <- factor(idx$contrast, levels = names(CON))

ytop <- max(bars$mean + bars$se); ybot <- min(bars$mean - bars$se)
span <- ytop - ybot
idx$lab <- sprintf("0-1 vs 6+ and 4-5:~italic(P)==%s",
                   vapply(idx$p, function(x) fmt_p_math(x), character(1)))
idx$lab <- vapply(idx$p, function(x) sprintf("italic(P)==%s", fmt_p_math(x)), character(1))

## Saturation is an ORDERED variable, so it gets a sequential ramp: one hue,
## light to dark with increasing unsaturation. The hue is purple because red and
## blue already carry meaning in this manuscript (variant and alcohol) and must
## not be reused for a third thing. Worst adjacent pair separates by dE 15.1
## under protanopia and 16.5 with normal vision, against a floor of 8.
BIN_COL <- c("0-1" = "#D9CCEC", "2-3" = "#A98CD6", "4-5" = "#7A4FBF", "6+" = "#452A80")

e <- ggplot(bars, aes(contrast, mean, fill = ubin)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NC_INK) +
  geom_col(position = position_dodge(0.8), width = 0.72,
           colour = NC_INK, linewidth = 0.15) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = ubin),
                position = position_dodge(0.8), width = 0.2,
                colour = NC_INK, linewidth = 0.25) +
  scale_fill_manual(values = BIN_COL, name = "Double bonds") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = NULL, y = expression("Mean log"[2]~"fold-change")) +
  guides(fill = guide_legend(nrow = 1, title.position = "left",
                             keywidth  = grid::unit(7, "pt"),
                             keyheight = grid::unit(7, "pt"),
                             override.aes = list(colour = NC_INK, linewidth = 0.15))) +
  theme_nc(grid_y = TRUE) +
  theme(legend.position = "top",
        legend.title = element_text(size = NC_BASE),
        axis.text.x = element_text(size = NC_BASE)) +
  nc_tag("e")

save_nc(e, file.path(OUT, "Fig2e_PC_unsaturation"), width_mm = 120, height_mm = 62)

## ---- panel f: the index itself, one point per mouse ------------------------
fd <- data.frame(G = factor(unname(G_WRAP[as.character(d$G)]), levels = unname(G_WRAP)),
                 alc = grepl("alc", as.character(d$G)), index = index)
mu <- do.call(rbind, lapply(split(fd, fd$G), function(z)
  data.frame(G = z$G[1], m = mean(z$index), se = sd(z$index)/sqrt(nrow(z)))))
ytop2 <- max(fd$index); span2 <- diff(range(fd$index))
p_wt  <- idx$p[idx$contrast == "Alcohol on wild type"]
p_mut <- idx$p[idx$contrast == "Alcohol on variant"]

f <- ggplot(fd, aes(G, index)) +
  geom_errorbar(data = mu, aes(G, y = m, ymin = m - se, ymax = m + se),
                width = 0.18, linewidth = 0.25, colour = NC_INK, inherit.aes = FALSE) +
  geom_crossbar(data = mu, aes(G, y = m, ymin = m, ymax = m), width = 0.42,
                linewidth = 0.3, colour = NC_INK, fatten = 0, inherit.aes = FALSE) +
  geom_point(aes(fill = alc), shape = 21, size = 1.3, stroke = 0.2,
             colour = NC_INK, position = position_jitter(width = 0.11, height = 0, seed = 1)) +
  annotate("segment", x = 1, xend = 3, y = ytop2 + span2 * 0.16, yend = ytop2 + span2 * 0.16,
           linewidth = 0.25, colour = NC_INK) +
  annotate("text", x = 2, y = ytop2 + span2 * 0.21, parse = TRUE,
           label = fmt_p_expr(p_wt), size = NC_P / .pt, colour = NC_INK) +
  annotate("segment", x = 2, xend = 4, y = ytop2 + span2 * 0.03, yend = ytop2 + span2 * 0.03,
           linewidth = 0.25, colour = NC_INK) +
  annotate("text", x = 3, y = ytop2 + span2 * 0.08, parse = TRUE,
           label = fmt_p_expr(p_mut), size = NC_P / .pt, colour = NC_INK) +
  scale_fill_manual(values = c("FALSE" = NC_GREY, "TRUE" = NC_BLUE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  ## "4+", not ">=4": the plotmath >= is a real Unicode glyph and drops out as a
  ## missing-character box in Illustrator, the same failure as Fig 1c's pool
  ## labels. It also matches panel e's own bin names, which are already 4-5/6+.
  labs(x = NULL, y = expression("Log"[2]~"(PC 0-1 db / PC 4+ db)")) +
  theme_nc() +
  nc_tag("f")

save_nc(f, file.path(OUT, "Fig2f_PC_saturation_index"), width_mm = 60, height_mm = 62)

## ---- Source Data + the numbers for the legend ------------------------------
utils::write.csv(data.frame(sample_id = d[[1]], group = as.character(d$G), tot,
                            saturation_index = index, check.names = FALSE),
                 file.path(OUT, "SourceData_Fig2e-f.csv"), row.names = FALSE)
utils::write.csv(bars, file.path(OUT, "SourceData_Fig2e_bars.csv"), row.names = FALSE)

cat("\n-- panel e/f statistics, all two-tailed Welch t-tests on mice --\n")
for (i in seq_len(nrow(idx)))
  cat(sprintf("  %-22s diff = %+5.2f  t = %5.2f  d.f. = %4.1f  P = %-8s  n = %d vs %d\n",
      as.character(idx$contrast[i]), idx$diff[i], idx$t[i], idx$df[i],
      fmt_p(idx$p[i]), idx$n_ref[i], idx$n_trt[i]))
