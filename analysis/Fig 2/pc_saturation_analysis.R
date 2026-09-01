## ════════════════════════════════════════════════════════════════════
##  PHOSPHATIDYLCHOLINE BY SATURATION  —  answers Dr. Sanyal's C13/C14
##
##  C13: "PC is generally considered hepatoprotective."
##  C14: "Have we looked at this carefully — most folk think this is
##        decreased in alcohol. Does the transcriptomic data support
##        less utilization or more synthesis?"
##
##  The manuscript currently reports "total phosphatidylcholine" going
##  up, inside a panel labelled as lipids with hepatotoxic activity.
##  This splits the diacyl-PC species by degree of unsaturation and
##  asks whether the pool expands or whether its composition shifts.
##
##  Three outputs:
##    1. per-species table   log2FC, p, q for all 4 contrasts
##    2. group table         SFA / MUFA / PUFA index per contrast
##    3. trend table         slope of log2FC per double bond  <- the
##                           single statistic that answers the question
##
##  HOW TO USE: edit the PATHS block, then Source.
## ════════════════════════════════════════════════════════════════════

## ── PATHS ───────────────────────────────────────────────────────────
DATA_FILE <- "~/Desktop/Lipidomics project with Dr Smirnova/Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx"
OUT_DIR   <- "~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 2/pc_saturation_outs"
CSV_FALLBACK <- ""          # tidy per-sample csv: cod, group, sample_id, value

PC_CLASS  <- "PC_DAPC"      # diacyl-PC (SIMPLIFIED.NAME). Ether species excluded.
PUFA_MIN  <- 2              # >= this many double bonds counts as polyunsaturated
ACCENT    <- "#1565C0"      # house blue, matches volcano_standalone.R

## ════════════════════════════════════════════════════════════════════
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load_from_excel <- function(path) {
  stopifnot(requireNamespace("readxl", quietly = TRUE))
  data <- readxl::read_excel(path, sheet = "Data per metabolite", skip = 16)
  cls  <- readxl::read_excel(path, sheet = "Data per metabolite", skip = 0, n_max = 16)
  names(data) <- names(cls)
  cm <- cls[, -c(2:4)]; tmp <- t(cm[-1]); colnames(tmp) <- cm[["COD.OWL"]]; cm <- data.frame(tmp)
  mets <- names(data)[-c(1:4)]
  for (m in mets) data[[m]][is.na(data[[m]])] <- min(data[[m]], na.rm = TRUE) / 5
  g <- c("Male 24w CD"="CD","AAV Luc"="HF","PNPLA3 MUT"="HF + PNPLA3",
         "PNPLA3 MUT + EtOH"="HF+PNPLA3+ALC","PNPLA3 Luc-WD + EtOH"="HF+ALC")
  list(mat = as.matrix(data[, mets]), group = unname(g[as.character(data$Sample)]),
       simp = setNames(as.character(cm$SIMPLIFIED.NAME), rownames(cm)),
       notn = setNames(as.character(cm$INDIVIDUAL.NOTATION), rownames(cm)))
}
load_from_csv <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  w <- reshape(d[, c("cod","group","sample_id","value")],
               idvar = c("group","sample_id"), timevar = "cod", direction = "wide")
  mets <- sub("^value\\.", "", grep("^value\\.", names(w), value = TRUE))
  mat <- as.matrix(w[, grep("^value\\.", names(w))]); colnames(mat) <- mets
  ## impute exactly as the Excel path does, so the two sources agree
  for (j in seq_len(ncol(mat)))
    mat[is.na(mat[, j]), j] <- min(mat[, j], na.rm = TRUE) / 5
  map <- c(CTL_NAL="HF", MUT_NAL="HF + PNPLA3", MUT_PAL="HF+PNPLA3+ALC", CTL_PAL="HF+ALC")
  a <- unique(d[, c("cod","simp","notation")]); rownames(a) <- a$cod
  list(mat = mat, group = unname(map[w$group]),
       simp = setNames(a$simp, a$cod), notn = setNames(a$notation, a$cod))
}
src <- if (nzchar(CSV_FALLBACK) && file.exists(CSV_FALLBACK)) load_from_csv(CSV_FALLBACK) else
       load_from_excel(path.expand(DATA_FILE))

## ── select diacyl-PC and read total double bonds from the notation ──
db_of <- function(x) {
  m <- regmatches(x, regexpr("\\(\\d+:\\d+\\)", x))
  out <- rep(NA_integer_, length(x))
  ok <- nzchar(m) & !is.na(m)
  out[which(ok)] <- as.integer(sub(".*:(\\d+)\\)", "\\1", m[ok]))
  out
}
cods <- names(src$simp)[src$simp %in% PC_CLASS]
cods <- cods[cods %in% colnames(src$mat)]
DB   <- db_of(src$notn[cods]); names(DB) <- cods
cods <- cods[!is.na(DB)]; DB <- DB[cods]

## NB: do NOT drop species that have a missing value in some other group --
## a species can be fully measured in the two groups of the contrast at hand.
## Missing values are imputed upstream (min/5), matching the rest of the pipeline.
M <- src$mat[, cods, drop = FALSE]
n_imp <- colSums(!is.finite(M))
if (any(n_imp > 0))
  cat("species carrying imputed values:",
      paste(sprintf("%s(%d)", cods[n_imp > 0], n_imp[n_imp > 0]), collapse = ", "), "\n")
drop <- !is.finite(colSums(M))
if (any(drop)) { cat("unusable:", paste(cods[drop], collapse = ", "), "\n")
                 M <- M[, !drop, drop = FALSE]; cods <- cods[!drop]; DB <- DB[!drop] }
Y <- log2(M)
grp <- src$group
cat(sprintf("diacyl-PC species analysed: %d  (double bonds %d-%d)\n", ncol(Y), min(DB), max(DB)))

CONTRASTS <- list(
  "C1 alcohol on mutant"  = c("HF+PNPLA3+ALC", "HF + PNPLA3"),
  "C2 alcohol on WT"      = c("HF+ALC",        "HF"),
  "C3 variant alone"      = c("HF + PNPLA3",   "HF"),
  "C4 variant on alcohol" = c("HF+PNPLA3+ALC", "HF+ALC"))

welch <- function(a, b) stats::t.test(a, b, var.equal = FALSE)$p.value

## ── 1. per species ──────────────────────────────────────────────────
per <- do.call(rbind, lapply(names(CONTRASTS), function(cn) {
  tr <- CONTRASTS[[cn]][1]; rf <- CONTRASTS[[cn]][2]
  i <- grp == tr; j <- grp == rf
  data.frame(contrast = cn, cod = cods, notation = src$notn[cods], double_bonds = DB,
             log2FC = colMeans(Y[i, , drop=FALSE]) - colMeans(Y[j, , drop=FALSE]),
             p = vapply(seq_along(cods), function(k) welch(Y[i, k], Y[j, k]), 0),
             n_imputed_in_contrast = vapply(seq_along(cods), function(k)
               sum(!is.finite(src$mat[i | j, cods[k]])), 0L),
             row.names = NULL, stringsAsFactors = FALSE)
}))
per$q <- stats::p.adjust(per$p, "BH")
utils::write.csv(per, file.path(OUT_DIR, "pc_by_species.csv"), row.names = FALSE)

## ── 2. saturation groups (index = mean of member species) ───────────
grp_def <- list("SFA (0 db)" = DB == 0, "MUFA (1 db)" = DB == 1,
                "PUFA (>=2 db)" = DB >= PUFA_MIN)
grpres <- do.call(rbind, lapply(names(CONTRASTS), function(cn) {
  tr <- CONTRASTS[[cn]][1]; rf <- CONTRASTS[[cn]][2]
  do.call(rbind, lapply(names(grp_def), function(gn) {
    sel <- grp_def[[gn]]
    v <- rowMeans(2^Y[, sel, drop = FALSE])          # index on the linear scale
    a <- v[grp == tr]; b <- v[grp == rf]
    data.frame(contrast = cn, group = gn, n_species = sum(sel),
               fold_change = mean(a)/mean(b), p = welch(a, b),
               row.names = NULL)
  }))
}))
utils::write.csv(grpres, file.path(OUT_DIR, "pc_by_saturation_group.csv"), row.names = FALSE)

## ── 3. the trend: does log2FC depend on degree of unsaturation? ─────
trend <- do.call(rbind, lapply(names(CONTRASTS), function(cn) {
  s <- per[per$contrast == cn, ]
  f <- summary(stats::lm(log2FC ~ double_bonds, data = s))$coefficients["double_bonds", ]
  data.frame(contrast = cn, slope_per_double_bond = f[["Estimate"]], se = f[["Std. Error"]],
             t = f[["t value"]], p = f[["Pr(>|t|)"]], row.names = NULL)
}))
utils::write.csv(trend, file.path(OUT_DIR, "pc_saturation_trend.csv"), row.names = FALSE)

cat("\n== saturation groups ==\n");  print(grpres, digits = 3, row.names = FALSE)
cat("\n== trend (log2FC per double bond) ==\n"); print(trend, digits = 3, row.names = FALSE)

## ── figure: log2FC vs double bonds, one panel per contrast ──────────
## Position on x already encodes saturation, so the marks stay one ink
## colour — colouring them by a binned x would be redundant.
lab <- setNames(sprintf("slope %+.3f per double bond, p = %.3f",
                        trend$slope_per_double_bond, trend$p), trend$contrast)
per$panel <- paste0(per$contrast, "\n", lab[per$contrast])

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- ggplot(per, aes(double_bonds, log2FC)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = .4) +
    geom_smooth(method = "lm", formula = y ~ x, colour = ACCENT,
                fill = ACCENT, alpha = .13, linewidth = .8) +
    geom_point(shape = 21, size = 2.1, stroke = .45,
               colour = "grey25", fill = "grey80", alpha = .85) +
    facet_wrap(~ panel, nrow = 2) +
    scale_x_continuous(breaks = seq(0, max(per$double_bonds), 2)) +
    labs(x = "total double bonds per phosphatidylcholine species",
         y = expression(log[2]~fold~change),
         title = "Phosphatidylcholine composition shifts with saturation",
         subtitle = sprintf("%d diacyl-PC species; a negative slope means alcohol favours the saturated species",
                            length(cods))) +
    theme_classic(base_size = 12) +
    theme(plot.title.position = "plot",
          plot.subtitle = element_text(size = 9.5, colour = "grey30"),
          strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold", lineheight = 1.15),
          panel.spacing = unit(1.1, "lines"))
  ggsave(file.path(OUT_DIR, "pc_saturation.png"), p, width = 8.0, height = 6.4, dpi = 300)
  ggsave(file.path(OUT_DIR, "pc_saturation.pdf"), p, width = 8.0, height = 6.4)
} else {
  png(file.path(OUT_DIR, "pc_saturation.png"), width = 2400, height = 1900, res = 300)
  op <- par(mfrow = c(2,2), mar = c(4,4.2,3.2,1))
  for (cn in names(CONTRASTS)) {
    s <- per[per$contrast == cn, ]
    plot(s$double_bonds, s$log2FC, pch = 21, bg = "grey80", col = "grey25", bty = "l",
         xlab = "total double bonds", ylab = "log2 fold change", main = cn, cex.main = .95)
    abline(h = 0, col = "grey55"); abline(stats::lm(log2FC ~ double_bonds, data = s), col = ACCENT, lwd = 2)
    mtext(lab[[cn]], side = 3, line = .2, cex = .62)
  }
  par(op); dev.off()
  pdf(file.path(OUT_DIR, "pc_saturation.pdf"), width = 8.0, height = 6.4)
  op <- par(mfrow = c(2,2), mar = c(4,4.2,3.2,1))
  for (cn in names(CONTRASTS)) {
    s <- per[per$contrast == cn, ]
    plot(s$double_bonds, s$log2FC, pch = 21, bg = "grey80", col = "grey25", bty = "l",
         xlab = "total double bonds", ylab = "log2 fold change", main = cn, cex.main = .95)
    abline(h = 0, col = "grey55"); abline(stats::lm(log2FC ~ double_bonds, data = s), col = ACCENT, lwd = 2)
    mtext(lab[[cn]], side = 3, line = .2, cex = .62)
  }
  par(op); dev.off()
}
cat(sprintf("\nWrote to %s\n  pc_by_species.csv\n  pc_by_saturation_group.csv\n  pc_saturation_trend.csv\n  pc_saturation.png/.pdf\n", OUT_DIR))
