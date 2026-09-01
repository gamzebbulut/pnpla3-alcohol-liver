# =============================================================================
# PC log2FC resolved by total unsaturation, 4 contrasts
#   Companion panel to suppfig3/TG_by_unsaturation_4contrasts.R -- same input,
#   same binning, same statistics, same visual grammar, PC instead of TG.
#
#   Answers Dr. Sanyal's C13/C14: phosphatidylcholine is hepatoprotective, so
#   does alcohol expand the pool or change its composition? Composition. The
#   rise is carried by the saturated/monounsaturated species; the
#   polyunsaturated ones do not move.
#
#   HEADLINE: mean log2FC of the 0-1 double-bond species minus that of the
#   >=4 double-bond species, per contrast (Welch t on per-species logFC).
#
# Input: Fig 3/Lipid_LSEA_Q3_vs_Q4/de_lipids_Q2_Q5.csv  -- lipidr de_analysis,
#        limma moderated statistics, isomeric duplicates already averaged.
#        NB this is why there are 43 PC species and not 55: the provider file
#        lists 12 of them twice as unresolved isomers.
# =============================================================================

ROOT <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova")
IN   <- file.path(ROOT, "writing skeleton/figure making/Fig 3/Lipid_LSEA_Q3_vs_Q4/de_lipids_Q2_Q5.csv")
OUT  <- file.path(ROOT, "writing skeleton/figure making/Fig 2/pc_saturation_outs")
## optional: the raw workbook enables the per-animal confirmation at the end
DATA_FILE <- file.path(ROOT, "Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx")
FONT <- "Arial"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## contrast keys -> labels and colours, identical to the TG panel
qlab <- c("Q2"="C3 mutation","Q5"="C4 mut (alc bg)",
          "Q3"="C1 alcohol (mut)","Q4"="C2 alcohol (WT)")
qord <- c("C3 mutation","C4 mut (alc bg)","C1 alcohol (mut)","C2 alcohol (WT)")
cols <- c("C3 mutation"="#C62828","C4 mut (alc bg)"="#E57373",
          "C1 alcohol (mut)"="#1565C0","C2 alcohol (WT)"="#64B5F6")

## ---- load, keep PC, bin exactly as the TG panel does -----------------------
d <- utils::read.csv(IN, stringsAsFactors = FALSE)
d <- d[d$Class == "PC", ]
d$contrastC <- unname(qlab[sub(" .*", "", d$contrast)])
d$db <- suppressWarnings(as.integer(d$total_cs))
d$cl <- suppressWarnings(as.integer(d$total_cl))
d$logFC <- suppressWarnings(as.numeric(d$logFC))
d <- d[!is.na(d$db) & !is.na(d$logFC) & !is.na(d$contrastC), ]
d$ubin <- cut(d$db, breaks = c(-1, 1, 3, 5, Inf), labels = c("0-1","2-3","4-5","6+"))
bins <- levels(d$ubin)

cat(sprintf("PC species: %d unique across %d contrasts (max %d double bonds)\n",
            length(unique(d$Molecule)), length(unique(d$contrastC)), max(d$db)))

se <- function(x) stats::sd(x) / sqrt(length(x))

## ---- per-bin means, and whether each bin differs from zero -----------------
agg <- do.call(rbind, lapply(qord, function(cn) do.call(rbind, lapply(bins, function(b) {
  v <- d$logFC[d$contrastC == cn & d$ubin == b]
  p <- if (length(v) > 1) stats::t.test(v)$p.value else NA_real_
  data.frame(contrastC = cn, ubin = b, n = length(v), mean = mean(v), se = se(v),
             p_vs0 = p, sig = ifelse(!is.na(p) & p < 0.05, "*", ""),
             stringsAsFactors = FALSE)
}))))
utils::write.csv(agg, file.path(OUT, "pc_by_unsaturation_bins.csv"), row.names = FALSE)

## ---- HEADLINE: saturated/MUFA minus polyunsaturated ------------------------
head_stat <- do.call(rbind, lapply(qord, function(cn) {
  lo <- d$logFC[d$contrastC == cn & d$db <= 1]
  hi <- d$logFC[d$contrastC == cn & d$db >= 4]
  tt <- stats::t.test(lo, hi, var.equal = FALSE)
  data.frame(contrastC = cn, n_0_1 = length(lo), n_4up = length(hi),
             mean_0_1 = mean(lo), mean_4up = mean(hi),
             difference = mean(lo) - mean(hi),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
             stringsAsFactors = FALSE)
}))
utils::write.csv(head_stat, file.path(OUT, "pc_saturated_vs_polyunsaturated.csv"), row.names = FALSE)
utils::write.csv(d[, c("Molecule","contrastC","db","cl","logFC","t","P.Value","adj.P.Val")],
                 file.path(OUT, "pc_by_species_lipidr.csv"), row.names = FALSE)

cat("\n== per-bin mean log2FC ==\n");           print(agg, digits = 3, row.names = FALSE)
cat("\n== HEADLINE: 0-1 db minus >=4 db ==\n"); print(head_stat, digits = 3, row.names = FALSE)

## ---- optional per-animal confirmation, independent of everything above -----
## Share of the PC pool carried by species with >=3 double bonds. >=3 guarantees
## at least one polyunsaturated chain; a total of 2 could be two MUFA chains.
## One number per mouse, so n is 19 animals rather than 43 species.
if (file.exists(DATA_FILE) && requireNamespace("readxl", quietly = TRUE)) {
  raw <- readxl::read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 16)
  cls <- readxl::read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 0, n_max = 16)
  names(raw) <- names(cls)
  cm <- cls[, -c(2:4)]; tmp <- t(cm[-1]); colnames(tmp) <- cm[["COD.OWL"]]; cm <- data.frame(tmp)
  mets <- names(raw)[-c(1:4)]
  for (m in mets) raw[[m]][is.na(raw[[m]])] <- min(raw[[m]], na.rm = TRUE) / 5
  g <- c("AAV Luc"="HF", "PNPLA3 MUT"="HF + PNPLA3",
         "PNPLA3 MUT + EtOH"="HF+PNPLA3+ALC", "PNPLA3 Luc-WD + EtOH"="HF+ALC")
  grp  <- unname(g[as.character(raw$Sample)])
  notn <- stats::setNames(as.character(cm$INDIVIDUAL.NOTATION), rownames(cm))
  simp <- stats::setNames(as.character(cm$SIMPLIFIED.NAME),     rownames(cm))
  pcc  <- intersect(names(simp)[simp == "PC_DAPC"], mets)
  m    <- regmatches(notn[pcc], regexpr("\\(\\d+:\\d+\\)", notn[pcc]))
  dbv  <- rep(NA_integer_, length(pcc))
  dbv[nzchar(m)] <- as.integer(sub(".*:(\\d+)\\)", "\\1", m[nzchar(m)]))
  pcc <- pcc[!is.na(dbv)]; dbv <- dbv[!is.na(dbv)]
  keep <- !is.na(grp)
  M <- as.matrix(raw[keep, pcc, drop = FALSE]); gg <- grp[keep]
  share <- rowSums(M[, dbv >= 3, drop = FALSE]) / rowSums(M)
  ## per-animal values, one row per mouse -- panel f of Figure 2 plots these
  utils::write.csv(data.frame(sample_id = raw[[1]][keep], group = gg,
                              pufa_share_pct = 100 * share),
                   file.path(OUT, "pc_pufa_share_values.csv"), row.names = FALSE)
  CN <- list("C3 mutation"      = c("HF + PNPLA3","HF"),
             "C4 mut (alc bg)"  = c("HF+PNPLA3+ALC","HF+ALC"),
             "C1 alcohol (mut)" = c("HF+PNPLA3+ALC","HF + PNPLA3"),
             "C2 alcohol (WT)"  = c("HF+ALC","HF"))
  sh <- do.call(rbind, lapply(qord, function(cn) {
    a <- share[gg == CN[[cn]][1]]; b <- share[gg == CN[[cn]][2]]
    data.frame(contrastC = cn, n_before = length(b), n_after = length(a),
               before_pct = 100*mean(b), after_pct = 100*mean(a),
               change_pts = 100*(mean(a) - mean(b)),
               p = stats::t.test(a, b, var.equal = FALSE)$p.value,
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(sh, file.path(OUT, "pc_pufa_share_per_animal.csv"), row.names = FALSE)
  cat("\n== per-animal confirmation: share of PC pool with >=3 double bonds ==\n")
  print(sh, digits = 3, row.names = FALSE)
} else {
  cat("\n(raw workbook not found - skipped the per-animal share confirmation)\n")
}

## ---- figure, mirroring the TG panel ----------------------------------------
nl <- agg[agg$contrastC == qord[1], c("ubin","n")]
agg$lab <- factor(paste0(agg$ubin, "\n(n=", nl$n[match(agg$ubin, nl$ubin)], ")"),
                  levels = paste0(nl$ubin, "\n(n=", nl$n, ")"))
agg$contrastC <- factor(agg$contrastC, levels = qord)
sub <- paste0("saturated/MUFA minus polyunsaturated: ",
              paste(sprintf("%s %+.2f (p=%.3f)", sub(" .*", "", head_stat$contrastC),
                            head_stat$difference, head_stat$p), collapse = "   ·   "))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- ggplot(agg, aes(lab, mean, fill = contrastC)) +
    geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey25", linewidth = 0.2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.8), width = 0.2, colour = "grey25") +
    geom_text(aes(label = sig, y = mean + sign(mean) * (se + 0.06)),
              position = position_dodge(0.8), size = 5, colour = "grey25", vjust = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
    scale_fill_manual(values = cols, name = NULL, breaks = qord) +
    labs(x = "Phosphatidylcholine total acyl-chain double bonds",
         y = expression(mean~log[2]~fold-change),
         title = "Alcohol raises phosphatidylcholine only in the saturated species",
         subtitle = sub) +
    theme_classic(base_size = 12, base_family = FONT) +
    theme(legend.position = "top", plot.title.position = "plot",
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
  okpdf <- TRUE
  tryCatch(ggsave(file.path(OUT, "PC_by_unsaturation_4contrasts.pdf"), p,
                  width = 7.4, height = 4.8, device = grDevices::cairo_pdf),
           error = function(e) okpdf <<- FALSE)
  if (!okpdf) ggsave(file.path(OUT, "PC_by_unsaturation_4contrasts.pdf"), p, width = 7.4, height = 4.8)
  ggsave(file.path(OUT, "PC_by_unsaturation_4contrasts.png"), p, width = 7.4, height = 4.8, dpi = 300)
} else {
  mat <- matrix(agg$mean, nrow = length(qord), byrow = TRUE, dimnames = list(qord, bins))
  sem <- matrix(agg$se,   nrow = length(qord), byrow = TRUE)
  draw <- function() {
    op <- par(mar = c(4.6, 4.6, 4.4, 1), xpd = NA)
    bp <- barplot(mat, beside = TRUE, col = cols[qord], border = "grey25",
                  ylim = range(c(mat - sem, mat + sem, 0)) * 1.35,
                  xlab = "Phosphatidylcholine total acyl-chain double bonds",
                  ylab = "mean log2 fold-change")
    arrows(bp, mat - sem, bp, mat + sem, angle = 90, code = 3, length = .02, col = "grey25")
    abline(h = 0, col = "grey30")
    text(bp, mat + sign(mat) * (sem + .05), matrix(agg$sig, nrow = length(qord), byrow = TRUE), cex = 1.1)
    title("Alcohol raises phosphatidylcholine only in the saturated species", adj = 0, cex.main = 1.05)
    mtext(sub, side = 3, line = .3, adj = 0, cex = .62, col = "grey30")
    legend("topright", legend = qord, fill = cols[qord], border = "grey25", bty = "n", cex = .75)
    par(op)
  }
  png(file.path(OUT, "PC_by_unsaturation_4contrasts.png"), width = 2220, height = 1440, res = 300); draw(); dev.off()
  pdf(file.path(OUT, "PC_by_unsaturation_4contrasts.pdf"), width = 7.4, height = 4.8);            draw(); dev.off()
}
cat(sprintf("\nwrote PC_by_unsaturation_4contrasts.{pdf,png} + csv files to\n  %s\n", OUT))
