# =============================================================================
#  Figure 2, panel d — lipid species with established hepatotoxic activity
#  Alcohol on the PNPLA3-I148M variant background. Nature Communications style.
#
#  CHANGED FROM fig2_C1_slim_boxplots.R (revision task B1):
#    - "Total PC index" REMOVED. Phosphatidylcholine is hepatoprotective, its
#      rise is not variant-specific, and it now has its own panels e and f.
#      Seven sub-panels become six.
#    - Nature style: no titles, lower-case bold panel letter, 5-7 pt Arial,
#      exact italic P (task J8), vector PDF at 180 mm.
#    - Group labels are plain language, never C1-C4 (task A2).
#
#  Exports at 180 mm so it sits directly above panels e + f (120 + 60 mm).
# =============================================================================

ROOT      <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT   <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
DATA_FILE <- Sys.getenv("OWL_XLSX",   file.path(ROOT, "Data/Data_22_0018_VCU_PNPLA3 EtOH.xlsx"))
OUT       <- Sys.getenv("FIG2_OUT",   file.path(FIGROOT, "Fig 2/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

## ---- load -------------------------------------------------------------------
## OWL_CACHE is an optional speed-up: two CSVs holding exactly what the two
## read_excel() calls below return. Unset on a normal run.
CACHE <- Sys.getenv("OWL_CACHE", "")
if (nzchar(CACHE)) {
  data          <- utils::read.csv(file.path(CACHE, "owl_data.csv"),  check.names = FALSE, stringsAsFactors = FALSE)
  metabol_class <- utils::read.csv(file.path(CACHE, "owl_class.csv"), check.names = FALSE, stringsAsFactors = FALSE)
} else {
  suppressPackageStartupMessages(library(readxl))
  data          <- readxl::read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 16)
  metabol_class <- readxl::read_excel(DATA_FILE, sheet = "Data per metabolite", skip = 0, n_max = 16)
}
names(data) <- names(metabol_class)
metabol_class.m <- local({
  cm <- metabol_class[, -c(2:4)]
  tmp <- t(cm[-1]); colnames(tmp) <- cm[["COD.OWL"]]; data.frame(tmp)
})
metabol_names <- names(data)[-c(1:4)]
for (col in metabol_names)
  data[[col]][is.na(data[[col]])] <- min(data[[col]], na.rm = TRUE) / 5
GRP_MAP <- c("Male 24w CD" = "CD", "AAV Luc" = "HF", "PNPLA3 MUT" = "HF + PNPLA3",
             "PNPLA3 MUT + EtOH" = "HF+PNPLA3+ALC", "PNPLA3 Luc-WD + EtOH" = "HF+ALC")
data$Group <- factor(unname(GRP_MAP[as.character(data$Sample)]), levels = unname(GRP_MAP))

## Scoped names. Bare REF/ALT/TRT must not be left in the global environment:
## Fig 4/fig4_consolidated_pathview_boxplots.R starts with
##   REF <- if (exists("REF")) REF else "MUT_WSD_NAL"
## so a leftover REF from another script is silently adopted and its
## makeContrasts() then fails with "object 'HF...PNPLA3' not found".
F2_REF <- "HF + PNPLA3"; F2_TRT <- "HF+PNPLA3+ALC"

## ---- the six panels (Total PC index removed) --------------------------------
simp <- setNames(metabol_class.m$SIMPLIFIED.NAME, rownames(metabol_class.m))
notn <- setNames(metabol_class.m$INDIVIDUAL.NOTATION, rownames(metabol_class.m))
ms <- function(codes) names(simp)[simp %in% codes]
mn <- function(x)     names(notn)[notn %in% x]

PANELS <- list(
  list(name = "Lyso-PC index",     cols = ms(c("LPC_1MAPC","LPC_2MAPC","LPC_1MEPC"))),
  list(name = "Free cholesterol",  cols = ms("Cho")),
  list(name = "Cholesteryl ester", cols = ms("ChoE")),
  list(name = "LPI(22:6)",         cols = mn("LPI(22:6)")),
  list(name = "LPG(22:6)",         cols = mn("LPG(22:6)")),
  list(name = "LPC(20:5)",         cols = mn("PC(0:0/20:5)")))

panel_df <- function(cols) {
  d <- data[data$Group %in% c(F2_REF, F2_TRT), c("Group", cols), drop = FALSE]
  val <- if (length(cols) == 1) d[[cols]] else rowMeans(d[, cols, drop = FALSE], na.rm = TRUE)
  data.frame(grp = ifelse(d$Group == F2_REF, "Variant", "Variant\n+ alcohol"), val = val)
}

## Long form for a faceted panel: one facet per readout, free y scale.
dfl <- do.call(rbind, lapply(PANELS, function(P) {
  z <- panel_df(P$cols); z$readout <- P$name; z }))
dfl$readout <- factor(dfl$readout, levels = vapply(PANELS, `[[`, "", "name"))
dfl$grp     <- factor(dfl$grp, levels = c("Variant", "Variant\n+ alcohol"))

## Exact two-tailed Welch t per readout, with n and d.f. (task J8).
stats <- do.call(rbind, lapply(split(dfl, dfl$readout), function(z) {
  a <- z$val[z$grp == "Variant\n+ alcohol"]; b <- z$val[z$grp == "Variant"]
  tt <- stats::t.test(a, b, var.equal = FALSE)
  rng <- range(z$val); data.frame(readout = z$readout[1], n_ref = length(b), n_trt = length(a),
    t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
    ytop = rng[2], span = diff(rng)) }))
stats$lab <- vapply(stats$p, fmt_p_expr, character(1))

d_panel <- ggplot(dfl, aes(grp, val)) +
  geom_boxplot(aes(fill = grp), width = 0.55, outlier.shape = NA,
               linewidth = 0.25, fatten = 1.2, colour = NC_INK) +
  geom_jitter(width = 0.10, height = 0, size = 0.8, shape = 21, stroke = 0.2,
              fill = "white", colour = NC_INK) +
  geom_text(data = stats, aes(x = 1.5, y = ytop + span * 0.30, label = lab),
            parse = TRUE, size = NC_P / .pt, colour = NC_INK, inherit.aes = FALSE) +
  facet_wrap(~ readout, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Variant" = NC_GREY, "Variant\n+ alcohol" = NC_BLUE),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.34))) +
  labs(x = NULL, y = "Relative abundance") +
  theme_nc() +
  nc_tag("d")

save_nc(d_panel, file.path(OUT, "Fig2d_lipotoxic"), width_mm = 180, height_mm = 42)

cat("\n-- exact statistics for the legend (all two-tailed Welch t-tests) --\n")
for (i in seq_len(nrow(stats))) cat(sprintf("  %-18s n = %d vs %d, t = %.2f, d.f. = %.1f, P = %s\n",
  as.character(stats$readout[i]), stats$n_ref[i], stats$n_trt[i],
  stats$t[i], stats$df[i], fmt_p(stats$p[i])))
