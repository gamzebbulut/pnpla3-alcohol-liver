# =============================================================================
#  Combined lipid supplement, panel b — unsaturation slopes, both pools,
#  all four contrasts
#
#  REPLACES the mirror scatter + slope forest pair. Those were two orthogonal
#  slices of one fit -- two pools on one contrast, and one pool on four -- and
#  they overlapped in exactly one number. Fitting the whole grid instead is
#  fewer panels AND more information.
#
#  WHAT IT SHOWS THAT THE PAIR COULD NOT
#  The phospholipid slope is not flat everywhere. Under alcohol on the wild type
#  it is significantly negative (-0.084, P = 0.0013): BOTH pools lose their
#  polyunsaturated members, which is what peroxidation of polyunsaturated acyl
#  chains looks like wherever those chains sit. Under the variant, only
#  triglyceride moves. So pool-specificity is itself a property of the
#  mechanism, not a fixed feature -- and the old pair hid that by showing the
#  phospholipid comparison on one contrast only, the contrast where it happens
#  to be flat.
#
#  Point and 95% interval, not bars: a bar draws a rectangle from zero to the
#  estimate and asserts the distance from zero is the quantity of interest.
#  Here the quantity is a slope that sits either side of zero and the interval
#  matters more than the bar.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
W_MM <- 88
H_MM <- 58
TG_CLASSES <- "TG"
PL_CLASSES <- c("PC","PE","PI","PS","PG","LPC","LPE","LPI","LPG","LPS")
QLAB <- c(Q4 = "Alcohol on wild type", Q3 = "Alcohol on variant",
          Q2 = "Variant alone",       Q5 = "Variant on alcohol")
## Pool colours are the Figure 1 class palette, NOT red/blue -- red means
## increased and blue decreased everywhere in this manuscript, and this panel's
## x axis already carries a direction.
POOL_COL <- c("Triglyceride" = "#E65100", "Phospholipid" = "#555555")
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("SUPP_DE",   file.path(FIGROOT, "Fig 3/Lipid_LSEA_Q3_vs_Q4",
                                             "de_lipids_Q2_Q5.csv"))
OUT     <- Sys.getenv("SUPP_OUT",  file.path(FIGROOT, "supp_lipid_story/nc_panels"))

source(file.path(FIGROOT, "nc_style.R"))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nc_style_report()

if (!file.exists(IN)) stop("Missing ", IN)
d <- utils::read.csv(IN, stringsAsFactors = FALSE)
d$db    <- suppressWarnings(as.numeric(d$total_cs))
d$logFC <- suppressWarnings(as.numeric(d$logFC))
d$arm   <- unname(QLAB[sub(" .*", "", d$contrast)])
d$pool  <- ifelse(d$Class %in% TG_CLASSES, "Triglyceride",
           ifelse(d$Class %in% PL_CLASSES, "Phospholipid", NA))
d <- d[!is.na(d$db) & !is.na(d$logFC) & !is.na(d$arm) & !is.na(d$pool), ]
stopifnot("no rows after filtering" = nrow(d) > 0)

ARMS <- unname(QLAB[c("Q4","Q3","Q2","Q5")])
sl <- do.call(rbind, lapply(ARMS, function(a) do.call(rbind, lapply(names(POOL_COL), function(pl) {
  s <- d[d$arm == a & d$pool == pl, ]
  if (nrow(s) < 3) return(NULL)
  co <- summary(stats::lm(logFC ~ db, data = s))$coefficients["db", ]
  data.frame(arm = a, pool = pl, slope = co[1], se = co[2], p = co[4], n = nrow(s),
             stringsAsFactors = FALSE) }))))
## Reading order runs top to bottom; the first factor level sits at the bottom
## of a discrete y axis, so the levels are reversed.
sl$arm  <- factor(sl$arm,  levels = rev(ARMS))
sl$pool <- factor(sl$pool, levels = names(POOL_COL))

dodge <- position_dodge(width = 0.55)
p <- ggplot(sl, aes(slope, arm, colour = pool, group = pool)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = NC_INK) +
  geom_errorbarh(aes(xmin = slope - 1.96 * se, xmax = slope + 1.96 * se),
                 height = 0, linewidth = 0.5, position = dodge) +
  geom_point(size = 1.8, shape = 16, position = dodge) +
  scale_colour_manual(values = POOL_COL, name = NULL) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  labs(x = "Slope of log2 fold-change\non acyl-chain double bonds", y = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8))) +
  theme_nc() +
  theme(legend.position = "top", legend.justification = "left",
        legend.key.size = grid::unit(6, "pt"),
        legend.text = element_text(size = NC_BIG),
        axis.title = element_text(size = NC_BIG),
        axis.text  = element_text(size = NC_BIG),
        plot.margin = ggplot2::margin(2, 4, 2, 2, "pt")) +
  nc_tag("b")

save_nc(p, file.path(OUT, "Supp_lipid_b_slopes"), W_MM, H_MM)
utils::write.csv(sl, file.path(OUT, "SourceData_Supp_lipid_b_slopes.csv"), row.names = FALSE)

cat(sprintf("\n  %-22s %-14s %8s %10s %8s %5s\n",
            "contrast", "pool", "slope", "95% CI", "P", "n"))
for (i in seq_len(nrow(sl)))
  cat(sprintf("  %-22s %-14s %+8.3f  %+.3f/%+.3f %8s %5d\n",
              sl$arm[i], sl$pool[i], sl$slope[i],
              sl$slope[i] - 1.96 * sl$se[i], sl$slope[i] + 1.96 * sl$se[i],
              fmt_p(sl$p[i]), sl$n[i]))
cat("\n  LEGEND MATERIAL\n")
cat("  (1) The triglyceride slope flips sign by mechanism: negative under\n")
cat("      alcohol, positive under the variant, with intervals nowhere near\n")
cat("      overlapping between the two.\n")
cat("  (2) The phospholipid slope is negative under alcohol on the WILD TYPE\n")
cat("      but not on the variant. Do not describe phospholipids as 'flat' --\n")
cat("      they are flat on three contrasts of four, and the exception is the\n")
cat("      informative one.\n")
