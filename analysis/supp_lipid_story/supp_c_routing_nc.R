# =============================================================================
#  Combined lipid supplement, panel c — lipid routing
#
#  WHERE DOES THE FAT GO. Panels a, b and d establish that polyunsaturated
#  storage lipid is retained under the variant and depleted under alcohol. This
#  panel asks the complementary question: what is BUILT. Pools are class crossed
#  with saturation, the resolution a class-level view cannot reach -- "TG went
#  down" and "polyunsaturated TG went down while saturated TG did not" are
#  different claims and only the second is true.
#
#  BARS, AND WHY NOT DUMBBELLS (tried and rejected 30 Aug)
#  The four contrasts are two natural pairs, so paired dumbbells looked like the
#  compact answer. They were not. A dumbbell means "this quantity moved from
#  here to there" -- it carries a before-and-after semantics. These two ends are
#  not a movement; they are two independent experiments on two backgrounds, and
#  drawing a connector between them asserts a transition that was never
#  measured. It also read badly: the alcohol pair and the variant pair sit on
#  opposite sides of zero, so the two dumbbells of one pool ended up far apart
#  and the panel read as ten scattered objects rather than five pools.
#
#  Grouped bars from zero say exactly what is true and nothing more: length is
#  magnitude, side is direction, adjacency is "same pool". They need the page
#  width to do it, which is why this panel is a full-width row.
#
#  Restyled 30 Aug from suppfig3/lipid_routing_panel.R. Pools, means and
#  standard errors are unchanged; C1-C4 labels and the bar geometry are gone.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
## 88 x 58 -- c shares the middle row with the slope forest. The bars are thin
## at this width but the panel only needs five group labels, and the alternative
## (another full-width band) makes the page read as a stack of ribbons.
W_MM   <- 88
H_MM   <- 58
QLAB <- c(Q4 = "Alcohol on wild type", Q3 = "Alcohol on variant",
          Q2 = "Variant alone",       Q5 = "Variant on alcohol")
POOLS <- list(
  list(key = "TG 6+",  compartment = "Storage",  cls = "TG",  test = function(db) db >= 6),
  list(key = "TG 0-1", compartment = "Storage",  cls = "TG",  test = function(db) db <= 1),
  list(key = "PC 0-1", compartment = "Membrane", cls = "PC",  test = function(db) db <= 1),
  list(key = "PE 0-1", compartment = "Membrane", cls = "PE",  test = function(db) db <= 1),
  list(key = "LPC",    compartment = "Membrane", cls = "LPC", test = function(db) rep(TRUE, length(db))))
## ASCII only in these keys -- the first draft used a >= sign and an en dash and
## both rendered as missing-glyph boxes. Same failure as Figure 5's triangles.
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
d$db    <- suppressWarnings(as.integer(d$total_cs))
d$logFC <- suppressWarnings(as.numeric(d$logFC))
d$arm   <- unname(QLAB[sub(" .*", "", d$contrast)])
d <- d[!is.na(d$db) & !is.na(d$logFC) & !is.na(d$arm), ]

d$pool <- NA_character_; d$compartment <- NA_character_
for (P in POOLS) {
  hit <- d$Class == P$cls & P$test(d$db)
  if (any(hit & !is.na(d$pool)))
    warning("species matched two pools: ",
            paste(unique(d$Molecule[hit & !is.na(d$pool)]), collapse = ", "))
  d$pool[hit] <- P$key; d$compartment[hit] <- P$compartment
}
d <- d[!is.na(d$pool), ]
stopifnot("no species fell into any pool" = nrow(d) > 0)

KEYS <- vapply(POOLS, `[[`, "", "key")
COMP <- vapply(POOLS, `[[`, "", "compartment")
se <- function(x) stats::sd(x) / sqrt(length(x))
g <- interaction(d$pool, d$arm, drop = TRUE)
agg <- do.call(rbind, lapply(split(seq_len(nrow(d)), g), function(i) data.frame(
  pool = d$pool[i][1], arm = d$arm[i][1], mean = mean(d$logFC[i]),
  se = se(d$logFC[i]), n = length(i), stringsAsFactors = FALSE)))

agg$pool <- factor(agg$pool, levels = KEYS)
agg$compartment <- factor(unname(stats::setNames(COMP, KEYS)[as.character(agg$pool)]),
                          levels = c("Storage", "Membrane"))
agg$arm <- factor(agg$arm, levels = unname(QLAB[c("Q4","Q3","Q2","Q5")]))

dodge <- position_dodge(width = 0.78)
p <- ggplot(agg, aes(pool, mean, fill = arm, group = arm)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NC_INK) +
  geom_col(position = dodge, width = 0.72, colour = NC_INK, linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = dodge,
                width = 0.18, linewidth = 0.25, colour = NC_INK) +
  facet_grid(~ compartment, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = CONTRAST_COL, name = NULL,
                    breaks = unname(QLAB[c("Q4","Q3","Q2","Q5")])) +
  labs(x = NULL, y = expression("Mean log"[2]~"fold-change")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_nc(grid_y = TRUE) +
  theme(legend.position = "top", legend.justification = "left",
        legend.key.size = grid::unit(6, "pt"),
        legend.text = element_text(size = NC_BIG),
        strip.text = element_text(size = NC_BIG),
        panel.spacing = grid::unit(5, "pt"),
        axis.title = element_text(size = NC_BIG),
        axis.text  = element_text(size = NC_BIG)) +
  nc_tag("c")

save_nc(p, file.path(OUT, "Supp_lipid_c_routing"), W_MM, H_MM)
utils::write.csv(agg[, c("pool","arm","mean","se","n")],
                 file.path(OUT, "SourceData_Supp_lipid_c_routing.csv"), row.names = FALSE)

cat("\n  species per pool\n")
for (k in KEYS) cat(sprintf("    %-8s n = %d\n", k, agg$n[agg$pool == k][1]))
cat("\n  LEGEND MATERIAL: alcohol depletes polyunsaturated storage and builds\n")
cat("  saturated membrane and lyso species; the variant does the storage half\n")
cat("  without the membrane half. Give n per pool -- TG 6+ (n = 5) and PE 0-1\n")
cat("  (n = 4) are small and carry the largest effects on the panel. Error bars\n")
cat("  are standard errors ACROSS SPECIES within a pool, not across animals.\n")
