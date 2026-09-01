# =============================================================================
# Lipid Set Enrichment Analysis (lipidr) — alcohol response, mutant vs WT background
#   Q3 = MUT_PAL - MUT_NAL   (alcohol on PNPLA3-mutant background, the lead)
#   Q4 = CTL_PAL - CTL_NAL   (alcohol on WT/Luc background)
# Enriches the lipidome by CLASS, total CHAIN LENGTH, and UNSATURATION — the axes
# KEGG cannot see. Reads the OWL matrix (todos.xlsx), reformats OWL lipid names to
# lipidr shorthand (with a free-fatty-acid rescue so the PUFA signal is included).
# =============================================================================
set.seed(1)
suppressWarnings(suppressMessages({
  need <- c("lipidr","readxl","dplyr","tidyr","ggplot2","limma")
  miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
}))
if (length(miss))
  stop("Install first:\n  install.packages('BiocManager'); BiocManager::install('lipidr')\n",
       "  install.packages(c('readxl','dplyr','tidyr','ggplot2'))\nMissing: ", paste(miss, collapse=", "))
library(lipidr); library(readxl); library(dplyr); library(tidyr); library(ggplot2)

## ----------------------------- CONFIG ---------------------------------------
XLSX    <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/Data/todos.xlsx")
SHEET   <- "Data per metabolite"
OUT_DIR <- path.expand("~/Desktop/Lipidomics project with Dr Smirnova/writing skeleton/figure making/Fig 3/outs")
MIN_SET <- 4        # min lipids per set for LSEA
NPERM   <- 10000
SIG_FDR <- 0.10     # significance threshold for highlighting/reporting (matches the integrated pipeline)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_MAP <- c("AAV Luc"="CTL_NAL", "PNPLA3 MUT"="MUT_NAL",
               "PNPLA3 MUT + EtOH"="MUT_PAL", "PNPLA3 Luc-WD + EtOH"="CTL_PAL")  # chow excluded

## ----------------------- 1. read the transposed sheet -----------------------
raw   <- as.data.frame(readxl::read_excel(XLSX, sheet = SHEET, col_names = FALSE))
labs  <- as.character(raw[[1]])                       # row labels in col 1
r_cls <- which(labs == "CLASS")[1]
r_not <- which(labs == "INDIVIDUAL NOTATION")[1]
grp_v <- as.character(raw[[3]])                        # col 3 = Sample/group
samp  <- which(grp_v %in% names(GROUP_MAP))           # sample rows (the 4 WSD groups -> n=19)
mcol  <- 4:ncol(raw)                                  # metabolite columns

classes  <- as.character(unlist(raw[r_cls, mcol]))
notation <- as.character(unlist(raw[r_not, mcol]))
vals     <- suppressWarnings(sapply(raw[samp, mcol], function(x) as.numeric(as.character(x))))
mat      <- t(vals)                                   # metabolites x samples
colnames(mat) <- paste0("S", raw[samp, 2])
grp      <- unname(GROUP_MAP[grp_v[samp]])
cat(sprintf("Samples: %d  (%s)\n", length(grp), paste(names(table(grp)), table(grp), collapse=", ")))

## ------------- 2. OWL name -> lipidr shorthand  (+ FA rescue) ----------------
convert_name <- function(nm, cls) {
  if (is.na(nm)) return(NA_character_)
  nm <- trimws(nm)
  # free fatty acids: OWL writes them with no class prefix (e.g. "18:2n-6", "16:0")
  if (grepl("Fatty acid", cls, ignore.case = TRUE) && grepl("^[0-9]+:[0-9]+", nm)) {
    core <- trimws(sub("n-?[0-9]+.*$", "", nm))
    return(if (grepl("^[0-9]+:[0-9]+$", core)) paste("FA", core) else NA_character_)
  }
  if (grepl("\\+", nm))          return(NA_character_)   # sum species (can't parse)
  if (tolower(nm) == "cholesterol") return(NA_character_)
  m <- regmatches(nm, regexec("^([A-Za-z]+)\\((.*)\\)$", nm))[[1]]
  if (length(m) < 3) return(NA_character_)
  cl <- m[2]; inner <- m[3]
  inner <- gsub("n-?[0-9]+", "", inner)                 # strip omega position
  inner <- gsub("([0-9]+:[0-9]+)e", "O-\\1", inner)      # ether 'e' suffix -> 'O-' prefix
  lyso <- grepl("0:0/|/0:0", inner)
  if (lyso) { p <- strsplit(inner, "/")[[1]]; p <- p[p != "0:0"]; if (length(p)) inner <- paste(p, collapse="/") }
  cl <- c(ChoE="CE", CMH="HexCer", AC="CAR")[cl] %>% (function(x) if (is.na(x)) m[2] else unname(x))
  if (lyso && cl %in% c("PC","PE","PG","PI","PS","PA")) cl <- paste0("L", cl)
  paste(cl, inner)
}
lipid_cls <- grepl("phospho|glycero|sphingo|sterol|fatty", classes, ignore.case = TRUE)
clean <- mapply(convert_name, notation, classes, USE.NAMES = FALSE)
keep  <- !is.na(clean)
cat(sprintf("Metabolite columns: %d | lipid-class: %d | mapped to lipidr: %d\n",
            length(clean), sum(lipid_cls), sum(keep)))
drop_lip <- notation[lipid_cls & !keep]
cat(sprintf("Dropped lipid species (%d): %s\n", length(drop_lip),
            paste(drop_lip, collapse = ", ")))

mat        <- mat[keep, , drop = FALSE]
clean_keep <- clean[keep]
mat        <- log2(mat)                                # OWL is linear & all > 0 -> log2 fold changes
# Several OWL entries map to the same lipidr name (e.g. sn-isomeric lyso species).
# Average them into ONE row per unique lipid -- no make.unique suffix, so names stay parseable.
mat        <- limma::avereps(mat, ID = clean_keep)
clean      <- rownames(mat)
cat(sprintf("Unique lipid species after merging duplicates: %d\n", nrow(mat)))

# manual class/length/unsat annotation (backup, also handy for sanity checks)
tot <- t(sapply(clean, function(s) {
  pr <- regmatches(s, gregexpr("([0-9]+):([0-9]+)", s))[[1]]
  cc <- sum(as.numeric(sub(":.*","",pr))); db <- sum(as.numeric(sub(".*:","",pr)))
  c(length = cc, unsat = db)
}))
annot_manual <- data.frame(Molecule = clean, Class = sub(" .*","",clean),
                           total_length = tot[,"length"], total_unsat = tot[,"unsat"])
readr::write_csv(annot_manual, file.path(OUT_DIR, "lipid_annotation_manual.csv"))

## ------------------------- 3. build lipidr object ---------------------------
le <- as_lipidomics_experiment(mat)   # matrix: rownames = lipid molecule names, columns = samples
cat(sprintf("lipidr experiment built: %d lipids x %d samples\n", nrow(le), ncol(le)))
cat("  samples in experiment: ", paste(colnames(le), collapse = ", "), "\n")
# align Group to whatever samples lipidr kept: by name first, positional fallback
g <- grp[match(colnames(le), colnames(mat))]
if (anyNA(g) && ncol(le) == length(grp)) g <- grp          # names altered but order intact
if (anyNA(g)) stop("Could not map Group to all samples (see sample names printed above).")
SummarizedExperiment::colData(le)$Group <- g
# NOTE: as_lipidomics_experiment() already annotated class/chain/unsaturation into rowData.
rd        <- as.data.frame(SummarizedExperiment::rowData(le))
cat("\nrowData columns: ", paste(names(rd), collapse = ", "), "\n")
class_col <- intersect(c("Class","class","lipid_class"), names(rd))[1]
if (!is.na(class_col)) {
  cat("lipidr parsed classes:\n"); print(table(rd[[class_col]], useNA = "ifany"))
  not_parsed <- sum(is.na(rd[[class_col]]))
  if (not_parsed > 0) message(not_parsed, " molecules unparsed by lipidr (mostly ether/plasmalogen; ",
                              "kept in matrix; use lipid_annotation_manual.csv for a fallback fgsea).")
}

## --------------------------- 4. differential -------------------------------
meas <- SummarizedExperiment::assayNames(le)[1]
de   <- de_analysis(le,
                    MUT_NAL - CTL_NAL,   # Q2  mutation, no alcohol
                    MUT_PAL - MUT_NAL,   # Q3  alcohol on mutant background (lead)
                    CTL_PAL - CTL_NAL,   # Q4  alcohol on WT background
                    MUT_PAL - CTL_PAL,   # Q5  mutation on alcohol background
                    measure = meas, group_col = "Group")
desc  <- c("MUT_NAL - CTL_NAL" = "Q2 (mutation | no-alcohol)",
           "MUT_PAL - MUT_NAL" = "Q3 (alcohol | mutant bg)",
           "CTL_PAL - CTL_NAL" = "Q4 (alcohol | WT bg)",
           "MUT_PAL - CTL_PAL" = "Q5 (mutation | alcohol bg)")
de$contrast <- dplyr::recode(de$contrast, !!!desc)
short <- c("Q2 (mutation | no-alcohol)"="Q2", "Q3 (alcohol | mutant bg)"="Q3",
           "Q4 (alcohol | WT bg)"="Q4", "Q5 (mutation | alcohol bg)"="Q5")
readr::write_csv(de, file.path(OUT_DIR, "de_lipids_Q2_Q5.csv"))

## ----------------------------- 5. LSEA -------------------------------------
enr <- lsea(de, rank.by = "logFC", min_size = MIN_SET, nperm = NPERM)
readr::write_csv(as.data.frame(enr), file.path(OUT_DIR, "lsea_results_Q2_Q5.csv"))
cat("\nSignificant lipid sets per contrast (FDR <", SIG_FDR, "):\n")
print(significant_lipidsets(enr, p.cutoff = SIG_FDR))

# per-contrast enrichment plots (class + unsaturation)
for (cn in unique(de$contrast)) {
  sub <- de[de$contrast == cn, ]
  sig <- tryCatch(significant_lipidsets(enr[enr$contrast == cn, ], p.cutoff = SIG_FDR),
                  error = function(e) NULL)
  for (an in c("class","unsat","length")) {
    p <- tryCatch(plot_enrichment(sub, sig, annotation = an), error = function(e) NULL)
    if (!is.null(p)) {
      f <- file.path(OUT_DIR, sprintf("enrich_%s_%s", gsub("[^A-Za-z0-9]+","_",cn), an))
      ggsave(paste0(f,".pdf"), p, width = 7, height = 5)
      ggsave(paste0(f,".png"), p, width = 7, height = 5, dpi = 300)
    }
  }
}

## --------------- 6. cross-contrast comparison across Q2-Q5 ------------------
ew       <- as.data.frame(enr)
ew$Q     <- short[ew$contrast]
stat_col <- intersect(c("NES","Enrichment","enrichment","logFC"), names(ew))[1]
set_col  <- intersect(c("set","Set","pathway","lipid_set"), names(ew))[1]
pcol     <- intersect(c("padj","FDR","adj.P.Val","p.adj","pval"), names(ew))[1]

tidy <- ew %>% dplyr::transmute(set = .data[[set_col]], Q,
                                NES = .data[[stat_col]], padj = .data[[pcol]])
wide <- tidy %>% tidyr::pivot_wider(names_from = Q, values_from = c(NES, padj))
readr::write_csv(wide, file.path(OUT_DIR, "lsea_ALLQ_comparison.csv"))   # all four, side by side

# (a) lipid-CLASS enrichment heatmap: classes x Q2..Q5, fill = NES, * = FDR < SIG_FDR
cls <- tidy %>% dplyr::filter(grepl("^Class_", set)) %>%
  dplyr::mutate(set  = sub("^Class_", "", set),
                Q    = factor(Q, levels = c("Q2","Q3","Q4","Q5")),
                star = ifelse(!is.na(padj) & padj < SIG_FDR, "*", ""))
ord <- cls %>% dplyr::group_by(set) %>% dplyr::summarise(m = mean(NES, na.rm = TRUE)) %>%
  dplyr::arrange(m) %>% dplyr::pull(set)
cls$set <- factor(cls$set, levels = ord)
hm <- ggplot(cls, aes(Q, set, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 4.5, vjust = 0.78) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                       midpoint = 0, limits = c(-3, 3), oob = scales::squish, name = "NES") +
  labs(x = NULL, y = "Lipid class", title = "Lipid-class enrichment across contrasts (LSEA)",
       subtitle = paste0("* FDR < ", SIG_FDR,
                         "  |  Q2 mutation(no-alc) \u00b7 Q3 alcohol(mut bg) \u00b7 Q4 alcohol(WT bg) \u00b7 Q5 mutation(alc bg)")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "lsea_class_heatmap_Q2_Q5.pdf"), hm, width = 6, height = 6)
ggsave(file.path(OUT_DIR, "lsea_class_heatmap_Q2_Q5.png"), hm, width = 6, height = 6, dpi = 300)

# (b) interaction scatters: two complementary views of the same genotype x alcohol interaction
scatter_pair <- function(xq, yq, xlab, ylab, file, subtitle) {
  xs <- paste0("NES_", xq); ys <- paste0("NES_", yq)
  if (!all(c(xs, ys) %in% names(wide))) { message("skip scatter ", xq, " vs ", yq); return(invisible()) }
  d <- as.data.frame(wide); d <- d[stats::complete.cases(d[, c(xs, ys)]), ]
  d$kind <- ifelse(grepl("^Class_", d$set), "Class", "total")
  g <- ggplot(d, aes(.data[[xs]], .data[[ys]], label = set)) +
    geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey70") +
    geom_point(aes(color = kind), size = 2.6, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 2.6, max.overlaps = 20, show.legend = FALSE) +
    labs(x = xlab, y = ylab, color = "set type",
         title = "Lipid-set enrichment comparison", subtitle = subtitle) +
    theme_minimal(base_size = 11)
  ggsave(paste0(file, ".pdf"), g, width = 7.5, height = 6.5)
  ggsave(paste0(file, ".png"), g, width = 7.5, height = 6.5, dpi = 300)
}
scatter_pair("Q4", "Q3",
  "Enrichment in Q4  (alcohol, WT background)",
  "Enrichment in Q3  (alcohol, mutant background)",
  file.path(OUT_DIR, "lsea_Q3_vs_Q4_scatter"),
  "Off-diagonal = PNPLA3 modifies the ALCOHOL response")
scatter_pair("Q2", "Q5",
  "Enrichment in Q2  (mutation, no-alcohol)",
  "Enrichment in Q5  (mutation, alcohol background)",
  file.path(OUT_DIR, "lsea_Q5_vs_Q2_scatter"),
  "Off-diagonal = ALCOHOL modifies the GENOTYPE response (same interaction, other view)")

# (c) ranked lipid-set interaction (lipid-level double-difference; both views of beta4)
inter <- as.data.frame(wide) %>%
  dplyr::mutate(delta_alc = NES_Q3 - NES_Q4,     # genotype x alcohol via alcohol response (mut - WT)
                delta_gen = NES_Q5 - NES_Q2) %>%  # same interaction via mutation response (alc - no-alc)
  dplyr::select(set, NES_Q2, NES_Q3, NES_Q4, NES_Q5, delta_alc, delta_gen) %>%
  dplyr::arrange(dplyr::desc(abs(delta_alc)))
readr::write_csv(inter, file.path(OUT_DIR, "lsea_interaction_ranked.csv"))
cat("\nTop lipid sets by |Q3 - Q4| interaction:\n"); print(utils::head(inter, 12))

cat("\nDone. Outputs in:\n  ", OUT_DIR, "\n")
