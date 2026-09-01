# =============================================================================
#  Figure 4 / supplement — pathview export, organised and named
#
#  WHAT THIS REPLACES
#  Section 5 of fig4_consolidated_pathview_boxplots.R wrote maps into the middle
#  of fig4_outs/ as "mmu00590.Q3_6hep.png", beside the KEGG downloads
#  ("mmu00590.png", "mmu00590.xml") and every CSV the pipeline produces. Nobody
#  could tell which file was which map, or which was a download. This does only
#  the export, into folders, with names that say what they are.
#
#  WHAT CHANGED SCIENTIFICALLY -- ONE THING
#  gene.data now comes from the ANIMAL-LEVEL limma, not the population-stacked
#  model. The two give identical log2 fold-changes (verified to two decimals
#  across the 24 Cyps), so no map changes colour. The point is that Figure 4
#  then rests on one model end to end and the Methods does not have to describe
#  two. The replication-unit problem was never in these colours -- it was in the
#  FDR, and the FDR is not what colours a pathview.
#
#  Metabolite colours are read from Figure 1's own table, so the maps and
#  Figure 1 cannot disagree about a fold-change. (See HANDOFF.md 8: the amino
#  acid panel now reads P and q from that same file, for the same reason.)
#
#  A NOTE ON SVG, BECAUSE IT KEEPS COMING UP
#  kegg.native = TRUE draws onto KEGG's own PNG of the map. That background is a
#  raster, so the output is a PNG at KEGG's fixed resolution and NO vector
#  export exists -- which is why these were unreadable at 180 mm. The only
#  vector route is kegg.native = FALSE, which re-lays the pathway out with
#  Graphviz: fully editable text, but it does not look like the KEGG map a
#  reader recognizes. This script writes BOTH so the trade can be judged by
#  looking. For the supplement the fix is one map per page at full size.
# =============================================================================

## ---------------------------------------------------------------------------
##  TUNING
## ---------------------------------------------------------------------------
QUESTION   <- "Q3"
CELL_LABEL <- "hepatocytes"   # goes in every filename; change if the input does
FDR_PATH   <- 0.10            # pathway selection, on FDR_Pmax from the integration
LIMIT      <- 3               # colour limit, genes and compounds
CEX        <- 0.25
WRITE_NATIVE   <- TRUE        # KEGG-style PNG, raster, familiar
WRITE_GRAPHVIZ <- TRUE        # Graphviz PDF, vector, editable text, different look
## ---------------------------------------------------------------------------

ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",   file.path(ROOT, "writing skeleton/figure making"))
OUTS    <- Sys.getenv("FIG4_OUTS", file.path(FIGROOT, "Fig 4/fig4_outs"))
MET_CSV <- Sys.getenv("FIG4_MET",  file.path(FIGROOT, "Fig 1",
                      sprintf("FC_%s_metabolites.csv", QUESTION)))
PMAX_CSV <- Sys.getenv("FIG4_PMAX", file.path(ROOT, "GB_metabolomics_20260628",
                      "Q3_StoufferPMAX_pooled", "modalities_pmax_gene_plus_metDec12.csv"))
EXPORT  <- file.path(OUTS, "pathview_export")

suppressPackageStartupMessages({
  library(pathview); library(org.Mm.eg.db); library(AnnotationDbi) })

.nc <- file.path(FIGROOT, "nc_style.R")
if (file.exists(.nc)) source(.nc) else {
  NC_RED_TINT <- "#E59E9E"; NC_BLUE_TINT <- "#95BAE2" }

DIR_PNG <- file.path(EXPORT, "native_png")
DIR_PDF <- file.path(EXPORT, "vector_pdf")
DIR_CAC <- file.path(EXPORT, "kegg_cache")
for (d in c(DIR_PNG, DIR_PDF, DIR_CAC)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- gene.data : animal-level log2 fold-change, ENTREZ-named ---------------
af <- file.path(OUTS, sprintf("animal_level_limma_%s.csv", QUESTION))
if (!file.exists(af)) stop("Run fig4_animal_level_limma.R first -- missing ", af)
an <- utils::read.csv(af, stringsAsFactors = FALSE)
map <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(an$SYMBOL),
                             keytype = "SYMBOL", columns = "ENTREZID")
map <- map[!is.na(map$ENTREZID) & !duplicated(map$SYMBOL), ]
ge  <- merge(an[, c("SYMBOL", "logFC_animal")], map, by = "SYMBOL")
gene.data <- stats::setNames(ge$logFC_animal, ge$ENTREZID)
gene.data <- gene.data[!is.na(names(gene.data)) & !is.na(gene.data)]
cat(sprintf("  gene.data  : %d ENTREZ-mapped genes (animal-level logFC)\n", length(gene.data)))
stopifnot("gene.data empty -- ENTREZ mapping failed" = length(gene.data) > 0)

## ---- cpd.data : metabolite log2 fold-change, from Figure 1's own table -----
cpd.data <- NULL
if (file.exists(MET_CSV)) {
  m <- utils::read.csv(MET_CSV, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("KEGG", "log2_FC") %in% names(m))) {
    k <- trimws(as.character(m$KEGG))
    ok <- !is.na(k) & nzchar(k) & k != "NA" & !is.na(m$log2_FC)
    ## a metabolite can carry several KEGG ids; take the first, and average any
    ## id that ends up claimed by more than one measured species
    k1 <- sub("[;, ].*$", "", k[ok])
    cpd.data <- tapply(m$log2_FC[ok], k1, mean)
    cpd.data <- stats::setNames(as.numeric(cpd.data), names(cpd.data))
    cat(sprintf("  cpd.data   : %d KEGG compounds (from %s)\n",
                length(cpd.data), basename(MET_CSV)))
  }
}
if (is.null(cpd.data)) cat("  cpd.data   : none -- maps will be drawn from genes only\n")

## ---- pathway list ----------------------------------------------------------
if (!file.exists(PMAX_CSV)) stop("Pathway table not found: ", PMAX_CSV,
                                 "\n  set FIG4_PMAX to the current one")
sig <- utils::read.csv(PMAX_CSV, stringsAsFactors = FALSE)
sig <- sig[!is.na(sig$FDR_Pmax) & sig$FDR_Pmax < FDR_PATH &
             sig$concordance != "One-missing", ]
sig <- sig[order(sig$FDR_Pmax), ]
cat(sprintf("  pathways   : %d at FDR_Pmax < %.2f\n\n", nrow(sig), FDR_PATH))
stopifnot("No pathways selected" = nrow(sig) > 0)

slug <- function(x) {
  x <- gsub("\\s*-\\s*Mus musculus.*$", "", x)
  x <- gsub("[^A-Za-z0-9]+", "-", x)
  gsub("^-|-$", "", substr(x, 1, 60))
}

## pathview always writes into the working directory and re-downloads the KEGG
## files if they are not there, so run it inside the cache folder and move the
## finished maps out. on.exit guarantees the working directory is restored even
## if a map fails.
old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(DIR_CAC)

IDX <- list()
for (i in seq_len(nrow(sig))) {
  kid <- sig$KEGG_id[i]; pid <- sub("^mmu", "", kid); ttl <- sig$KEGG_title[i]
  base <- sprintf("Fig4_pathview_%s_%s_%s_%s", kid, slug(ttl), CELL_LABEL, QUESTION)
  row <- list(KEGG_id = kid, KEGG_title = ttl, FDR_Pmax = sig$FDR_Pmax[i],
              concordance = sig$concordance[i], native_png = NA_character_,
              vector_pdf = NA_character_)

  run <- function(native) {
    sfx <- if (native) "native" else "graphviz"
    args <- list(gene.data = gene.data, pathway.id = pid, species = "mmu",
                 out.suffix = sfx, kegg.native = native, same.layer = !native,
                 cex = CEX, limit = list(gene = LIMIT, cpd = LIMIT),
                 low  = list(gene = NC_BLUE_TINT, cpd = NC_BLUE_TINT),
                 mid  = list(gene = "#FFFFFF",    cpd = "#FFFFFF"),
                 high = list(gene = NC_RED_TINT,  cpd = NC_RED_TINT))
    if (!is.null(cpd.data)) args$cpd.data <- cpd.data
    tryCatch({ do.call(pathview, args); TRUE },
             error = function(e) { message("  ", kid, " ", sfx, " FAILED: ",
                                           conditionMessage(e)); FALSE })
  }

  if (WRITE_NATIVE && run(TRUE)) {
    src <- sprintf("mmu%s.native.png", pid)
    dst <- file.path(DIR_PNG, paste0(base, ".png"))
    if (file.exists(src) && file.rename(src, dst)) row$native_png <- basename(dst)
  }
  if (WRITE_GRAPHVIZ && run(FALSE)) {
    src <- sprintf("mmu%s.graphviz.pdf", pid)
    dst <- file.path(DIR_PDF, paste0(base, ".pdf"))
    if (file.exists(src) && file.rename(src, dst)) row$vector_pdf <- basename(dst)
  }
  IDX[[length(IDX) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
  cat(sprintf("  %-10s %-46s %s %s\n", kid, substr(slug(ttl), 1, 46),
              if (!is.na(row$native_png)) "png" else "---",
              if (!is.na(row$vector_pdf)) "pdf" else "---"))
}
setwd(old)

idx <- do.call(rbind, IDX)
utils::write.csv(idx, file.path(EXPORT, "INDEX_pathviews.csv"), row.names = FALSE)

cat(sprintf("\n  written to %s\n", EXPORT))
cat("    native_png/   KEGG-style maps, raster (place one per supplement page)\n")
cat("    vector_pdf/   Graphviz layout, editable text, different look\n")
cat("    kegg_cache/   KEGG downloads -- regenerable, do NOT ship to GitHub\n")
cat("    INDEX_pathviews.csv   pathway id, title, FDR_Pmax, and both filenames\n")
cat(sprintf("\n  %d of %d maps have a raster version, %d have a vector version\n",
            sum(!is.na(idx$native_png)), nrow(idx), sum(!is.na(idx$vector_pdf))))
