# =============================================================================
# supp_table3_joint_pathways.R
#
# Supplementary Table 3: KEGG pathways perturbed in BOTH the gene arm and the
# metabolite arm under alcohol on the PNPLA3-I148M background (14 pathways).
# This is the list Figure 4a reports and the Results cite.
#
# Rule (as stated in Fig. 4a): gene-arm Benjamini-Hochberg FDR < 0.10 AND
# metabolite-arm P < 0.05, among pathways measured in both arms.
#
# Earlier this lived in SuppTables_joint_pathways_14_and_13.docx (31 Aug), built
# by hand outside R, next to a second 13-pathway list (FDR on Pmax). That second
# list is now a column here instead of a separate table: 13 of the 14 pass
# FDR(Pmax) < 0.10; the one that does not is shaded.
#
# Input: the integration pipeline output, read verbatim -- nothing is recomputed.
#   GB_metabolomics_20260628/Q3_StoufferPMAX_pooled/modalities_pmax_gene_plus_metDec12.csv
#
# Outputs (written to OUT_DIR):
#   Supplementary_Table_3.csv / .xlsx / .docx
#
# Paths follow the repo convention: run setup.R first from a clone, or leave the
# defaults to run from the Desktop project folder.
# =============================================================================

## ---- tuning block -----------------------------------------------------------
ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
IN      <- Sys.getenv("FIG4_IN",    file.path(ROOT, "GB_metabolomics_20260628/Q3_StoufferPMAX_pooled"))
OUT_DIR <- Sys.getenv("SUPPTAB_OUT", file.path(FIGROOT, "supp_tables/out"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENE_FDR_CUT <- 0.10
MET_P_CUT    <- 0.05
PMAX_FDR_CUT <- 0.10
TABLE_NO     <- "Supplementary Table 3"
N_EXPECTED   <- 14L      # what Fig. 4a and the Results report

## UTF-8 so the minus, times and arrow glyphs survive into the Word table
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))
if (!l10n_info()$`UTF-8`) invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))

## ---- read and select --------------------------------------------------------
m <- utils::read.csv(file.path(IN, "modalities_pmax_gene_plus_metDec12.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
both_arms <- m[!is.na(m$Z_gene) & !is.na(m$Z_met), ]
hit <- both_arms[both_arms$FDR_gene < GENE_FDR_CUT & both_arms$P_met < MET_P_CUT, ]
hit <- hit[order(hit$P_met), ]
cat(sprintf("pathways measured in both arms: %d | joint hits: %d | of which FDR(Pmax) < %.2f: %d\n",
            nrow(both_arms), nrow(hit), PMAX_FDR_CUT, sum(hit$FDR_Pmax < PMAX_FDR_CUT)))
if (nrow(hit) != N_EXPECTED)
  stop(sprintf("Found %d joint pathways; the manuscript reports %d.", nrow(hit), N_EXPECTED))

concord <- ifelse(sign(hit$Z_gene) != sign(hit$Z_met), "Discordant",
                  ifelse(hit$Z_gene > 0, "Concordant, increased", "Concordant, decreased"))

tab <- data.frame(
  Rank         = seq_len(nrow(hit)),
  `KEGG ID`    = hit$KEGG_id,
  Pathway      = hit$KEGG_title,
  `Z gene`     = round(hit$Z_gene, 2),
  `FDR gene`   = signif(hit$FDR_gene, 2),
  `Z met`      = round(hit$Z_met, 2),
  `P met`      = signif(hit$P_met, 2),
  `FDR (Pmax)` = signif(hit$FDR_Pmax, 2),
  Concordance  = concord,
  check.names = FALSE, stringsAsFactors = FALSE)
print(tab, row.names = FALSE)

n_conc <- sum(grepl("^Concordant", tab$Concordance))
n_up   <- sum(tab$Concordance == "Concordant, increased")
n_dn   <- sum(tab$Concordance == "Concordant, decreased")
n_disc <- sum(tab$Concordance == "Discordant")
n_pmax <- sum(hit$FDR_Pmax < PMAX_FDR_CUT)
miss_pmax <- tab$Pathway[hit$FDR_Pmax >= PMAX_FDR_CUT]
pmax_is_pmet <- all(abs(hit$Pmax - hit$P_met) < 1e-12)

utils::write.csv(tab, file.path(OUT_DIR, "Supplementary_Table_3.csv"), row.names = FALSE)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Supp table")
  openxlsx::writeData(wb, 1, sprintf("%s | Pathways jointly perturbed in the gene and metabolite arms under alcohol on the PNPLA3 I148M background (n = %d)",
                                     TABLE_NO, nrow(tab)), startRow = 1)
  openxlsx::writeData(wb, 1, tab, startRow = 3,
                      headerStyle = openxlsx::createStyle(textDecoration = "bold", border = "Bottom"))
  openxlsx::addStyle(wb, 1, openxlsx::createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  openxlsx::setColWidths(wb, 1, cols = 1:9, widths = c(6, 11, 46, 9, 10, 9, 9, 11, 22))
  openxlsx::saveWorkbook(wb, file.path(OUT_DIR, "Supplementary_Table_3.xlsx"), overwrite = TRUE)
}

## ---- Word table (Tables 1-3 share this style) ------------------------------
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  library(flextable); library(officer)
  FONT <- "Times New Roman"
  fmt_z <- function(x) ifelse(x > 0, sprintf("+%.2f", x), sub("-", "−", sprintf("%.2f", x)))
  sci_parts <- function(x) {
    e <- floor(log10(x)); m <- signif(x / 10^e, 2)
    list(mant = sprintf("%.1f × 10", m), exp = sub("-", "−", as.character(e)))
  }
  pq_par <- function(x) {
    if (x >= 0.001) return(as_paragraph(as_chunk(formatC(signif(x, 2), format = "fg", digits = 2, flag = "#"))))
    s <- sci_parts(x); as_paragraph(as_chunk(s$mant), as_sup(s$exp))
  }
  w <- data.frame(Rank = tab$Rank, KEGG = tab$`KEGG ID`, Pathway = tab$Pathway,
                  zg = fmt_z(tab$`Z gene`), fg = "", zm = fmt_z(tab$`Z met`), pm = "", fp = "",
                  Conc = sub("Concordant, increased", "Concordant ↑",
                             sub("Concordant, decreased", "Concordant ↓", tab$Concordance)),
                  stringsAsFactors = FALSE)
  ft <- flextable(w)
  ft <- set_header_labels(ft, Rank = "Rank", KEGG = "KEGG ID", Pathway = "Pathway",
                          zg = "", fg = "", zm = "", pm = "", fp = "", Conc = "Concordance")
  ft <- compose(ft, part = "header", j = "zg", value = as_paragraph(as_i("Z")))
  ft <- compose(ft, part = "header", j = "fg", value = as_paragraph("FDR"))
  ft <- compose(ft, part = "header", j = "zm", value = as_paragraph(as_i("Z")))
  ft <- compose(ft, part = "header", j = "pm", value = as_paragraph(as_i("P")))
  ft <- compose(ft, part = "header", j = "fp", value = as_paragraph("FDR (", as_i("P"), as_sub("max"), ")"))
  ft <- add_header_row(ft, colwidths = c(1, 1, 1, 2, 2, 1, 1),
                       values = c("", "", "", "Gene arm", "Metabolite arm", "", ""))
  for (j in c(1, 2, 3, 8, 9)) ft <- merge_at(ft, i = 1:2, j = j, part = "header")
  ft <- compose(ft, part = "header", i = 1, j = 1, value = as_paragraph("Rank"))
  ft <- compose(ft, part = "header", i = 1, j = 2, value = as_paragraph("KEGG ID"))
  ft <- compose(ft, part = "header", i = 1, j = 3, value = as_paragraph("Pathway"))
  ft <- compose(ft, part = "header", i = 1, j = 8, value = as_paragraph("Joint FDR (", as_i("P"), as_sub("max"), ")"))
  ft <- compose(ft, part = "header", i = 1, j = 9, value = as_paragraph("Concordance"))
  for (i in seq_len(nrow(w))) {
    ft <- compose(ft, i = i, j = "fg", value = pq_par(tab$`FDR gene`[i]))
    ft <- compose(ft, i = i, j = "pm", value = pq_par(tab$`P met`[i]))
    ft <- compose(ft, i = i, j = "fp", value = pq_par(tab$`FDR (Pmax)`[i]))
  }
  ## colour Z by direction, as in the earlier table: red increase, blue decrease
  for (k in c("zg", "zm")) {
    zz <- if (k == "zg") tab$`Z gene` else tab$`Z met`
    ft <- color(ft, i = which(zz > 0), j = k, color = "#C62828")
    ft <- color(ft, i = which(zz < 0), j = k, color = "#1565C0")
  }
  ## shade the row that fails the joint Pmax test
  shade <- which(hit$FDR_Pmax >= PMAX_FDR_CUT)
  ft <- font(ft, fontname = FONT, part = "all")
  ft <- fontsize(ft, size = 8.5, part = "all")
  ft <- bold(ft, part = "header")
  ft <- bg(ft, bg = "#E8EAED", part = "header")
  if (length(shade)) ft <- bg(ft, i = shade, bg = "#F5F6F7", part = "body")
  ft <- align(ft, j = c(1, 4:8), align = "center", part = "all")
  ft <- align(ft, j = c(2, 3, 9), align = "left", part = "all")
  ft <- valign(ft, valign = "center", part = "header")
  ft <- border_remove(ft)
  rule <- fp_border(color = "#999999", width = 1)
  hair <- fp_border(color = "#DDDDDD", width = 0.25)
  ft <- hline_top(ft, border = rule, part = "header")
  ft <- hline(ft, i = 1, j = 4:7, border = hair, part = "header")
  ft <- hline_bottom(ft, border = hair, part = "header")
  ft <- hline(ft, border = hair, part = "body")
  ft <- hline_bottom(ft, border = rule, part = "body")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, padding.left = 3.5,
                padding.right = 3.5, part = "all")
  ft <- width(ft, j = 1:9, width = c(0.45, 0.8, 2.9, 0.6, 0.9, 0.6, 0.7, 0.95, 1.2))
  ft <- set_table_properties(ft, layout = "fixed")

  tp <- function(size, ...) fp_text(font.family = FONT, font.size = size, ...)
  ch <- function(txt, size, ...) ftext(txt, tp(size, ...))
  pp <- function(...) fpar(..., fp_p = fp_par(padding.bottom = 7))
  title <- pp(ch(TABLE_NO, 12, bold = TRUE),
              ch(" | Pathways jointly perturbed in the gene and metabolite arms under alcohol on the PNPLA3", 12),
              ch("I148M", 12, vertical.align = "superscript"), ch(" background.", 12))
  intro <- pp(ch("A KEGG pathway qualifies when the gene arm reaches Benjamini–Hochberg FDR < ", 9.5),
              ch(sprintf("%.2f and the metabolite arm reaches ", GENE_FDR_CUT), 9.5),
              ch("P", 9.5, italic = TRUE),
              ch(sprintf(paste0(" < %.2f, among the %d pathways measured in both arms; this is the ",
                                "rule stated in Fig. 4a. Rows are ordered by metabolite "),
                         MET_P_CUT, nrow(both_arms)), 9.5),
              ch("P", 9.5, italic = TRUE), ch(".", 9.5))
  note <- pp(ch("Z", 9, italic = TRUE),
             ch(" gene, whole-liver cell-type-weighted Stouffer ", 9), ch("z", 9, italic = TRUE),
             ch(" across all 17 cell types; ", 9), ch("Z", 9, italic = TRUE),
             ch(" met, metabolite ssGSEA and limma ", 9), ch("z", 9, italic = TRUE),
             ch(" (whole liver, ", 9), ch("n", 9, italic = TRUE),
             ch(" = 9 mice); red, increased; blue, decreased. ", 9),
             ch("P", 9, italic = TRUE), ch("max", 9, vertical.align = "subscript"),
             ch(" is the larger of the two arms' ", 9), ch("P", 9, italic = TRUE),
             ch(sprintf(paste0(" values, and joint FDR is Benjamini–Hochberg on it. %d of the %d ",
                               "also pass joint FDR < %.2f; the shaded row (%s) is the one that does not. %d pathways are ",
                               "concordant (%d increased, %d decreased) and %d discordant."),
                        n_pmax, nrow(tab), PMAX_FDR_CUT, paste(miss_pmax, collapse = ", "),
                        n_conc, n_up, n_dn, n_disc), 9),
             if (pmax_is_pmet) ch(paste0(" The metabolite arm carries the larger P value on ",
                                         "every row, so it sets the order."), 9) else ch("", 9))
  src <- pp(ch("Source data: ", 9),
            ftext("modalities_pmax_gene_plus_metDec12.csv", fp_text(font.family = "Consolas", font.size = 9)),
            ch(" (integration output for alcohol on the variant background).", 9))
  landscape <- prop_section(page_size = page_size(orient = "landscape", width = 8.5, height = 11),
                            page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75))
  doc <- read_docx()
  doc <- body_add_fpar(doc, title)
  doc <- body_add_fpar(doc, intro)
  doc <- body_add_flextable(doc, ft, align = "left")
  doc <- body_add_fpar(doc, fpar(ftext("")))
  doc <- body_add_fpar(doc, note)
  doc <- body_add_fpar(doc, src)
  doc <- body_set_default_section(doc, landscape)
  print(doc, target = file.path(OUT_DIR, "Supplementary_Table_3.docx"))
}
