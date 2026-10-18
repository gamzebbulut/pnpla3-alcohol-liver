# =============================================================================
# supp_table1_alcohol_both_backgrounds.R
#
# Supplementary Table 1: the 9 metabolites altered by alcohol on BOTH the
# PNPLA3-I148M and the wild-type background (Sanyal: "Only 9 ... Which 9?").
# Supplementary Table 2 (supp_table2_alcohol_variant_only.R) lists the other 98.
#
# Intersects the two alcohol comparisons:
#   FC_Q3_metabolites.csv  alcohol on the PNPLA3-I148M background (HFSW + variant + alcohol vs HFSW + variant)
#   FC_Q4_metabolites.csv  alcohol on the wild-type background    (HFSW + alcohol vs HFSW)
#
# Significance rule = the one the manuscript uses for the 107 / 32 counts:
#   two-sided Welch P < 0.05 AND Storey q < 0.10, within each comparison.
#   (Raw P < 0.05 alone gives 23 shared, which is NOT what the text reports.)
#
# Outputs (written to OUT_DIR):
#   Supplementary_Table_1.csv    the 9-metabolite table, plain
#   Supplementary_Table_1.xlsx   formatted table + Venn counts + full per-arm lists
#   Supplementary_Table_1.docx   Word table for submission
#
# Paths follow the repo convention: run setup.R first from a clone, or leave the
# defaults to run from the Desktop project folder.
# =============================================================================

## ---- tuning block -----------------------------------------------------------
ROOT    <- Sys.getenv("LIPID_ROOT", path.expand("~/Desktop/Lipidomics project with Dr Smirnova"))
FIGROOT <- Sys.getenv("FIGROOT",    file.path(ROOT, "writing skeleton/figure making"))
FC_DIR  <- Sys.getenv("FC_DIR",     file.path(FIGROOT, "Fig 1"))
OUT_DIR <- Sys.getenv("SUPPTAB_OUT", file.path(FIGROOT, "supp_tables/out"))
Q3_FILE <- file.path(FC_DIR, "FC_Q3_metabolites.csv")
Q4_FILE <- file.path(FC_DIR, "FC_Q4_metabolites.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

P_CUT <- 0.05
Q_CUT <- 0.10
TABLE_NO <- "Supplementary Table 1"

## Approved contrast labels (Supp Fig 1)
LAB_Q3 <- "Alcohol (PNPLA3 I148M bg)"
LAB_Q4 <- "Alcohol alone"

## UTF-8 so the minus, times and arrow glyphs survive into the Word table
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))
if (!l10n_info()$`UTF-8`) invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))

## ---- read -------------------------------------------------------------------
q3 <- read.csv(Q3_FILE, check.names = FALSE, stringsAsFactors = FALSE)
q4 <- read.csv(Q4_FILE, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(!anyDuplicated(q3$variable), !anyDuplicated(q4$variable),
          setequal(q3$variable, q4$variable))

is_sig <- function(d) !is.na(d$p) & !is.na(d$qvalue) & d$p < P_CUT & d$qvalue < Q_CUT
sig3 <- q3$variable[is_sig(q3)]
sig4 <- q4$variable[is_sig(q4)]
both <- intersect(sig3, sig4)

venn <- data.frame(
  Set   = c(paste(LAB_Q3, "(total)"), paste(LAB_Q4, "(total)"),
            "Both backgrounds", paste(LAB_Q3, "only"), paste(LAB_Q4, "only")),
  Count = c(length(sig3), length(sig4), length(both),
            length(setdiff(sig3, sig4)), length(setdiff(sig4, sig3))))
print(venn, row.names = FALSE)

# Guard: the manuscript text says 107 / 32 / 9 / 98. Stop loudly if inputs drift.
if (!identical(venn$Count[1:4], c(107L, 32L, 9L, 98L)))
  stop("Counts no longer match the manuscript text (107 / 32 / 9 / 98).")

## ---- display names ----------------------------------------------------------
display_name <- function(d) {
  nm <- d[["INDIVIDUAL.NOTATION"]]
  nm[is.na(nm) | nm == ""] <- d[["COMMON.NAME"]][is.na(nm) | nm == ""]
  nm[is.na(nm) | nm == ""] <- d$variable[is.na(nm) | nm == ""]
  # Lyso species: PC(0:0/20:5) or PC(20:5/0:0) -> LPC(20:5), as in the text/Fig 2d
  lyso <- grepl("^Monoacyl|^[12]-Monoacyl", d[["SUBCLASS.A"]])
  nm[lyso] <- sub("^P([A-Z])\\((?:0:0/)?([0-9]+:[0-9]+)(?:/0:0)?\\)$",
                  "LP\\1(\\2)", nm[lyso], perl = TRUE)
  nm
}
class_label <- function(x) {
  x <- ifelse(x %in% c("Amino acids", "Amino Acids and Derivatives",
                       "Amino acid derivatives"), "Amino acids and derivatives", x)
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))   # sentence case
}

## ---- build the table --------------------------------------------------------
keep <- c("variable", "log2_FC", "p", "qvalue")
a <- q3[q3$variable %in% both, ]
b <- q4[q4$variable %in% both, keep]
names(b)[-1] <- paste0(names(b)[-1], "_wt")
t9 <- merge(a, b, by = "variable")

direction <- ifelse(sign(t9$log2_FC) == sign(t9$log2_FC_wt),
                    ifelse(t9$log2_FC > 0, "Increased on both", "Decreased on both"),
                    "Opposite")

tab <- data.frame(
  Metabolite        = display_name(t9),
  Class             = class_label(t9$CLASS),
  HMDB              = ifelse(is.na(t9$HMDB), "", t9$HMDB),
  q3_fc = round(t9$log2_FC, 2),    q3_p = signif(t9$p, 2),    q3_q = signif(t9$qvalue, 2),
  q4_fc = round(t9$log2_FC_wt, 2), q4_p = signif(t9$p_wt, 2), q4_q = signif(t9$qvalue_wt, 2),
  Direction         = direction,
  `Feature ID`      = t9$variable,
  check.names = FALSE, stringsAsFactors = FALSE)

# Amino acids first, then lipids; within class by effect in the lead comparison
cls_order <- c("Amino acids and derivatives", sort(setdiff(unique(tab$Class),
                                                           "Amino acids and derivatives")))
tab <- tab[order(match(tab$Class, cls_order), -tab$q3_fc), ]

hdr <- c(paste(LAB_Q3, c("log2 FC", "P", "q")), paste(LAB_Q4, c("log2 FC", "P", "q")))
names(tab)[4:9] <- hdr
print(tab, row.names = FALSE)

write.csv(tab, file.path(OUT_DIR, "Supplementary_Table_1.csv"),
          row.names = FALSE, na = "")

## ---- full per-arm lists (for the Venn / for Isabel) -------------------------
arm_list <- function(d, ids, set_lab) {
  s <- d[d$variable %in% ids, ]
  data.frame(Metabolite = display_name(s), Class = class_label(s$CLASS),
             `log2 FC` = round(s$log2_FC, 2), P = signif(s$p, 2), q = signif(s$qvalue, 2),
             Set = set_lab, `Feature ID` = s$variable,
             check.names = FALSE)[order(class_label(s$CLASS), -s$log2_FC), ]
}
l3 <- arm_list(q3, sig3, NA); l3$Set <- ifelse(l3$`Feature ID` %in% both, "Both", paste(LAB_Q3, "only"))
l4 <- arm_list(q4, sig4, NA); l4$Set <- ifelse(l4$`Feature ID` %in% both, "Both", paste(LAB_Q4, "only"))

## ---- formatted xlsx ---------------------------------------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, "Supp table")
  title <- sprintf(paste0("Supplementary Table 1 | Metabolites altered by alcohol on both the ",
                          "PNPLA3 I148M and wild-type backgrounds (n = %d)"), nrow(tab))
  note <- paste0("Significance in each comparison: two-sided Welch's t-test P < ", P_CUT,
                 " and Storey q < ", Q_CUT, ". log2 FC is the difference in mean log2 abundance, ",
                 "alcohol minus no alcohol (positive = higher with alcohol). ",
                 "n = 5 mice per group, except variant plus alcohol (n = 4). ",
                 "Of ", length(sig3), " metabolites altered on the variant background and ",
                 length(sig4), " on wild type, ", length(both), " were shared; all changed in the ",
                 "same direction on both backgrounds.")
  writeData(wb, 1, title, startRow = 1)
  # two-level header: comparison spans over log2 FC / P / q
  writeData(wb, 1, t(c(rep("", 3), LAB_Q3, "", "", LAB_Q4, "", "", "", "")),
            startRow = 3, colNames = FALSE)
  mergeCells(wb, 1, cols = 4:6, rows = 3); mergeCells(wb, 1, cols = 7:9, rows = 3)
  sub_hdr <- c(names(tab)[1:3], "log2 FC", "P", "q", "log2 FC", "P", "q", names(tab)[10:11])
  body <- tab; names(body) <- sub_hdr
  writeData(wb, 1, body, startRow = 4, headerStyle =
              createStyle(textDecoration = "bold", border = "Bottom", halign = "center"))
  addStyle(wb, 1, createStyle(textDecoration = "bold", halign = "center", border = "Bottom"),
           rows = 3, cols = 4:9, gridExpand = TRUE)
  addStyle(wb, 1, createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  addStyle(wb, 1, createStyle(halign = "center"), rows = 5:(4 + nrow(tab)), cols = 4:9,
           gridExpand = TRUE, stack = TRUE)
  writeData(wb, 1, note, startRow = 6 + nrow(tab))
  addStyle(wb, 1, createStyle(wrapText = TRUE, valign = "top", fontSize = 9),
           rows = 6 + nrow(tab), cols = 1)
  mergeCells(wb, 1, cols = 1:11, rows = 6 + nrow(tab)); setRowHeights(wb, 1, 6 + nrow(tab), 48)
  setColWidths(wb, 1, cols = 1:11, widths = c(16, 26, 11, 9, 9, 9, 9, 9, 9, 18, 11))
  freezePane(wb, 1, firstActiveRow = 5, firstActiveCol = 2)

  addWorksheet(wb, "Venn counts");               writeData(wb, "Venn counts", venn)
  addWorksheet(wb, "Variant bg - all 107");      writeData(wb, "Variant bg - all 107", l3)
  addWorksheet(wb, "Wild type - all 32");        writeData(wb, "Wild type - all 32", l4)
  for (s in 2:4) setColWidths(wb, s, cols = 1:7, widths = "auto")
  saveWorkbook(wb, file.path(OUT_DIR, "Supplementary_Table_1.xlsx"), overwrite = TRUE)
}

## ---- formatted Word table (Tables 1-3 share this style) ---------------------
## Times New Roman; grey header band (#E8EAED); grey top/bottom rules, hairline
## row rules; italic P / q / n; landscape Letter.

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  library(flextable); library(officer)

  FONT <- "Times New Roman"
  fmt_fc <- function(x) ifelse(x > 0, sprintf("+%.2f", x), sub("-", "−", sprintf("%.2f", x)))
  # P / q: 2 significant figures; below 0.001 as "a × 10^b" (superscript built later)
  sci_parts <- function(x) {
    e <- floor(log10(x)); m <- signif(x / 10^e, 2)
    list(mant = sprintf("%.1f × 10", m), exp = sub("-", "−", as.character(e)))
  }
  pq_par <- function(x) {
    if (x >= 0.001) return(as_paragraph(as_chunk(formatC(signif(x, 2), format = "fg", digits = 2, flag = "#"))))
    s <- sci_parts(x); as_paragraph(as_chunk(s$mant), as_sup(s$exp))
  }

  w <- data.frame(Metabolite = tab$Metabolite, Class = tab$Class,
                  fc3 = tab[[4]], p3 = tab[[5]], q3 = tab[[6]],
                  fc4 = tab[[7]], p4 = tab[[8]], q4 = tab[[9]],
                  Direction = ifelse(tab$Direction == "Increased on both", "Increased ↑",
                              ifelse(tab$Direction == "Decreased on both", "Decreased ↓", "Opposite")),
                  stringsAsFactors = FALSE)
  w$fc3 <- fmt_fc(w$fc3); w$fc4 <- fmt_fc(w$fc4)
  num_raw <- tab[, 4:9]   # keep numeric P / q for the superscript formatting

  ft <- flextable(w)
  ft <- set_header_labels(ft, Metabolite = "Metabolite", Class = "Class",
                          fc3 = "", p3 = "", q3 = "", fc4 = "", p4 = "", q4 = "",
                          Direction = "Direction")
  # second-level header: log2 FC / P / q under each comparison
  for (k in c("fc3", "fc4")) ft <- compose(ft, part = "header", j = k,
                                           value = as_paragraph("log", as_sub("2"), " FC"))
  for (k in c("p3", "p4"))   ft <- compose(ft, part = "header", j = k, value = as_paragraph(as_i("P")))
  for (k in c("q3", "q4"))   ft <- compose(ft, part = "header", j = k, value = as_paragraph(as_i("q")))
  ft <- add_header_row(ft, colwidths = c(1, 1, 3, 3, 1), values = c("", "", "", "", ""))
  ft <- compose(ft, part = "header", i = 1, j = 3,
                value = as_paragraph("Alcohol (PNPLA3", as_sup("I148M"), " bg)"))
  ft <- compose(ft, part = "header", i = 1, j = 6, value = as_paragraph("Alcohol alone"))
  ft <- compose(ft, part = "header", i = 1, j = 1, value = as_paragraph("Metabolite"))
  ft <- compose(ft, part = "header", i = 1, j = 2, value = as_paragraph("Class"))
  ft <- compose(ft, part = "header", i = 1, j = 9, value = as_paragraph("Direction"))
  ft <- merge_at(ft, i = 1:2, j = 1, part = "header")
  ft <- merge_at(ft, i = 1:2, j = 2, part = "header")
  ft <- merge_at(ft, i = 1:2, j = 9, part = "header")

  # body P / q with superscript exponents where needed
  cols_pq <- c(p3 = 2, q3 = 3, p4 = 5, q4 = 6)
  for (k in names(cols_pq)) for (i in seq_len(nrow(w)))
    ft <- compose(ft, i = i, j = k, value = pq_par(num_raw[i, cols_pq[[k]]]))

  # styling
  ft <- font(ft, fontname = FONT, part = "all")
  ft <- fontsize(ft, size = 8.5, part = "all")
  ft <- bold(ft, part = "header")
  ft <- bg(ft, bg = "#E8EAED", part = "header")
  ft <- align(ft, j = 3:8, align = "center", part = "all")
  ft <- align(ft, j = c(1, 2, 9), align = "left", part = "all")
  ft <- valign(ft, valign = "center", part = "header")
  ft <- border_remove(ft)
  rule  <- fp_border(color = "#999999", width = 1)
  hair  <- fp_border(color = "#DDDDDD", width = 0.25)
  ft <- hline_top(ft, border = rule, part = "header")
  ft <- hline(ft, i = 1, j = 3:8, border = hair, part = "header")
  ft <- hline_bottom(ft, border = hair, part = "header")
  ft <- hline(ft, border = hair, part = "body")
  ft <- hline_bottom(ft, border = rule, part = "body")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, padding.left = 3.5,
                padding.right = 3.5, part = "all")
  ft <- width(ft, j = 1:9, width = c(1.1, 1.9, 0.8, 0.95, 0.7, 0.8, 0.95, 0.7, 1.1))

  # text pieces
  tp  <- function(size, ...) fp_text(font.family = FONT, font.size = size, ...)
  ch  <- function(txt, size, ...) ftext(txt, tp(size, ...))
  pp  <- function(...) fpar(..., fp_p = fp_par(padding.bottom = 7))
  title <- pp(ch(TABLE_NO, 12, bold = TRUE),
              ch(" | Metabolites altered by alcohol on both the wild-type and the PNPLA3", 12),
              ch("I148M", 12, vertical.align = "superscript"), ch(" backgrounds.", 12))
  intro <- pp(ch("A metabolite qualifies when it reaches ", 9.5), ch("P", 9.5, italic = TRUE),
              ch(sprintf(" < %.2f and ", P_CUT), 9.5), ch("q", 9.5, italic = TRUE),
              ch(sprintf(paste0(" < %.2f in both alcohol comparisons. %s of the %d metabolites ",
                                "altered on the variant background and %d altered on wild type ",
                                "qualify; the other %d altered on the variant background are in Supplementary Table 2. ",
                                "Rows are grouped by class and ordered by log"),
                         Q_CUT, c("Zero","One","Two","Three","Four","Five","Six","Seven","Eight",
                                  "Nine","Ten")[nrow(tab) + 1], length(sig3), length(sig4),
                         length(setdiff(sig3, sig4))), 9.5),
              ch("2", 9.5, vertical.align = "subscript"),
              ch(" fold-change on the variant background.", 9.5))
  n_up <- sum(w$Direction == "Increased ↑"); n_dn <- sum(w$Direction == "Decreased ↓")
  note <- pp(ch("Log", 9), ch("2", 9, vertical.align = "subscript"),
             ch(paste0(" FC is the difference in mean log", ""), 9),
             ch("2", 9, vertical.align = "subscript"),
             ch(" abundance, alcohol minus no alcohol (positive = higher with alcohol). ", 9),
             ch("P", 9, italic = TRUE), ch(" is from a two-sided Welch's ", 9),
             ch("t", 9, italic = TRUE), ch("-test and ", 9), ch("q", 9, italic = TRUE),
             ch(" is the Storey ", 9), ch("q", 9, italic = TRUE), ch(" value within each comparison; ", 9),
             ch("n", 9, italic = TRUE),
             ch(" = 5 mice per group, except variant plus alcohol (", 9),
             ch("n", 9, italic = TRUE),
             ch(sprintf(paste0(" = 4). All %d metabolites changed in the same direction on both ",
                               "backgrounds — %d increased and %d decreased."),
                        nrow(tab), n_up, n_dn), 9))
  src  <- pp(ch("Source data: ", 9),
             ftext("FC_Q3_metabolites.csv", fp_text(font.family = "Consolas", font.size = 9)),
             ch(" and ", 9),
             ftext("FC_Q4_metabolites.csv", fp_text(font.family = "Consolas", font.size = 9)),
             ch(".", 9))

  landscape <- prop_section(page_size = page_size(orient = "landscape", width = 8.5, height = 11),
                            page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75,
                                                    right = 0.75))
  doc <- read_docx()
  doc <- body_add_fpar(doc, title)
  doc <- body_add_fpar(doc, intro)
  doc <- body_add_flextable(doc, ft, align = "left")
  doc <- body_add_fpar(doc, fpar(ftext("")))
  doc <- body_add_fpar(doc, note)
  doc <- body_add_fpar(doc, src)
  doc <- body_set_default_section(doc, landscape)
  print(doc, target = file.path(OUT_DIR, "Supplementary_Table_1.docx"))
}
