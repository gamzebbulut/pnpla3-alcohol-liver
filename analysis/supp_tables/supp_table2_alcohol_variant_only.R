# =============================================================================
# supp_table2_alcohol_variant_only.R
#
# Supplementary Table 2: the 98 metabolites altered by alcohol on the
# PNPLA3-I148M background but NOT on wild type -- the other side of the Venn
# whose overlap (9) is Supplementary Table 1. Requested by Katia.
#
# Manuscript sentence this table backs (Results, para 1):
#   "Only 9 metabolites responded to alcohol on both backgrounds ... whereas 98
#    responded only on the mutant background."
#
# Same two comparisons, same significance rule as Supplementary Table 1:
#   FC_Q3_metabolites.csv  alcohol on the PNPLA3-I148M background (HFSW + variant + alcohol vs HFSW + variant)
#   FC_Q4_metabolites.csv  alcohol on the wild-type background    (HFSW + alcohol vs HFSW)
#   significant = two-sided Welch P < 0.05 AND Storey q < 0.10, within each comparison.
#
# "Only on the variant background" means it missed P < 0.05 / q < 0.10 on wild
# type. It does NOT mean wild type had no effect: the wild-type log2 FC, P and q
# are shown for every row, and the footnote counts how many changed in the same
# direction on wild type and how many reached nominal P < 0.05 there.
#
# The 5 unannotated features (FFA42, ChoE_Col_02, ChoE_Col_03, MAPE09, MAPE36)
# all fall in this set. They are counted in the 98 in the text. INCLUDE_UNANNOTATED
# decides whether they get rows here (TRUE -> 98 rows) or a footnote (FALSE -> 93).
#
# Outputs (written to OUT_DIR):
#   Supplementary_Table_2.csv / .xlsx / .docx
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

TABLE_NO            <- "Supplementary Table 2"   # 1 = the 9 shared; 3 = joint pathways
SHARED_TABLE        <- "Supplementary Table 1"
INCLUDE_UNANNOTATED <- TRUE                      # TRUE: 98 rows; FALSE: 93 rows + footnote

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
only <- setdiff(sig3, sig4)

cat(sprintf("variant bg %d | wild type %d | both %d | variant bg only %d\n",
            length(sig3), length(sig4), length(both), length(only)))
# Guard: the manuscript text says 107 / 32 / 9 / 98. Stop loudly if inputs drift.
if (!identical(c(length(sig3), length(sig4), length(both), length(only)),
               c(107L, 32L, 9L, 98L)))
  stop("Counts no longer match the manuscript text (107 / 32 / 9 / 98).")

## ---- display names (same rules as the 9-metabolite table) -------------------
display_name <- function(d) {
  nm <- d[["INDIVIDUAL.NOTATION"]]
  nm[is.na(nm) | nm == ""] <- d[["COMMON.NAME"]][is.na(nm) | nm == ""]
  unann <- is.na(nm) | nm == ""
  nm[unann] <- paste0("Unannotated (", d$variable[unann], ")")
  # Lyso species: PC(0:0/20:5) or PC(20:5/0:0) -> LPC(20:5), as in the text/Fig 2d
  lyso <- grepl("^Monoacyl|^[12]-Monoacyl", d[["SUBCLASS.A"]])
  nm[lyso] <- sub("^P([A-Z])\\((?:0:0/)?([0-9]+:[0-9]+)(?:/0:0)?\\)$",
                  "LP\\1(\\2)", nm[lyso], perl = TRUE)
  nm
}
is_unannotated <- function(d) {
  nm <- d[["INDIVIDUAL.NOTATION"]]
  nm[is.na(nm) | nm == ""] <- d[["COMMON.NAME"]][is.na(nm) | nm == ""]
  is.na(nm) | nm == ""
}
class_label <- function(x) {
  x <- ifelse(x %in% c("Amino acids", "Amino acid", "Amino Acids and Derivatives",
                       "Amino acid derivatives"), "Amino acids and derivatives", x)
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))   # sentence case
}

## ---- build the table --------------------------------------------------------
keep <- c("variable", "log2_FC", "p", "qvalue")
a <- q3[q3$variable %in% only, ]
b <- q4[q4$variable %in% only, keep]
names(b)[-1] <- paste0(names(b)[-1], "_wt")
t98 <- merge(a, b, by = "variable")
t98$unann <- is_unannotated(t98)

n_unann <- sum(t98$unann)
if (!INCLUDE_UNANNOTATED) t98 <- t98[!t98$unann, ]

tab <- data.frame(
  Metabolite        = display_name(t98),
  Class             = class_label(t98$CLASS),
  HMDB              = ifelse(is.na(t98$HMDB), "", t98$HMDB),
  q3_fc = round(t98$log2_FC, 2),    q3_p = signif(t98$p, 2),    q3_q = signif(t98$qvalue, 2),
  q4_fc = round(t98$log2_FC_wt, 2), q4_p = signif(t98$p_wt, 2), q4_q = signif(t98$qvalue_wt, 2),
  Direction         = ifelse(t98$log2_FC > 0, "Increased", "Decreased"),
  `Wild type`       = ifelse(sign(t98$log2_FC_wt) == sign(t98$log2_FC), "Same direction", "Opposite"),
  `Feature ID`      = t98$variable,
  check.names = FALSE, stringsAsFactors = FALSE)

# Amino acids first, then the other classes alphabetically; within class by
# effect in the lead comparison (largest increase first); unannotated last in class
cls_order <- c("Amino acids and derivatives", sort(setdiff(unique(tab$Class),
                                                           "Amino acids and derivatives")))
tab <- tab[order(match(tab$Class, cls_order), grepl("^Unannotated", tab$Metabolite),
                 -tab$q3_fc), ]
rownames(tab) <- NULL

## footnote numbers, computed from the rows actually shown
n_rows    <- nrow(tab)
n_up      <- sum(tab$Direction == "Increased")
n_dn      <- sum(tab$Direction == "Decreased")
n_same_wt <- sum(tab$`Wild type` == "Same direction")
n_nom_wt  <- sum(!is.na(tab$q4_p) & tab$q4_p < P_CUT)
cat(sprintf("rows %d | up %d | down %d | same direction on WT %d | nominal P<%.2f on WT %d\n",
            n_rows, n_up, n_dn, n_same_wt, P_CUT, n_nom_wt))

hdr <- c(paste(LAB_Q3, c("log2 FC", "P", "q")), paste(LAB_Q4, c("log2 FC", "P", "q")))
csv_out <- tab; names(csv_out)[4:9] <- hdr
write.csv(csv_out, file.path(OUT_DIR, "Supplementary_Table_2.csv"),
          row.names = FALSE, na = "")

## ---- formatted xlsx ---------------------------------------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, "Supp table")
  writeData(wb, 1, sprintf(paste0("%s | Metabolites altered by alcohol on the PNPLA3 I148M ",
                                  "background but not on wild type (n = %d)"), TABLE_NO, n_rows),
            startRow = 1)
  writeData(wb, 1, t(c(rep("", 3), LAB_Q3, "", "", LAB_Q4, "", "", "", "", "")),
            startRow = 3, colNames = FALSE)
  mergeCells(wb, 1, cols = 4:6, rows = 3); mergeCells(wb, 1, cols = 7:9, rows = 3)
  body <- tab
  names(body) <- c(names(tab)[1:3], "log2 FC", "P", "q", "log2 FC", "P", "q", names(tab)[10:12])
  writeData(wb, 1, body, startRow = 4, headerStyle =
              createStyle(textDecoration = "bold", border = "Bottom", halign = "center"))
  addStyle(wb, 1, createStyle(textDecoration = "bold", halign = "center", border = "Bottom"),
           rows = 3, cols = 4:9, gridExpand = TRUE)
  addStyle(wb, 1, createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  setColWidths(wb, 1, cols = 1:12, widths = c(22, 26, 11, 9, 9, 9, 9, 9, 9, 11, 15, 13))
  freezePane(wb, 1, firstActiveRow = 5, firstActiveCol = 2)
  saveWorkbook(wb, file.path(OUT_DIR, "Supplementary_Table_2.xlsx"), overwrite = TRUE)
}

## ---- Word table (Tables 1-3 share this style) ------------------------------
## Times New Roman; grey header band (#E8EAED); grey top/bottom rules, hairline
## row rules; italic P / q / n; landscape Letter. Class is a spanning row, not a
## column, so 98 rows stay readable.
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  library(flextable); library(officer)

  FONT <- "Times New Roman"
  fmt_fc <- function(x) ifelse(x > 0, sprintf("+%.2f", x), sub("-", "−", sprintf("%.2f", x)))
  sci_parts <- function(x) {
    e <- floor(log10(x)); m <- signif(x / 10^e, 2)
    list(mant = sprintf("%.1f × 10", m), exp = sub("-", "−", as.character(e)))
  }
  pq_par <- function(x) {
    if (is.na(x)) return(as_paragraph(as_chunk("")))
    if (x >= 0.001) return(as_paragraph(as_chunk(formatC(signif(x, 2), format = "fg", digits = 2, flag = "#"))))
    s <- sci_parts(x); as_paragraph(as_chunk(s$mant), as_sup(s$exp))
  }

  ## interleave one class-heading row before each class block
  rows <- list(); cls_row <- logical(0); num <- list()
  for (cl in unique(tab$Class)) {
    blk <- tab[tab$Class == cl, ]
    rows[[length(rows) + 1]] <- data.frame(Metabolite = sprintf("%s (%d)", cl, nrow(blk)),
      fc3 = "", p3 = "", q3 = "", fc4 = "", p4 = "", q4 = "", Direction = "", WT = "",
      stringsAsFactors = FALSE)
    cls_row <- c(cls_row, TRUE); num[[length(num) + 1]] <- rep(NA_real_, 6)
    for (i in seq_len(nrow(blk))) {
      rows[[length(rows) + 1]] <- data.frame(Metabolite = blk$Metabolite[i],
        fc3 = fmt_fc(blk$q3_fc[i]), p3 = "", q3 = "",
        fc4 = fmt_fc(blk$q4_fc[i]), p4 = "", q4 = "",
        Direction = ifelse(blk$Direction[i] == "Increased", "Increased ↑", "Decreased ↓"),
        WT = ifelse(blk$`Wild type`[i] == "Same direction", "Same", "Opposite"),
        stringsAsFactors = FALSE)
      cls_row <- c(cls_row, FALSE)
      num[[length(num) + 1]] <- unlist(blk[i, c("q3_p", "q3_q", "q4_p", "q4_q")])
    }
  }
  w <- do.call(rbind, rows)

  ft <- flextable(w)
  ft <- set_header_labels(ft, Metabolite = "Metabolite", fc3 = "", p3 = "", q3 = "",
                          fc4 = "", p4 = "", q4 = "", Direction = "Direction",
                          WT = "Wild type")
  for (k in c("fc3", "fc4")) ft <- compose(ft, part = "header", j = k,
                                           value = as_paragraph("log", as_sub("2"), " FC"))
  for (k in c("p3", "p4"))   ft <- compose(ft, part = "header", j = k, value = as_paragraph(as_i("P")))
  for (k in c("q3", "q4"))   ft <- compose(ft, part = "header", j = k, value = as_paragraph(as_i("q")))
  ft <- add_header_row(ft, colwidths = c(1, 3, 3, 1, 1), values = c("", "", "", "", ""))
  ft <- compose(ft, part = "header", i = 1, j = 2,
                value = as_paragraph("Alcohol (PNPLA3", as_sup("I148M"), " bg)"))
  ft <- compose(ft, part = "header", i = 1, j = 5, value = as_paragraph("Alcohol alone"))
  ft <- compose(ft, part = "header", i = 1, j = 1, value = as_paragraph("Metabolite"))
  ft <- compose(ft, part = "header", i = 1, j = 8, value = as_paragraph("Direction"))
  ft <- compose(ft, part = "header", i = 1, j = 9, value = as_paragraph("Wild type"))
  for (j in c(1, 8, 9)) ft <- merge_at(ft, i = 1:2, j = j, part = "header")

  ## P / q with superscript exponents where needed
  pq_cols <- c(p3 = 1, q3 = 2, p4 = 3, q4 = 4)
  for (i in which(!cls_row)) for (k in names(pq_cols))
    ft <- compose(ft, i = i, j = k, value = pq_par(num[[i]][pq_cols[[k]]]))

  ## class rows span the table
  for (i in which(cls_row)) ft <- merge_at(ft, i = i, j = 1:9, part = "body")

  ft <- font(ft, fontname = FONT, part = "all")
  ft <- fontsize(ft, size = 8.5, part = "all")
  ft <- bold(ft, part = "header")
  ft <- bg(ft, bg = "#E8EAED", part = "header")
  ft <- bold(ft, i = which(cls_row), part = "body")
  ft <- italic(ft, i = which(cls_row), part = "body")
  ft <- bg(ft, i = which(cls_row), bg = "#F5F6F7", part = "body")
  ft <- align(ft, j = 2:7, align = "center", part = "all")
  ft <- align(ft, j = c(1, 8, 9), align = "left", part = "all")
  ft <- valign(ft, valign = "center", part = "header")
  ft <- border_remove(ft)
  rule <- fp_border(color = "#999999", width = 1)
  hair <- fp_border(color = "#DDDDDD", width = 0.25)
  ft <- hline_top(ft, border = rule, part = "header")
  ft <- hline(ft, i = 1, j = 2:7, border = hair, part = "header")
  ft <- hline_bottom(ft, border = hair, part = "header")
  ft <- hline(ft, border = hair, part = "body")
  ft <- hline_bottom(ft, border = rule, part = "body")
  ft <- padding(ft, padding.top = 1.5, padding.bottom = 1.5, padding.left = 3.5,
                padding.right = 3.5, part = "all")
  ft <- width(ft, j = 1:9, width = c(1.9, 0.8, 0.95, 0.7, 0.8, 0.95, 0.7, 1.0, 0.8))
  ft <- set_table_properties(ft, layout = "fixed")

  tp <- function(size, ...) fp_text(font.family = FONT, font.size = size, ...)
  ch <- function(txt, size, ...) ftext(txt, tp(size, ...))
  pp <- function(...) fpar(..., fp_p = fp_par(padding.bottom = 7))
  words <- function(n) if (n <= 10) c("Zero","One","Two","Three","Four","Five","Six",
                                      "Seven","Eight","Nine","Ten")[n + 1] else as.character(n)

  title <- pp(ch(TABLE_NO, 12, bold = TRUE),
              ch(" | Metabolites altered by alcohol on the PNPLA3", 12),
              ch("I148M", 12, vertical.align = "superscript"),
              ch(" background but not on wild type.", 12))
  intro <- pp(ch("A metabolite qualifies when it reaches ", 9.5), ch("P", 9.5, italic = TRUE),
              ch(sprintf(" < %.2f and ", P_CUT), 9.5), ch("q", 9.5, italic = TRUE),
              ch(sprintf(paste0(" < %.2f with alcohol on the variant background but not on wild type: ",
                                "%d of the %d metabolites altered on the variant background, the ",
                                "remaining %d being shared with wild type (%s). ",
                                "Rows are grouped by class and ordered by log"),
                         Q_CUT, length(only), length(sig3), length(both), SHARED_TABLE), 9.5),
              ch("2", 9.5, vertical.align = "subscript"),
              ch(" fold-change on the variant background.", 9.5))
  unann_txt <- if (INCLUDE_UNANNOTATED)
    sprintf(" %s rows are unannotated features, listed by feature ID.", words(n_unann))
  else sprintf(paste0(" %s unannotated features are counted in the %d but not listed ",
                      "(%s)."), words(n_unann), length(only),
               paste(sort(q3$variable[q3$variable %in% only & is_unannotated(q3)]),
                     collapse = ", "))
  note <- pp(ch("Log", 9), ch("2", 9, vertical.align = "subscript"),
             ch(" FC is the difference in mean log", 9),
             ch("2", 9, vertical.align = "subscript"),
             ch(" abundance, alcohol minus no alcohol (positive = higher with alcohol). ", 9),
             ch("P", 9, italic = TRUE), ch(" is from a two-sided Welch's ", 9),
             ch("t", 9, italic = TRUE), ch("-test and ", 9), ch("q", 9, italic = TRUE),
             ch(" is the Storey ", 9), ch("q", 9, italic = TRUE),
             ch(" value within each comparison; ", 9), ch("n", 9, italic = TRUE),
             ch(" = 5 mice per group, except variant plus alcohol (", 9),
             ch("n", 9, italic = TRUE),
             ch(sprintf(paste0(" = 4). On the variant background, %d increased and %d decreased. ",
                               "Wild type gives the direction of the same metabolite under alcohol ",
                               "on wild type: %d of %d moved the same way there, and %d reached "),
                        n_up, n_dn, n_same_wt, n_rows, n_nom_wt), 9),
             ch("P", 9, italic = TRUE), ch(sprintf(" < %.2f but not ", P_CUT), 9),
             ch("q", 9, italic = TRUE),
             ch(sprintf(" < %.2f there.%s", Q_CUT, unann_txt), 9))
  src <- pp(ch("Source data: ", 9),
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
  print(doc, target = file.path(OUT_DIR, "Supplementary_Table_2.docx"))
}
