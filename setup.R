## Point the analysis scripts at this repository.
##
## Every script resolves its inputs from LIPID_ROOT and FIGROOT, with all other
## paths derived from those two by default. Setting them here is what lets the
## scripts run unmodified from a fresh clone.
##
##   source("setup.R")   # then source any script under analysis/
##
## Works whether you source it by name from the repository root or by full path
## from anywhere: it locates itself rather than trusting the working directory.

.repo_root <- local({
  ## When sourced, R records the file's own path in the calling frame.
  f <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- of; break }
  }
  if (!is.null(f)) dirname(normalizePath(f, mustWork = TRUE))
  else normalizePath(".", mustWork = TRUE)   # sourced some other way; assume cwd
})

if (!file.exists(file.path(.repo_root, "run_all.R")))
  stop("setup.R does not appear to be at the repository root.\n",
       "  looked in: ", .repo_root, "\n",
       "  either setwd() to the repository and source(\"setup.R\"), or\n",
       "  source it by its full path.")

Sys.setenv(LIPID_ROOT = .repo_root)
Sys.setenv(FIGROOT    = file.path(.repo_root, "analysis"))

## Figure 2d-f read two sheets of the raw metabolomics workbook. OWL_CACHE lets
## them read those two sheets as CSVs instead. If that cache is present in this
## repository, use it -- that is what makes 2d-f reproduce without the workbook.
.owl <- file.path(.repo_root, "analysis/Fig 2/_owl_cache")
if (dir.exists(.owl)) Sys.setenv(OWL_CACHE = .owl)

## The snRNA-seq pseudobulk objects are not distributed with this repository
## (DATA.md). Scripts that need them stop with a clear message rather than
## half-run; the panels that display them read pre-computed CSVs instead.
if (!dir.exists(file.path(.repo_root, "Data/pseudobulk_snRNA-seq_PNPLA3_alc")))
  message("note: pseudobulk snRNA-seq objects absent -- ",
          "Fig 4d/4e and Supp 2d read pre-computed CSVs instead. See DATA.md.")

invisible(list(LIPID_ROOT = .repo_root, FIGROOT = file.path(.repo_root, "analysis")))
