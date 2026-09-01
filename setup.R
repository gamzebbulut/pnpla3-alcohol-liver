## Point the analysis scripts at this repository.
##
## Every script resolves its inputs from LIPID_ROOT and FIGROOT, with all other
## paths derived from those two by default. Setting them here is what lets the
## scripts run unmodified from a fresh clone.
##
##   source("setup.R")   # then source any script under analysis/

repo <- normalizePath(".", mustWork = TRUE)
Sys.setenv(LIPID_ROOT = repo)
Sys.setenv(FIGROOT    = file.path(repo, "analysis"))

## The snRNA-seq pseudobulk objects are not distributed with this repo (DATA.md).
## Scripts that need them will stop with a clear message rather than half-run.
if (!dir.exists(file.path(repo, "Data/pseudobulk_snRNA-seq_PNPLA3_alc")))
  message("note: pseudobulk snRNA-seq objects absent -- ",
          "Fig 4d/4e and Supp 2d read pre-computed CSVs instead. See DATA.md.")

invisible(list(LIPID_ROOT = repo, FIGROOT = file.path(repo, "analysis")))
