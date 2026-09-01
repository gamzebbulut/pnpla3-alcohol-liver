# PNPLA3 I148M and alcohol in obese mouse liver

Analysis code and derived data for the manuscript. Every figure in the paper is
reproduced by the scripts here, from the data in this repository — no raw files
required.

## Reproduce the figures

```r
source("run_all.R")
```

Panels are written as SVG to `analysis/<Fig N>/nc_panels/`, and `sessionInfo.txt`
records the package versions of the run.

To rebuild one panel instead:

```r
source("setup.R")
source(file.path(Sys.getenv("FIGROOT"), "Fig 3/fig3_grid_nc.R"), chdir = TRUE)
```

## Layout

```
analysis/           figure scripts, one folder per figure, plus nc_style.R
  nc_style.R        shared colours, type sizes, contrast names, the SVG writer
GB_metabolomics_20260628/
                    the gene/metabolite integration, one folder per contrast
source_data/        Source Data files published with the paper
figures/            rendered panels, committed so the repo is checkable
docs/               the figure manifest and pipeline notes
setup.R             points the scripts at this repo
run_all.R           rebuilds everything in dependency order
```

`docs/MANIFEST_all_figures.md` maps every panel to the script that makes it,
the file it writes and its size in millimetres. Start there.

## Design

`nc_style.R` is a contract, not a helper: it holds every shared colour, type
size and contrast label, and every panel script sources it. Changing a colour
or a label in that one file changes it everywhere, which is how the figures
stay consistent across five main figures and three supplements.

Paths come from environment variables with sensible defaults, so nothing needs
editing to run somewhere else — `setup.R` sets the two that matter.

## Data

See `DATA.md`. Short version: the derived tables are here, the raw
metabolomics and snRNA-seq are not, and the figures do not need them.

## Statistics

Metabolites: two-tailed Welch t-tests with Storey q-values. Lipids:
limma-moderated statistics on log2 abundances via `lipidr`; lipid set enrichment
by Monte Carlo permutation, minimum set size 4, 10,000 permutations. Genes:
edgeR quasi-likelihood F tests on pseudobulk, one column per mouse, with
transcriptome-wide FDR. Pathway integration: cell-type-proportion-weighted
Stouffer z across 17 cell types for the gene arm, ssGSEA and limma for the
metabolite arm, combined by intersection-union on the larger of the two
p-values. Significance is Benjamini-Hochberg FDR < 0.10 throughout unless a
figure legend says otherwise.
