# Data

Three tiers of data sit behind this paper. **Two of them are in this
repository.** The third is not, and does not need to be for the figures to
reproduce.

## Tier 0 — raw. Not in this repository.

| what | file |
|---|---|
| Metabolomics and lipidomics, as delivered by the platform | `Data_22_0018_VCU_PNPLA3 EtOH.xlsx` |
| Pseudobulk snRNA-seq, one object per cell type (17) | `pseudobulk_snRNA-seq_PNPLA3_alc/*.RDS` |

These are excluded by `.gitignore` and will be deposited on publication.
Until then they are available from the corresponding author on request.

## Tier 1 — derived. In this repository.

The tables the figure scripts actually read: per-metabolite fold-changes, the
lipid differential-abundance and set-enrichment results, the per-population
gene-level results, and the five per-contrast integration folders under
`GB_metabolomics_20260628/`.

**This tier is why the repository is self-sufficient.** Every panel is built
from these files, so the figures reproduce from a clone with no raw data
present. Tier 0 is only needed to rebuild tier 1 itself.

## Tier 2 — published. In this repository.

`source_data/` holds the Source Data files that accompany the paper, one per
panel that carries data. `figures/` holds the rendered SVGs, committed so the
repository can be checked without running anything.

## What cannot be rebuilt from a clone

Two upstream steps need tier 0 and will stop with a clear message rather than
half-run:

- `Fig 3/lipid_LSEA_Q3_vs_Q4.R` — reads the raw workbook and writes the lipid
  results. Its outputs are committed.
- `Fig 4/fig4_single_population.R` and `supp_lipid_story/supp_d_enzymes_nc.R`
  — read the snRNA-seq objects. Their outputs are committed, and the panels
  that display them read the CSVs, not the objects.

Everything else runs from the clone.
