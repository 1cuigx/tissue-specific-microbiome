# Tissue-specific microbiome of Aiptasia

Source data and analysis code for the manuscript:

> **Algal symbiosis adaptively restructures the cnidarian microbiome to enhance holobiont metabolism**
> Guoxin Cui, David K. Ngugi, Migle Konciute, Si Cheng, Kit Xi Liew, Jianing Mi, Salim Al-Babili, Manuel Aranda.

This repository bundles every input data file, analysis script, and intermediate
result needed to regenerate the figures in the paper.

## Overview

The study dissects the tissue-resolved microbiome of the cnidarian model
Aiptasia under symbiotic vs. aposymbiotic states, combining:

- laser-microdissection of epidermis (Epi) and gastrodermis (Gas) tissues,
- tissue-resolved metatranscriptomics (kallisto → sleuth differential expression),
- KEGG pathway activity scoring (GSVA + limma),
- gnotobiotic hosts with 15N-labeled nitrite tracing of amino acids.

Each figure in the paper is reproduced by one R script in this repository.

## Repository layout

```
.
├── README.md
├── Fig1a_domain.{xlsx,R}                   # domain-level 16S composition
├── Fig1b_phyla.{xlsx,R}                    # phylum-level 16S composition
├── Fig1cde_phyla_boxplot.{xlsx,R}          # per-order boxplots + two-way ANOVA
├── Fig2abc_coding_taxonomy.{xlsx,R}        # mRNA-derived taxonomy boxplots
├── Fig2d_coding_phyla.{xlsx,R}             # mRNA phylum stacked bar + ANOVA
├── Fig3abcd_DE_PCA.Rmd                     # sleuth DE + PCA (4 pairwise contrasts)
├── Fig3abcd_RT_*.csv                       # sleuth results tables (4 files)
├── Fig3e_NiR.R                             # fold-change scatter w/ NirK/S hits
├── Fig3f_gsva_heatmap.R                    # GSVA heatmap + factorial limma
├── Fig3g_gsva_scatterplot.R                # GSVA Sym/Apo × Gas/Epi scatter
├── Fig4cde_AA.{xlsx,R}                     # 15N amino-acid tracing (rep 1)
├── FigS3_global_PCoA.R                     # global PCoA + PERMANOVA/PERMDISP
├── FigS5_AA_rep2.{xlsx,R}                  # 15N amino-acid tracing (rep 2)
├── kallisto/                               # kallisto pseudoalignment outputs
│   ├── ApoEpi1..4/  ApoGas1..4/
│   └── SymEpi1..4/  SymGas1..4/
├── aipBac.kegg.gs.all.RData                # KEGG gene-sets (`gs`) for GSVA
├── map_des.tsv                             # KEGG pathway descriptors
└── map_des_category.xlsx                   # descriptors + functional category
```

Figure-to-file mapping:

| Figure | Script                                   | Primary input(s)                                                   |
| -----: | ---------------------------------------- | ------------------------------------------------------------------ |
|   1a   | `Fig1a_domain.R`                         | `Fig1a_domain.xlsx`                                                |
|   1b   | `Fig1b_phyla.R`                          | `Fig1b_phyla.xlsx`                                                 |
|  1c–e  | `Fig1cde_phyla_boxplot.R`                | `Fig1cde_phyla_boxplot.xlsx`                                       |
|  2a–c  | `Fig2abc_coding_taxonomy.R`              | `Fig2abc_coding_taxonomy.xlsx`                                     |
|   2d   | `Fig2d_coding_phyla.R`                   | `Fig2d_coding_phyla.xlsx`                                          |
|  3a–d  | `Fig3abcd_sleuth_diff_expr.Rmd`          | `kallisto/`                                                        |
|   3e   | `Fig3e_expression_nitrite_reductase.R`   | `Fig3abcd_RT_*.csv` (from 3a–d)                                    |
|   3f   | `Fig3f_gsva_heatmap.R`                   | `kallisto/`, `aipBac.kegg.gs.all.RData`, `map_des.tsv`             |
|   3g   | `Fig3g_gsva_scatterplot.R`               | `kallisto/`, `aipBac.kegg.gs.all.RData`, `map_des_category.xlsx`   |
|  4c–e  | `Fig4cde_AA_15N_plot.R`                  | `Fig4cde_AA.xlsx`                                                  |
|  S3    | `FigS3_global_PCoA.R`                    | `kallisto/`                                                        |
|  S5    | `FigS5_AA_15N_plot_rep2.R`               | `FigS5_AA_rep2.xlsx`                                               |

## Requirements

- R ≥ 4.2
- CRAN packages: `tidyverse`, `readxl`, `ggpubr`, `car`, `emmeans`, `multcomp`,
  `rstatix`, `ggfortify`, `vegan`, `ggrepel`, `pheatmap`, `UpSetR`
- Bioconductor packages: `sleuth`, `GSVA`, `limma`, `ComplexHeatmap`

Quick install:

```r
install.packages(c("tidyverse", "readxl", "ggpubr", "car", "emmeans",
                   "multcomp", "rstatix", "ggfortify", "vegan", "ggrepel",
                   "pheatmap", "UpSetR"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("sleuth", "GSVA", "limma", "ComplexHeatmap"))
```

## Reproducing the figures

All scripts resolve paths relative to the repository root — clone and run
from the top-level directory:

```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>
```

Then, for any figure, either knit the Rmd or `source()` the `.R` file from R
launched in this directory, e.g.:

```r
source("Fig1cde_16S_phyla_boxplot.R")
rmarkdown::render("Fig3abcd_sleuth_diff_expr.Rmd")
```

`Fig3e_expression_nitrite_reductase.R` depends on the four
`Fig3abcd_RT_*.csv` tables; those are shipped in the repository and are also
regenerated by running `Fig3abcd_sleuth_diff_expr.Rmd`. All other scripts are
self-contained.

Each script writes its figure (PDF/SVG) and any accompanying Supplementary
Tables (CSV) to the working directory.

## Data availability

- **RNA-seq raw data**: NCBI Sequence Read Archive, BioProject accession
  `PRJNA631577`.
- **UHPLC-HR-MS metabolomics**: NIH Common Fund’s National Metabolomics Data
  Repository, accession `ST004292` (released upon publication).
- **Processed data and code**: this repository.

## Citation

If you use any material from this repository, please cite the paper above.

## Contact

- Guoxin Cui — guoxin.cui@kaust.edu.sa
- Manuel Aranda — manuel.aranda@kaust.edu.sa

## License

Released under the MIT License (see `LICENSE`).
