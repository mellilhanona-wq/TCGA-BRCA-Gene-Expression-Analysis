# TCGA-BRCA Gene Expression Analysis

## Overview

This project demonstrates a cancer bioinformatics workflow for analyzing RNA-seq gene expression data from The Cancer Genome Atlas (TCGA) Breast Cancer dataset (TCGA-BRCA).

The analysis compares gene expression between primary breast tumor samples and solid tissue normal samples using R and Bioconductor tools.

The workflow includes:

- RNA-seq expression data organization
- Sample metadata preparation
- Gene-level count matrix construction
- Low-count gene filtering
- Principal Component Analysis (PCA)
- Differential gene expression analysis using DESeq2
- Identification of the top differentially expressed genes
- Heatmap visualization
- Volcano plot visualization
- Export of analysis results for reproducibility

## Biological Objective

The main objective is to identify genes whose expression differs between breast tumor and normal tissue.

Differential gene expression analysis can help identify genes and biological processes associated with cancer and provides a foundation for downstream functional and pathway analysis.

## Dataset

The project uses RNA-seq data from TCGA-BRCA obtained through the Genomic Data Commons (GDC).

The original project query included:

- Primary Tumor: 1,111 samples
- Solid Tissue Normal: 113 samples
- Total: 1,224 samples

The raw GDC data are not included in this repository because of their large size.

## Analysis Workflow

```text
TCGA/GDC RNA-seq Data
        ↓
Sample Metadata
        ↓
Gene Count Matrix
        ↓
Low-count Filtering
        ↓
Exploratory Analysis
        ├── PCA
        ↓
DESeq2 Differential Expression
        ├── Volcano Plot
        ├── Top Differentially Expressed Genes
        ↓
Heatmap of Top Genes
```

## Tools and Technologies

- R
- Bioconductor
- DESeq2
- TCGAbiolinks
- ggplot2
- pheatmap
- Linux / WSL
- Git / GitHub

## Project Structure

```text
TCGA-BRCA-Gene-Expression-Analysis/
│
├── figures/
│   ├── Heatmap_Top50_DEGs.pdf
│   ├── PCA_Tumor_vs_Normal.pdf
│   └── Volcano_DESeq2_Tumor_vs_Normal.pdf
│
├── results/
│   ├── BRCA_DESeq2_results.tsv
│   ├── BRCA_sample_metadata.tsv
│   ├── Heatmap_Top50_expression.tsv
│   ├── PCA_coordinates.tsv
│   └── Top50_DE_genes.tsv
│
├── scripts/
│   ├── 01_build_expression_matrix.R
│   ├── 02_filter_counts.R
│   ├── 03_deseq2_setup.R
│   ├── 04_differential_expression.R
│   ├── 05_qc_pca.R
│   ├── 06_volcano_plot.R
│   └── 07_heatmap.R
│
├── MANIFEST.txt
├── TCGA_BRCA_1224_manifest.txt
├── TCGA_BRCA_1224_query.tsv
├── TCGA_BRCA_missing_414_manifest.txt
├── .gitignore
└── README.md
```

## Main Results

### Principal Component Analysis

PCA was used to explore global expression patterns and assess whether tumor and normal samples show separation based on their transcriptomic profiles.

See:

`figures/PCA_Tumor_vs_Normal.pdf`

### Differential Expression Analysis

DESeq2 identified strong differences in gene expression between primary tumor and solid tissue normal samples.

Among the most strongly upregulated genes were **COL10A1, MMP11, COL11A1, MMP13, and IBSP**, which are associated with extracellular matrix remodeling and tumor-associated tissue changes. Several cell-cycle and mitotic genes, including **NEK2, KIF4A, UBE2C, CDC25C, TPX2, CDK1, and NUF2**, were also strongly upregulated, consistent with increased proliferative activity in tumor samples.

The complete exported DESeq2 results are available in:

`results/BRCA_DESeq2_results.tsv`

### Biological Interpretation

The differential expression results show two prominent patterns: increased expression of extracellular matrix/remodeling genes and increased expression of genes involved in cell-cycle progression and mitosis.

For example, **MMP11, COL10A1, COL11A1, and MMP13** showed strong upregulation, while **NEK2, UBE2C, CDC25C, TPX2, and CDK1** showed strong upregulation and very small adjusted p-values. These findings are consistent with biological processes commonly studied in breast cancer, including tumor-associated extracellular matrix remodeling and increased cellular proliferation.

These observations represent transcriptomic associations from the tumor-versus-normal comparison and would require downstream pathway enrichment and additional experimental validation to determine their functional significance.

### Top Differentially Expressed Genes

The top 50 differentially expressed genes were extracted for visualization and further interpretation.

See:

`results/Top50_DE_genes.tsv`

### Heatmap

A heatmap was generated to visualize expression patterns of the top 50 differentially expressed genes across samples.

See:

`figures/Heatmap_Top50_DEGs.pdf`

### Volcano Plot

A volcano plot was generated to visualize the relationship between statistical significance and magnitude of gene-expression changes.

See:

`figures/Volcano_DESeq2_Tumor_vs_Normal.pdf`
## Reproducibility

The analysis scripts are provided in the `scripts/` directory.

The scripts are numbered according to the analysis workflow so that the project can be followed from data preparation through visualization.

Large raw data files and expression matrices are intentionally excluded from the repository.

The provided manifest and query files document the TCGA/GDC data retrieval used for the project.

## Skills Demonstrated

This project demonstrates practical experience with:

- Cancer transcriptomics
- RNA-seq data analysis
- Gene expression matrices
- Metadata handling
- Data preprocessing
- Differential expression analysis
- Statistical analysis
- PCA
- Data visualization
- R programming
- Bioconductor
- Linux/WSL
- Reproducible bioinformatics workflows
- Git and GitHub

## Author

**Hanane Mellil**

Bachelor's student in Biotechnology Engineering

Interested in Bioinformatics, Computational Biology, Genomics, and Cancer Research.
