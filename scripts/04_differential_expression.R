library(DESeq2)

project_dir <- "/mnt/c/Users/melli/Desktop/TCGA-BRCA-Gene-Expression-Analysis"

counts_file <- file.path(
  project_dir,
  "results",
  "BRCA_filtered_counts_matrix.tsv"
)

metadata_file <- file.path(
  project_dir,
  "results",
  "BRCA_sample_metadata.tsv"
)

output_file <- file.path(
  project_dir,
  "results",
  "BRCA_DESeq2_results.tsv"
)

# Read filtered count matrix
counts_df <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

gene_ids <- counts_df$gene_id
gene_names <- counts_df$gene_name

counts <- counts_df[, -(1:2), drop = FALSE]

rownames(counts) <- gene_ids

counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# Read metadata
metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Use unique GDC file IDs
rownames(metadata) <- metadata$file_id
colnames(counts) <- metadata$file_id

# Define experimental groups
metadata$sample_type <- factor(
  metadata$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ sample_type
)

# Run differential expression analysis
dds <- DESeq(dds)

# Extract Primary Tumor vs Solid Tissue Normal results
res <- results(
  dds,
  contrast = c(
    "sample_type",
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)

# Add gene information
res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

gene_info <- data.frame(
  gene_id = gene_ids,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

res_df <- merge(
  gene_info,
  res_df,
  by = "gene_id",
  all.y = TRUE,
  sort = FALSE
)

# Reorder columns
res_df <- res_df[, c(
  "gene_id",
  "gene_name",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

# Sort by adjusted p-value
res_df <- res_df[
  order(res_df$padj, na.last = TRUE),
]

# Save results
write.table(
  res_df,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Summary
cat("\nDESeq2 analysis completed successfully.\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("Results file:", output_file, "\n\n")

cat("Significant genes (padj < 0.05):",
    sum(res_df$padj < 0.05, na.rm = TRUE), "\n")

cat("Upregulated genes (padj < 0.05, log2FC > 1):",
    sum(res_df$padj < 0.05 &
        res_df$log2FoldChange > 1,
        na.rm = TRUE), "\n")

cat("Downregulated genes (padj < 0.05, log2FC < -1):",
    sum(res_df$padj < 0.05 &
        res_df$log2FoldChange < -1,
        na.rm = TRUE), "\n")
