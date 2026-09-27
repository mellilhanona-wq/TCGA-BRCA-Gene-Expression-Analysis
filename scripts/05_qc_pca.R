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

figure_dir <- file.path(project_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE)

# Read filtered counts
counts_df <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

gene_ids <- counts_df$gene_id

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

# Define sample groups
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

# Run normalization and dispersion estimation
dds <- estimateSizeFactors(dds)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA
pca_data <- plotPCA(
  vsd,
  intgroup = "sample_type",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar"),
  1
)

# Save PCA coordinates
write.table(
  pca_data,
  file.path(project_dir, "results", "PCA_coordinates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Create PCA plot
pdf(
  file.path(figure_dir, "PCA_Tumor_vs_Normal.pdf"),
  width = 8,
  height = 6
)

plot(
  pca_data$PC1,
  pca_data$PC2,
  xlab = paste0("PC1: ", percent_var[1], "% variance"),
  ylab = paste0("PC2: ", percent_var[2], "% variance"),
  main = "PCA of TCGA-BRCA Samples",
  pch = 19
)

legend(
  "topright",
  legend = levels(metadata$sample_type),
  pch = 19
)

dev.off()

cat("\nPCA QC completed successfully.\n")
cat("PC1 variance:", percent_var[1], "%\n")
cat("PC2 variance:", percent_var[2], "%\n")
cat("PCA coordinates:", file.path(project_dir, "results", "PCA_coordinates.tsv"), "\n")
cat("PCA plot:", file.path(figure_dir, "PCA_Tumor_vs_Normal.pdf"), "\n")

