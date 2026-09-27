library(DESeq2)
library(pheatmap)

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

de_file <- file.path(
  project_dir,
  "results",
  "BRCA_DESeq2_results.tsv"
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

# Estimate normalization factors
dds <- estimateSizeFactors(dds)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Read differential expression results
res <- read.delim(
  de_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Remove missing adjusted p-values
res <- res[!is.na(res$padj), ]

# Select top 50 genes by adjusted p-value
res <- res[order(res$padj), ]

top_genes <- head(res$gene_id, 50)

# Extract transformed expression values
heatmap_matrix <- assay(vsd)[top_genes, ]

# Convert gene IDs to gene names when available
gene_names <- res$gene_name[match(top_genes, res$gene_id)]

valid_names <- !is.na(gene_names) & gene_names != ""

rownames(heatmap_matrix)[valid_names] <- gene_names[valid_names]

# Create sample annotation
annotation_col <- data.frame(
  Sample_Type = metadata$sample_type
)

rownames(annotation_col) <- metadata$file_id

# Generate heatmap
pdf(
  file.path(
    figure_dir,
    "Heatmap_Top50_DEGs.pdf"
  ),
  width = 10,
  height = 12
)

pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

# Save heatmap matrix
write.table(
  heatmap_matrix,
  file.path(
    project_dir,
    "results",
    "Heatmap_Top50_expression.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

cat("\nHeatmap created successfully.\n")
cat("Genes displayed:", nrow(heatmap_matrix), "\n")
cat(
  "Heatmap:",
  file.path(
    figure_dir,
    "Heatmap_Top50_DEGs.pdf"
  ),
  "\n"
)
cat(
  "Expression matrix:",
  file.path(
    project_dir,
    "results",
    "Heatmap_Top50_expression.tsv"
  ),
  "\n"
)

