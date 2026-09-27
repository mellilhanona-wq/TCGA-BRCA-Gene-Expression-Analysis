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

# Read filtered count matrix
counts_df <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Extract gene IDs
gene_ids <- counts_df$gene_id

# Extract count columns
counts <- counts_df[, -(1:2)]

# Set gene IDs as row names
rownames(counts) <- gene_ids

# Convert to numeric matrix
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# Read sample metadata
metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Use sample IDs as metadata row names
rownames(metadata) <- metadata$file_id
colnames(counts) <- metadata$file_id

# Define the experimental condition
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

# Print setup information
cat("DESeq2 object created successfully.\n")
cat("Genes:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n\n")

cat("Sample groups:\n")
print(table(dds$sample_type))

cat("\nReference level:", levels(dds$sample_type)[1], "\n")
cat("Comparison level:", levels(dds$sample_type)[2], "\n")
