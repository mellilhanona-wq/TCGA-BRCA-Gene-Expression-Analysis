project_dir <- "/mnt/c/Users/melli/Desktop/TCGA-BRCA-Gene-Expression-Analysis"

input_file <- file.path(
  project_dir,
  "results",
  "BRCA_counts_matrix.tsv"
)

output_file <- file.path(
  project_dir,
  "results",
  "BRCA_filtered_counts_matrix.tsv"
)

counts_df <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

counts <- as.matrix(counts_df[, -(1:2)])

# Keep genes with at least 10 counts in at least 113 samples
min_count <- 10
min_samples <- 113

keep <- rowSums(counts >= min_count) >= min_samples

filtered_df <- counts_df[keep, ]

write.table(
  filtered_df,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Original genes:", nrow(counts_df), "\n")
cat("Genes retained:", nrow(filtered_df), "\n")
cat("Genes removed:", sum(!keep), "\n")
cat("Samples:", ncol(counts), "\n")
cat("Filtering rule: >=", min_count,
    "counts in >=", min_samples, "samples\n")
cat("Output:", output_file, "\n")
