project_dir <- "/mnt/c/Users/melli/Desktop/TCGA-BRCA-Gene-Expression-Analysis"

gdc_dir <- file.path(project_dir, "GDC_download")
manifest_file <- file.path(project_dir, "TCGA_BRCA_1224_manifest.txt")
query_file <- file.path(project_dir, "TCGA_BRCA_1224_query.tsv")

output_matrix <- file.path(project_dir, "results", "BRCA_counts_matrix.tsv")
output_metadata <- file.path(project_dir, "results", "BRCA_sample_metadata.tsv")

dir.create(file.path(project_dir, "results"), showWarnings = FALSE)

# Find all downloaded RNA-seq files
files <- list.files(
  gdc_dir,
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Number of downloaded files:", length(files), "\n")

# Read query metadata
metadata <- read.delim(
  query_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Match downloaded filenames to GDC metadata
file_names <- basename(files)

metadata <- metadata[match(file_names, metadata$file_name), ]

if (any(is.na(metadata$file_name))) {
  stop("Some downloaded files could not be matched to the GDC metadata.")
}

# Read one file to obtain gene information
first <- read.delim(
  files[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Keep gene-level rows only
first <- first[!grepl("^N_", first$gene_id), ]

gene_ids <- first$gene_id
gene_names <- first$gene_name

# Create count matrix
counts <- matrix(
  nrow = length(gene_ids),
  ncol = length(files)
)

rownames(counts) <- gene_ids
colnames(counts) <- metadata$sample.submitter_id

# Fill matrix one sample at a time
for (i in seq_along(files)) {

  cat("Reading", i, "of", length(files), "\n")

  dat <- read.delim(
    files[i],
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  dat <- dat[!grepl("^N_", dat$gene_id), ]

  counts[, i] <- dat$unstranded[match(gene_ids, dat$gene_id)]
}

# Convert to data frame
counts_df <- data.frame(
  gene_id = rownames(counts),
  gene_name = gene_names,
  counts,
  check.names = FALSE
)

# Save expression matrix
write.table(
  counts_df,
  output_matrix,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Prepare sample metadata
sample_metadata <- metadata[, c(
  "file_id",
  "file_name",
  "cases.submitter_id",
  "sample.submitter_id",
  "sample_type"
)]

write.table(
  sample_metadata,
  output_metadata,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nDONE\n")
cat("Expression matrix:", output_matrix, "\n")
cat("Sample metadata:", output_metadata, "\n")
