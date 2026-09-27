project_dir <- "/mnt/c/Users/melli/Desktop/TCGA-BRCA-Gene-Expression-Analysis"

results_file <- file.path(
  project_dir,
  "results",
  "BRCA_DESeq2_results.tsv"
)

figure_dir <- file.path(project_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE)

# Read DESeq2 results
res <- read.delim(
  results_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Remove rows with missing values
res <- res[
  !is.na(res$log2FoldChange) &
  !is.na(res$padj),
]

# Define significance categories
res$significance <- "Not significant"

res$significance[
  res$padj < 0.05 &
  res$log2FoldChange > 1
] <- "Upregulated"

res$significance[
  res$padj < 0.05 &
  res$log2FoldChange < -1
] <- "Downregulated"

# Avoid -log10(0)
res$plot_padj <- pmax(res$padj, .Machine$double.xmin)

# Create PDF
pdf(
  file.path(
    figure_dir,
    "Volcano_DESeq2_Tumor_vs_Normal.pdf"
  ),
  width = 8,
  height = 6
)

plot(
  res$log2FoldChange,
  -log10(res$plot_padj),
  pch = 20,
  cex = 0.6,
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "TCGA-BRCA Differential Expression"
)

abline(
  v = c(-1, 1),
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)

legend(
  "topright",
  legend = c(
    "Downregulated",
    "Not significant",
    "Upregulated"
  ),
  pch = 20
)

dev.off()

cat("\nVolcano plot created successfully.\n")
cat(
  "Output:",
  file.path(
    figure_dir,
    "Volcano_DESeq2_Tumor_vs_Normal.pdf"
  ),
  "\n"
)

