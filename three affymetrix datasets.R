# Set your dataset folder
setwd("/Users/abarnageethanjali/Documents/set/combine datasets/affymetrix datasets ")
# 2. CREATE OUTPUT DIRECTORIES
output_dir <- "AUD_Fresh_Analysis_Results"

dir.create(output_dir, showWarnings = FALSE)

dir.create(file.path(output_dir, "01_Data_QC"),
           showWarnings = FALSE)
dir.create(file.path(output_dir, "02_Gene_Mapping"),
           showWarnings = FALSE)
dir.create(file.path(output_dir, "03_PCA_Before_ComBat"),
           showWarnings = FALSE)
dir.create(file.path(output_dir, "04_PCA_After_ComBat"),
           showWarnings = FALSE)
dir.create(file.path(output_dir, "05_Differential_Expression"),
           showWarnings = FALSE)
dir.create(file.path(output_dir, "06_Machine_Learning"),
           showWarnings = FALSE)
dir.create(file.path(output_dir, "07_PPI_CytoHubba"),
           showWarnings = FALSE)

# 3. CHECK DOWNLOADED FILES
files_needed <- c(
  "GSE44456_series_matrix.txt.gz",
  "GSE53808_series_matrix.txt.gz",
  "GSE62699-GPL571_series_matrix.txt.gz"
)

file_check <- data.frame(
  File = files_needed,
  Exists = file.exists(files_needed),
  Size_MB = round(
    file.info(files_needed)$size / (1024^2), 2
  )
)

print(file_check)

if (!all(file_check$Exists) ||
    any(is.na(file_check$Size_MB)) ||
    any(file_check$Size_MB <= 0)) {
  stop("One or more GEO files are missing or empty. Check the filenames and working directory.")
}

write.csv(
  file_check,
  file.path(output_dir, "01_Data_QC", "Downloaded_File_Check.csv"),
  row.names = FALSE
)

# 4. INSTALL BIOCONDUCTOR IF NECESSARY
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 5. INSTALL REQUIRED PACKAGES IF MISSING
cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "digest"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_packages <- c(
  "GEOquery",
  "Biobase",
  "limma",
  "sva",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "hugene10sttranscriptcluster.db",
  "hgu219.db",
  "hgu133a2.db"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# 6. LOAD LIBRARIES
library(GEOquery)
library(Biobase)
library(limma)
library(sva)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(hugene10sttranscriptcluster.db)
library(hgu219.db)
library(hgu133a2.db)

library(ggplot2)
library(ggrepel)
library(pheatmap)

cat("\nAll required libraries loaded.\n")

# 7. CLEAR OLD ANALYSIS OBJECTS
# This is a fresh script, so remove old objects from the workspace.
# Do not run this if you have unsaved work in the current session.
rm(list = setdiff(
  ls(),
  c("output_dir", "files_needed", "file_check")
))

# Reload essential paths and file-check objects after cleanup
output_dir <- "AUD_Fresh_Analysis_Results"

files_needed <- c(
  "GSE44456_series_matrix.txt.gz",
  "GSE53808_series_matrix.txt.gz",
  "GSE62699-GPL571_series_matrix.txt.gz"
)

# 8. LOAD GEO SERIES MATRIX FILES
gse44456 <- getGEO(
  filename = files_needed[1]
)

gse53808 <- getGEO(
  filename = files_needed[2]
)

gse62699 <- getGEO(
  filename = files_needed[3]
)

# 9. EXTRACT EXPRESSION MATRICES
expr44456_raw <- exprs(gse44456)
expr53808_raw <- exprs(gse53808)
expr62699_raw <- exprs(gse62699)

# 10. EXTRACT SAMPLE METADATA
pheno44456 <- pData(gse44456)
pheno53808 <- pData(gse53808)
pheno62699 <- pData(gse62699)

# 11. REPORT DIMENSIONS
dataset_summary <- data.frame(
  Dataset = c("GSE44456", "GSE53808", "GSE62699"),
  Probes = c(
    nrow(expr44456_raw),
    nrow(expr53808_raw),
    nrow(expr62699_raw)
  ),
  Samples = c(
    ncol(expr44456_raw),
    ncol(expr53808_raw),
    ncol(expr62699_raw)
  ),
  Metadata_Rows = c(
    nrow(pheno44456),
    nrow(pheno53808),
    nrow(pheno62699)
  )
)

print(dataset_summary)

write.csv(
  dataset_summary,
  file.path(output_dir, "01_Data_QC", "Dataset_Dimensions.csv"),
  row.names = FALSE
)

# 12. VERIFY SAMPLE ALIGNMENT
stopifnot(
  identical(rownames(pheno44456), colnames(expr44456_raw)),
  identical(rownames(pheno53808), colnames(expr53808_raw)),
  identical(rownames(pheno62699), colnames(expr62699_raw))
)

# 13. SAVE RAW SAMPLE METADATA
write.csv(
  pheno44456,
  file.path(output_dir, "01_Data_QC", "GSE44456_Sample_Metadata.csv")
)

write.csv(
  pheno53808,
  file.path(output_dir, "01_Data_QC", "GSE53808_Sample_Metadata.csv")
)

write.csv(
  pheno62699,
  file.path(output_dir, "01_Data_QC", "GSE62699_Sample_Metadata.csv")
)

# 14. SAVE RAW EXPRESSION MATRICES
# These are the expression values as supplied in the series matrices.
write.csv(
  data.frame(ProbeID = rownames(expr44456_raw), expr44456_raw,
             check.names = FALSE),
  file.path(output_dir, "01_Data_QC", "GSE44456_Raw_SeriesMatrix_Expression.csv")
)

write.csv(
  data.frame(ProbeID = rownames(expr53808_raw), expr53808_raw,
             check.names = FALSE),
  file.path(output_dir, "01_Data_QC", "GSE53808_Raw_SeriesMatrix_Expression.csv")
)

write.csv(
  data.frame(ProbeID = rownames(expr62699_raw), expr62699_raw,
             check.names = FALSE),
  file.path(output_dir, "01_Data_QC", "GSE62699_Raw_SeriesMatrix_Expression.csv")
)

# 15. CHECK EXPRESSION VALUES
expression_summary <- data.frame(
  Dataset = c("GSE44456", "GSE53808", "GSE62699"),
  Minimum = c(
    min(expr44456_raw, na.rm = TRUE),
    min(expr53808_raw, na.rm = TRUE),
    min(expr62699_raw, na.rm = TRUE)
  ),
  Median = c(
    median(expr44456_raw, na.rm = TRUE),
    median(expr53808_raw, na.rm = TRUE),
    median(expr62699_raw, na.rm = TRUE)
  ),
  Maximum = c(
    max(expr44456_raw, na.rm = TRUE),
    max(expr53808_raw, na.rm = TRUE),
    max(expr62699_raw, na.rm = TRUE)
  ),
  Missing_Values = c(
    sum(is.na(expr44456_raw)),
    sum(is.na(expr53808_raw)),
    sum(is.na(expr62699_raw))
  )
)

print(expression_summary)

write.csv(
  expression_summary,
  file.path(output_dir, "01_Data_QC", "Expression_Value_Summary.csv"),
  row.names = FALSE
)

cat("\nSTEP 1 COMPLETED.\n")
cat("Results saved in:", file.path(getwd(), output_dir), "\n")

# STEP 2: Inspect sample annotations

# GSE44456
cat("\n========== GSE44456 ==========\n")
print(colnames(pd44456))
print(head(pd44456[, intersect(
  c("title", "source_name_ch1", "characteristics_ch1",
    "characteristics_ch1.1", "characteristics_ch1.2"),
  colnames(pd44456)
), drop = FALSE], 10))

# GSE53808
cat("\n========== GSE53808 ==========\n")
print(colnames(pd53808))
print(head(pd53808[, intersect(
  c("title", "source_name_ch1", "characteristics_ch1",
    "characteristics_ch1.1", "characteristics_ch1.2"),
  colnames(pd53808)
), drop = FALSE], 10))

# GSE62699
cat("\n========== GSE62699 ==========\n")
print(colnames(pd62699))
print(head(pd62699[, intersect(
  c("title", "source_name_ch1", "characteristics_ch1",
    "characteristics_ch1.1", "characteristics_ch1.2"),
  colnames(pd62699)
), drop = FALSE], 10))

# Check which dataset and metadata objects currently exist
ls()

# Find objects related to the three GEO datasets
grep("44456|53808|62699|pd|meta",
     ls(),
     value = TRUE,
     ignore.case = TRUE)

# Inspect phenotype column names
cat("\n========== GSE44456 ==========\n")
print(colnames(pheno44456))
print(head(pheno44456, 3))

cat("\n========== GSE53808 ==========\n")
print(colnames(pheno53808))
print(head(pheno53808, 3))

cat("\n========== GSE62699 ==========\n")
print(colnames(pheno62699))
print(head(pheno62699, 3))

# Display likely sample-identification and diagnosis fields
for (nm in c("pheno44456", "pheno53808", "pheno62699")) {
  cat("\n\n==========", nm, "==========\n")
  x <- get(nm)
  
  candidate_cols <- grep(
    "title|source_name|characteristics|diagnos|group|tissue",
    colnames(x),
    ignore.case = TRUE,
    value = TRUE
  )
  
  print(unique(x[, candidate_cols, drop = FALSE]))
}

# ============================================================
# STEP 2A: Create and verify sample information
# ============================================================

library(dplyr)

# 1. Extract groups from the correct metadata columns

group44456 <- ifelse(
  grepl("alcoholic", pheno44456$`phenotype:ch1`,
        ignore.case = TRUE),
  "AUD",
  ifelse(
    grepl("control", pheno44456$`phenotype:ch1`,
          ignore.case = TRUE),
    "Control", NA
  )
)

group53808 <- ifelse(
  grepl("^control$", pheno53808$`disease state:ch1`,
        ignore.case = TRUE),
  "Control",
  ifelse(
    grepl("^alcoholic$", pheno53808$`disease state:ch1`,
          ignore.case = TRUE),
    "AUD",
    ifelse(
      grepl("HE alcoholic",
            pheno53808$`disease state:ch1`,
            ignore.case = TRUE),
      "HE_alcoholic", NA
    )
  )
)

group62699 <- ifelse(
  grepl("control", pheno62699$`diagnosis:ch1`,
        ignore.case = TRUE),
  "Control",
  ifelse(
    grepl("alcohol dependence",
          pheno62699$`diagnosis:ch1`,
          ignore.case = TRUE),
    "AUD", NA
  )
)

# 2. Build one metadata table per dataset

info44456 <- data.frame(
  Sample = rownames(pheno44456),
  Dataset = "GSE44456",
  Group = group44456,
  Tissue = pheno44456$`tissue:ch1`,
  Platform = pheno44456$platform_id,
  stringsAsFactors = FALSE
)

info53808 <- data.frame(
  Sample = rownames(pheno53808),
  Dataset = "GSE53808",
  Group = group53808,
  Tissue = pheno53808$`tissue:ch1`,
  Platform = pheno53808$platform_id,
  stringsAsFactors = FALSE
)

info62699 <- data.frame(
  Sample = rownames(pheno62699),
  Dataset = "GSE62699",
  Group = group62699,
  Tissue = pheno62699$`tissue type:ch1`,
  Platform = pheno62699$platform_id,
  stringsAsFactors = FALSE
)

# 3. Combine the metadata
sample_info_all <- rbind(
  info44456,
  info53808,
  info62699
)

# Save all sample labels before filtering
write.csv(
  sample_info_all,
  file.path(output_dir, "01_Data_QC",
            "Sample_Information_All_Samples.csv"),
  row.names = FALSE
)

# 4. Inspect the group counts
print(with(sample_info_all, table(Dataset, Group)))

# 5. Filter out HE alcoholic samples from GSE53808
sample_info <- subset(
  sample_info_all,
  Group %in% c("Control", "AUD")
)

sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Control", "AUD")
)

# 6. Check sample counts after filtering
cat("\nFinal sample counts by dataset and group:\n")
print(with(sample_info, table(Dataset, Group)))

cat("\nTotal retained samples:", nrow(sample_info), "\n")

# 7. Save the retained sample information
write.csv(
  sample_info,
  file.path(output_dir, "01_Data_QC",
            "Sample_Information_Filtered.csv"),
  row.names = FALSE
)

# 8. Stop if any diagnosis labels could not be assigned
stopifnot(!anyNA(sample_info_all$Group))
print(with(sample_info_all, table(Dataset, Group)))
print(with(sample_info, table(Dataset, Group)))


# ============================================================
# STEP 3: Align expression matrices and sample information
# ============================================================

# Filter GSE53808 to retain only Control and AUD samples
keep53808 <- sample_info$Sample[
  sample_info$Dataset == "GSE53808"
]

# Retain samples in the order specified by sample_info
expr44456 <- expr44456_raw[
  , sample_info$Sample[sample_info$Dataset == "GSE44456"],
  drop = FALSE
]

expr53808 <- expr53808_raw[
  , keep53808,
  drop = FALSE
]

expr62699 <- expr62699_raw[
  , sample_info$Sample[sample_info$Dataset == "GSE62699"],
  drop = FALSE
]

# Check dimensions
cat("GSE44456:", dim(expr44456), "\n")
cat("GSE53808:", dim(expr53808), "\n")
cat("GSE62699:", dim(expr62699), "\n")

# Confirm sample order matches the metadata
stopifnot(identical(
  colnames(expr44456),
  sample_info$Sample[sample_info$Dataset == "GSE44456"]
))

stopifnot(identical(
  colnames(expr53808),
  sample_info$Sample[sample_info$Dataset == "GSE53808"]
))

stopifnot(identical(
  colnames(expr62699),
  sample_info$Sample[sample_info$Dataset == "GSE62699"]
))

# Confirm the expected sample totals
stopifnot(ncol(expr44456) == 39)
stopifnot(ncol(expr53808) == 24)
stopifnot(ncol(expr62699) == 36)
stopifnot(nrow(sample_info) == 99)

# Save the filtered expression matrices
write.csv(
  data.frame(ProbeID = rownames(expr44456), expr44456,
             check.names = FALSE),
  file.path(output_dir, "01_Data_QC",
            "GSE44456_Filtered_Expression.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(ProbeID = rownames(expr53808), expr53808,
             check.names = FALSE),
  file.path(output_dir, "01_Data_QC",
            "GSE53808_Filtered_Expression.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(ProbeID = rownames(expr62699), expr62699,
             check.names = FALSE),
  file.path(output_dir, "01_Data_QC",
            "GSE62699_Filtered_Expression.csv"),
  row.names = FALSE
)

cat("\nSUCCESS: Expression matrices and sample metadata are aligned.\n")

# ============================================================
# STEP 4B: Probe-to-gene mapping for all three datasets
# ============================================================

library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)
library(hgu219.db)
library(hgu133a2.db)

mapping_dir <- file.path(output_dir, "02_Gene_Mapping")
dir.create(mapping_dir, showWarnings = FALSE, recursive = TRUE)

# Function to map probes and collapse duplicate gene symbols
map_probes_to_genes <- function(expr, annotation_db, dataset_name) {
  
  cat("\nMapping", dataset_name, "...\n")
  
  # Check expression matrix
  stopifnot(!is.null(rownames(expr)))
  stopifnot(!anyDuplicated(rownames(expr)))
  
  probe_ids <- rownames(expr)
  
  # Confirm annotation database supports PROBEID and SYMBOL
  available_keys <- AnnotationDbi::keytypes(annotation_db)
  available_columns <- AnnotationDbi::columns(annotation_db)
  
  if (!"PROBEID" %in% available_keys) {
    stop("PROBEID keytype unavailable for ", dataset_name)
  }
  
  if (!"SYMBOL" %in% available_columns) {
    stop("SYMBOL column unavailable for ", dataset_name)
  }
  
  # Keep only probe IDs recognized by the annotation package
  valid_keys <- AnnotationDbi::keys(
    annotation_db,
    keytype = "PROBEID"
  )
  
  matched_probes <- intersect(probe_ids, valid_keys)
  
  if (length(matched_probes) == 0) {
    stop(
      "No probe IDs matched for ", dataset_name,
      ". Check the expression matrix row names."
    )
  }
  
  # Retrieve probe-to-symbol annotations
  annotation <- AnnotationDbi::select(
    annotation_db,
    keys = matched_probes,
    keytype = "PROBEID",
    columns = "SYMBOL"
  )
  
  # Remove missing or blank symbols
  annotation <- annotation[
    !is.na(annotation$SYMBOL) &
      nzchar(trimws(annotation$SYMBOL)),
    ,
    drop = FALSE
  ]
  
  # Use one symbol per probe.
  # If a probe has multiple annotation records, retain the
  # first distinct non-empty symbol and document this rule.
  annotation <- unique(annotation[, c("PROBEID", "SYMBOL")])
  annotation <- annotation[!duplicated(annotation$PROBEID), ]
  
  # Match annotations to expression rows
  idx <- match(annotation$PROBEID, rownames(expr))
  keep <- !is.na(idx)
  
  annotation <- annotation[keep, , drop = FALSE]
  idx <- idx[keep]
  
  mapped_expr <- expr[idx, , drop = FALSE]
  rownames(mapped_expr) <- annotation$SYMBOL
  
  # Average expression across probes assigned to the same gene
  gene_expr <- rowsum(
    as.matrix(mapped_expr),
    group = rownames(mapped_expr),
    reorder = FALSE
  )
  
  probe_counts <- table(rownames(mapped_expr))
  gene_expr <- sweep(
    gene_expr,
    1,
    as.numeric(probe_counts[rownames(gene_expr)]),
    FUN = "/"
  )
  
  # Save mapping report for every original probe
  mapping_report <- data.frame(
    ProbeID = probe_ids,
    GeneSymbol = NA_character_,
    stringsAsFactors = FALSE
  )
  
  report_idx <- match(annotation$PROBEID, mapping_report$ProbeID)
  mapping_report$GeneSymbol[report_idx] <- annotation$SYMBOL
  
  write.csv(
    mapping_report,
    file.path(mapping_dir,
              paste0(dataset_name, "_Probe_to_Gene_Mapping.csv")),
    row.names = FALSE
  )
  
  cat("Original probes:", nrow(expr), "\n")
  cat("Mapped probes:", nrow(mapped_expr), "\n")
  cat("Unique gene symbols:", nrow(gene_expr), "\n")
  cat("Unmapped probes:", sum(is.na(mapping_report$GeneSymbol)), "\n")
  
  return(list(
    expression = gene_expr,
    report = mapping_report
  ))
}

# Run the appropriate annotation package for each platform
res44456 <- map_probes_to_genes(
  expr44456,
  hugene10sttranscriptcluster.db,
  "GSE44456"
)

res53808 <- map_probes_to_genes(
  expr53808,
  hgu219.db,
  "GSE53808"
)

res62699 <- map_probes_to_genes(
  expr62699,
  hgu133a2.db,
  "GSE62699"
)

# Extract gene-level expression matrices
gene44456 <- res44456$expression
gene53808 <- res53808$expression
gene62699 <- res62699$expression

# Save gene-level expression matrices
write.csv(
  data.frame(GeneSymbol = rownames(gene44456),
             gene44456, check.names = FALSE),
  file.path(mapping_dir, "GSE44456_GeneLevel_Expression.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(GeneSymbol = rownames(gene53808),
             gene53808, check.names = FALSE),
  file.path(mapping_dir, "GSE53808_GeneLevel_Expression.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(GeneSymbol = rownames(gene62699),
             gene62699, check.names = FALSE),
  file.path(mapping_dir, "GSE62699_GeneLevel_Expression.csv"),
  row.names = FALSE
)

# Save mapping summary
mapping_summary <- data.frame(
  Dataset = c("GSE44456", "GSE53808", "GSE62699"),
  Original_Probes = c(
    nrow(expr44456), nrow(expr53808), nrow(expr62699)
  ),
  Mapped_Probes = c(
    sum(!is.na(res44456$report$GeneSymbol)),
    sum(!is.na(res53808$report$GeneSymbol)),
    sum(!is.na(res62699$report$GeneSymbol))
  ),
  Unique_Gene_Symbols = c(
    nrow(gene44456), nrow(gene53808), nrow(gene62699)
  ),
  Samples = c(
    ncol(gene44456), ncol(gene53808), ncol(gene62699)
  )
)

write.csv(
  mapping_summary,
  file.path(mapping_dir, "Gene_Mapping_Summary.csv"),
  row.names = FALSE
)

print(mapping_summary)

cat("\nSTEP 4 COMPLETED: Gene-level matrices saved to:\n")
cat(mapping_dir, "\n")
dim(gene44456)
dim(gene53808)
dim(gene62699)

stopifnot(ncol(gene44456) == 39)
stopifnot(ncol(gene53808) == 24)
stopifnot(ncol(gene62699) == 36)

stopifnot(!anyDuplicated(rownames(gene44456)))
stopifnot(!anyDuplicated(rownames(gene53808)))
stopifnot(!anyDuplicated(rownames(gene62699)))


# Find common gene symbols across all three datasets
common_genes <- Reduce(
  intersect,
  list(
    rownames(gene44456),
    rownames(gene53808),
    rownames(gene62699)
  )
)

# Sort alphabetically for easier reading
common_genes <- sort(common_genes)

# Display the total number of common genes
cat("Total common genes across all three datasets:",
    length(common_genes), "\n\n")

# Display every common gene name
print(common_genes)
# Create output folder path
common_gene_file <- file.path(
  output_dir,
  "02_Gene_Mapping",
  "Common_Genes_All_Three_Datasets.csv"
)

# Save common gene names
write.csv(
  data.frame(Gene_Symbol = common_genes),
  file = common_gene_file,
  row.names = FALSE
)

# Confirm the saved file and gene count
cat("Total common genes saved:", length(common_genes), "\n")
cat("File location:", normalizePath(common_gene_file), "\n")

# ============================================================
# STEP 1: COMBINE THE 12,373 COMMON GENES
# AUD study: GSE44456, GSE53808, GSE62699
# ============================================================

# Create separate output folders
pca_dir <- file.path(output_dir, "03_PCA_ComBat")
dge_dir <- file.path(output_dir, "04_Differential_Gene_Expression")

dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)

# Check common gene count
stopifnot(length(common_genes) == 12373)

# Check that common genes exist in all three matrices
stopifnot(
  all(common_genes %in% rownames(gene44456)),
  all(common_genes %in% rownames(gene53808)),
  all(common_genes %in% rownames(gene62699))
)

# Subset copies of the original matrices
common44456 <- gene44456[common_genes, , drop = FALSE]
common53808 <- gene53808[common_genes, , drop = FALSE]
common62699 <- gene62699[common_genes, , drop = FALSE]

# Combine samples across datasets
expr_combined <- cbind(
  common44456,
  common53808,
  common62699
)

# Verify dimensions and sample names
cat("Combined matrix dimensions:\n")
print(dim(expr_combined))

cat("\nSamples by dataset:\n")
print(c(
  GSE44456 = ncol(common44456),
  GSE53808 = ncol(common53808),
  GSE62699 = ncol(common62699)
))

stopifnot(
  nrow(expr_combined) == 12373,
  ncol(expr_combined) == 99,
  !anyDuplicated(rownames(expr_combined)),
  !anyDuplicated(colnames(expr_combined)),
  all(is.finite(expr_combined))
)

cat("\nStep 1 passed: 12,373 genes and 99 samples.\n")

# ============================================================
# STEP 2: ALIGN SAMPLE METADATA
# ============================================================

# Rebuild metadata in the same order as the combined matrix
sample_info_pca <- sample_info[
  match(colnames(expr_combined), sample_info$SampleID),
  ,
  drop = FALSE
]

# Check for missing or incorrectly matched samples
stopifnot(
  !anyNA(sample_info_pca$SampleID),
  identical(
    as.character(sample_info_pca$SampleID),
    colnames(expr_combined)
  ),
  !anyNA(sample_info_pca$Group),
  !anyNA(sample_info_pca$Dataset)
)

# Set the correct reference levels
sample_info_pca$Group <- factor(
  sample_info_pca$Group,
  levels = c("Control", "AUD")
)

sample_info_pca$Dataset <- factor(
  sample_info_pca$Dataset,
  levels = c("GSE44456", "GSE53808", "GSE62699")
)

# Verify sample counts
cat("Sample counts by dataset and group:\n")
print(with(sample_info_pca, table(Dataset, Group)))

cat("\nTotal samples:", nrow(sample_info_pca), "\n")
cat("Total genes:", nrow(expr_combined), "\n")

stopifnot(
  nrow(sample_info_pca) == 99,
  sum(sample_info_pca$Group == "Control") == 46,
  sum(sample_info_pca$Group == "AUD") == 53
)

cat("\nMetadata alignment verified successfully.\n")

# Diagnose sample ID mismatch

cat("Expression matrix sample names (first 10):\n")
print(head(colnames(expr_combined), 10))

cat("\nMetadata sample IDs (first 10):\n")
print(head(sample_info$SampleID, 10))

cat("\nNumber of expression samples:",
    ncol(expr_combined), "\n")
cat("Number of metadata rows:",
    nrow(sample_info), "\n")

cat("\nExpression samples missing from metadata:\n")
print(head(setdiff(colnames(expr_combined), sample_info$SampleID), 20))

cat("\nMetadata samples missing from expression matrix:\n")
print(head(setdiff(sample_info$SampleID, colnames(expr_combined)), 20))

cat("\nDuplicate expression sample names:",
    sum(duplicated(colnames(expr_combined))), "\n")
cat("Duplicate metadata sample IDs:",
    sum(duplicated(sample_info$SampleID)), "\n")

# Inspect metadata column names and structure
print(names(sample_info))
str(sample_info)
head(sample_info, 3)


# Align metadata to expression matrix column order
sample_info_pca <- sample_info[
  match(colnames(expr_combined), sample_info$Sample),
  ,
  drop = FALSE
]

# Verify exact sample alignment
stopifnot(
  !anyNA(sample_info_pca$Sample),
  identical(
    as.character(sample_info_pca$Sample),
    colnames(expr_combined)
  ),
  !anyNA(sample_info_pca$Group),
  !anyNA(sample_info_pca$Dataset)
)

# Set reference levels
sample_info_pca$Group <- factor(
  sample_info_pca$Group,
  levels = c("Control", "AUD")
)

sample_info_pca$Dataset <- factor(
  sample_info_pca$Dataset,
  levels = c("GSE44456", "GSE53808", "GSE62699")
)

# Confirm sample counts
cat("Sample counts by dataset and group:\n")
print(with(sample_info_pca, table(Dataset, Group)))

cat("\nTotal samples:", nrow(sample_info_pca), "\n")
cat("Total genes:", nrow(expr_combined), "\n")

stopifnot(
  nrow(sample_info_pca) == 99,
  sum(sample_info_pca$Group == "Control") == 46,
  sum(sample_info_pca$Group == "AUD") == 53
)

cat("\nSUCCESS: All samples are correctly aligned.\n")

library(ggplot2)

# PCA using samples as rows
pca_before <- prcomp(
  t(expr_combined),
  center = TRUE,
  scale. = FALSE
)

# Calculate variance explained
variance_before <- 100 *
  (pca_before$sdev^2 / sum(pca_before$sdev^2))

# Prepare plotting data in verified sample order
pca_before_df <- data.frame(
  Sample = rownames(pca_before$x),
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  Dataset = sample_info_pca$Dataset,
  Group = sample_info_pca$Group
)

# Plot
p1 <- ggplot(
  pca_before_df,
  aes(x = PC1, y = PC2, color = Dataset, shape = Group)
) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA Before ComBat",
    subtitle = "AUD study — 12,373 common genes",
    x = paste0("PC1 (", round(variance_before[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_before[2], 2), "%)"),
    color = "Dataset",
    shape = "Group"
  ) +
  theme_bw(base_size = 12)

print(p1)

# Save the plot and PCA results
ggsave(
  file.path(pca_dir, "PCA_Before_ComBat.png"),
  plot = p1,
  width = 9,
  height = 7,
  dpi = 300
)

write.csv(
  pca_before_df,
  file.path(pca_dir, "PCA_Before_ComBat_Coordinates.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Component = paste0("PC", seq_along(variance_before)),
    Variance_Explained_Percent = variance_before
  ),
  file.path(pca_dir, "PCA_Before_ComBat_Variance.csv"),
  row.names = FALSE
)

cat("PC1 variance explained:", round(variance_before[1], 2), "%\n")
cat("PC2 variance explained:", round(variance_before[2], 2), "%\n")
cat("PCA outputs saved in:", pca_dir, "\n")


# ============================================================
# STEP 4: COMBAT BATCH CORRECTION
# ============================================================

if (!requireNamespace("sva", quietly = TRUE)) {
  stop("Package 'sva' is not installed. Install it before continuing.")
}

library(sva)

# Verify alignment again
stopifnot(
  identical(colnames(expr_combined),
            as.character(sample_info_pca$Sample)),
  !anyNA(sample_info_pca$Group),
  !anyNA(sample_info_pca$Dataset),
  all(is.finite(expr_combined))
)

# Define factors
sample_info_pca$Group <- factor(
  sample_info_pca$Group,
  levels = c("Control", "AUD")
)

sample_info_pca$Dataset <- factor(
  sample_info_pca$Dataset,
  levels = c("GSE44456", "GSE53808", "GSE62699")
)

# Model protects the AUD/control group term during correction
mod_combat <- model.matrix(
  ~ Group,
  data = sample_info_pca
)

# Apply ComBat
expr_combat <- sva::ComBat(
  dat = as.matrix(expr_combined),
  batch = sample_info_pca$Dataset,
  mod = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)

# Verify corrected matrix
stopifnot(
  identical(dim(expr_combat), dim(expr_combined)),
  identical(rownames(expr_combat), rownames(expr_combined)),
  identical(colnames(expr_combat), colnames(expr_combined)),
  all(is.finite(expr_combat))
)

cat("ComBat completed successfully.\n")
cat("Corrected matrix dimensions:\n")
print(dim(expr_combat))

# Save without overwriting the original matrices
write.csv(
  expr_combat,
  file.path(pca_dir, "Common_Gene_Expression_ComBat.csv"),
  quote = FALSE
)

cat("Corrected expression matrix saved in:", pca_dir, "\n")

# ============================================================
# STEP 5: PCA AFTER COMBAT
# ============================================================

library(ggplot2)

# Verify corrected matrix and metadata alignment
stopifnot(
  nrow(expr_combat) == 12373,
  ncol(expr_combat) == 99,
  identical(colnames(expr_combat),
            as.character(sample_info_pca$Sample)),
  all(is.finite(expr_combat))
)

# Perform PCA on the corrected expression matrix
pca_after <- prcomp(
  t(expr_combat),
  center = TRUE,
  scale. = FALSE
)

# Calculate variance explained
variance_after <- 100 *
  (pca_after$sdev^2 / sum(pca_after$sdev^2))

# Prepare PCA coordinates
pca_after_df <- data.frame(
  Sample = rownames(pca_after$x),
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Dataset = sample_info_pca$Dataset,
  Group = sample_info_pca$Group
)

# Generate PCA plot
p2 <- ggplot(
  pca_after_df,
  aes(x = PC1, y = PC2, color = Dataset, shape = Group)
) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA After ComBat",
    subtitle = "AUD study — 12,373 common genes",
    x = paste0("PC1 (", round(variance_after[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_after[2], 2), "%)"),
    color = "Dataset",
    shape = "Group"
  ) +
  theme_bw(base_size = 12)

print(p2)

# Save PCA plot
ggsave(
  filename = file.path(pca_dir, "PCA_After_ComBat.png"),
  plot = p2,
  width = 9,
  height = 7,
  dpi = 300
)

# Save PCA coordinates
write.csv(
  pca_after_df,
  file.path(pca_dir, "PCA_After_ComBat_Coordinates.csv"),
  row.names = FALSE
)

# Save variance explained
write.csv(
  data.frame(
    Component = paste0("PC", seq_along(variance_after)),
    Variance_Explained_Percent = variance_after
  ),
  file.path(pca_dir, "PCA_After_ComBat_Variance.csv"),
  row.names = FALSE
)

# Report results
cat("PC1 variance explained after ComBat:",
    round(variance_after[1], 2), "%\n")
cat("PC2 variance explained after ComBat:",
    round(variance_after[2], 2), "%\n")
cat("Combined PC1 + PC2:",
    round(sum(variance_after[1:2]), 2), "%\n")

cat("PCA after ComBat saved in:", pca_dir, "\n")


# ============================================================
# STEP 6: DIFFERENTIAL GENE EXPRESSION (AUD vs CONTROL)
# Input: ComBat-corrected expression matrix
# Genes: 12,373 common genes
# Samples: 99
# ============================================================

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Install Bioconductor package 'limma' before continuing.")
}

library(limma)
library(ggplot2)

# ------------------------------------------------------------
# 1. Validate inputs and sample alignment
# ------------------------------------------------------------

stopifnot(
  nrow(expr_combat) == 12373,
  ncol(expr_combat) == 99,
  identical(
    colnames(expr_combat),
    as.character(sample_info_pca$Sample)
  ),
  all(is.finite(expr_combat))
)

sample_info_dge <- sample_info_pca

sample_info_dge$Group <- factor(
  sample_info_dge$Group,
  levels = c("Control", "AUD")
)

stopifnot(
  identical(
    levels(sample_info_dge$Group),
    c("Control", "AUD")
  ),
  sum(sample_info_dge$Group == "Control") == 46,
  sum(sample_info_dge$Group == "AUD") == 53
)

cat("Sample counts:\n")
print(table(sample_info_dge$Group))

# ------------------------------------------------------------
# 2. Build design matrix
# Control is the reference; GroupAUD is the AUD effect
# ------------------------------------------------------------

design_dge <- model.matrix(
  ~ Group,
  data = sample_info_dge
)

print(colnames(design_dge))

stopifnot("GroupAUD" %in% colnames(design_dge))

# ------------------------------------------------------------
# 3. Fit limma model
# ------------------------------------------------------------

fit_dge <- lmFit(
  expr_combat,
  design_dge
)

fit_dge <- eBayes(fit_dge)

# Extract results for AUD versus Control
dge_results <- topTable(
  fit_dge,
  coef = "GroupAUD",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

dge_results$Gene_Symbol <- rownames(dge_results)

# Put gene symbol first
dge_results <- dge_results[
  ,
  c("Gene_Symbol", setdiff(names(dge_results), "Gene_Symbol"))
]

# ------------------------------------------------------------
# 4. Classify genes using both significance criteria
# ------------------------------------------------------------

dge_results$Paper_Screen <- with(
  dge_results,
  P.Value < 0.05 & abs(logFC) > 1
)

dge_results$FDR_Significant <- with(
  dge_results,
  adj.P.Val < 0.05
)

dge_results$FDR_and_FC <- with(
  dge_results,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

dge_results$Direction <- ifelse(
  dge_results$logFC > 0,
  "Higher in AUD",
  "Lower in AUD"
)

# ------------------------------------------------------------
# 5. Summarize DGE results
# ------------------------------------------------------------

paper_deg <- subset(dge_results, Paper_Screen)

fdr_deg <- subset(dge_results, FDR_Significant)

fdr_fc_deg <- subset(dge_results, FDR_and_FC)

cat("\n================ DGE SUMMARY ================\n")
cat("Total genes tested:", nrow(dge_results), "\n")

cat(
  "Nominal P < 0.05 and |logFC| > 1:",
  nrow(paper_deg), "\n"
)

cat(
  "FDR-adjusted P < 0.05:",
  nrow(fdr_deg), "\n"
)

cat(
  "FDR-adjusted P < 0.05 and |logFC| > 1:",
  nrow(fdr_fc_deg), "\n"
)

cat("\nPaper-screened genes by direction:\n")
print(table(paper_deg$Direction))

cat("\nFDR-significant genes by direction:\n")
print(table(fdr_deg$Direction))

# ------------------------------------------------------------
# 6. Save result tables
# ------------------------------------------------------------

write.csv(
  dge_results,
  file.path(dge_dir, "AUD_vs_Control_All_12373_Genes.csv"),
  row.names = FALSE
)

write.csv(
  paper_deg,
  file.path(dge_dir, "AUD_vs_Control_Paper_Screen_DEGs.csv"),
  row.names = FALSE
)

write.csv(
  fdr_deg,
  file.path(dge_dir, "AUD_vs_Control_FDR_Significant_DEGs.csv"),
  row.names = FALSE
)

write.csv(
  fdr_fc_deg,
  file.path(dge_dir, "AUD_vs_Control_FDR_and_FC_DEGs.csv"),
  row.names = FALSE
)


ls()


# Install ggplot2 only if it is not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

# Replace deg_pooled with your full pooled limma results object
volcano_data <- as.data.frame(deg_pooled)

# Keep gene identifiers
volcano_data$Gene <- rownames(volcano_data)

# Remove rows with missing values
volcano_data <- volcano_data[
  !is.na(volcano_data$logFC) &
    !is.na(volcano_data$adj.P.Val), ]

# Classify genes using the FDR threshold
volcano_data$Regulation <- "Not significant"

volcano_data$Regulation[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC > 0
] <- "Upregulated"

volcano_data$Regulation[
  volcano_data$adj.P.Val < 0.05 &
    volcano_data$logFC < 0
] <- "Downregulated"

# Avoid log(0) if any adjusted p-values are zero
volcano_data$FDR_for_plot <- pmax(
  volcano_data$adj.P.Val,
  .Machine$double.xmin
)

# Create volcano plot
p <- ggplot(
  volcano_data,
  aes(x = logFC, y = -log10(FDR_for_plot),
      color = Regulation)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "#D73027",
    "Downregulated" = "#4575B4",
    "Not significant" = "grey75"
  )) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes in AUD",
    subtitle = "Pooled limma analysis of three GEO datasets",
    x = "log2 Fold Change (AUD vs Control)",
    y = "-log10 Adjusted P-value (FDR)",
    color = "Gene regulation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p)

# Save a high-resolution image for your research paper
ggsave(
  "AUD_Pooled_Volcano_Plot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# Verify the number of significant genes
table(volcano_data$Regulation)



# ------------------------------------------------------------
# 8. Save DGE summary
# ------------------------------------------------------------

dge_summary <- data.frame(
  Analysis = c(
    "All genes tested",
    "Nominal P < 0.05 and |logFC| > 1",
    "FDR-adjusted P < 0.05",
    "FDR-adjusted P < 0.05 and |logFC| > 1"
  ),
  Gene_Count = c(
    nrow(dge_results),
    nrow(paper_deg),
    nrow(fdr_deg),
    nrow(fdr_fc_deg)
  )
)

write.csv(
  dge_summary,
  file.path(dge_dir, "AUD_vs_Control_DGE_Summary.csv"),
  row.names = FALSE
)

cat("\nDGE completed.\n")
cat("Results saved in:", dge_dir, "\n")

# ============================================================
# STEP 7: HEATMAP OF TOP 50 DEGs
# ============================================================

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Install the pheatmap package before continuing.")
}

library(pheatmap)

# 1. Select top 50 genes ranked by nominal P-value
top50_genes <- head(
  dge_results$Gene_Symbol[
    order(dge_results$P.Value, na.last = NA)
  ],
  50
)

stopifnot(length(top50_genes) == 50)

# 2. Extract corrected expression values
heatmap_matrix <- expr_combat[
  top50_genes,
  ,
  drop = FALSE
]

# 3. Z-score each gene across samples
# This highlights relative expression patterns per gene
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

# Remove genes with undefined scaled values, if any
heatmap_matrix_scaled <- heatmap_matrix_scaled[
  apply(heatmap_matrix_scaled, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

# 4. Prepare sample annotations in expression-matrix order
annotation_col <- data.frame(
  Dataset = sample_info_pca$Dataset,
  Group = sample_info_pca$Group
)

rownames(annotation_col) <- sample_info_pca$Sample

stopifnot(
  identical(colnames(heatmap_matrix_scaled),
            rownames(annotation_col))
)

# 5. Generate heatmap
pheatmap(
  heatmap_matrix_scaled,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  border_color = NA,
  fontsize = 9,
  fontsize_row = 7,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  main = "Top 50 DEGs: AUD vs Control",
  filename = file.path(
    dge_dir,
    "AUD_vs_Control_Top50_DEG_Heatmap.png"
  ),
  width = 12,
  height = 9
)

# 6. Save the selected gene list
write.csv(
  data.frame(Gene_Symbol = rownames(heatmap_matrix_scaled)),
  file.path(dge_dir, "Top50_DEG_Heatmap_Genes.csv"),
  row.names = FALSE
)

cat("Top 50 genes selected:", length(top50_genes), "\n")
cat("Genes displayed after scaling:",
    nrow(heatmap_matrix_scaled), "\n")
cat("Heatmap saved in:", dge_dir, "\n")

# ============================================================
# STEP 8: LASSO FEATURE SELECTION
# ============================================================

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Install the glmnet package before continuing.")
}

library(glmnet)

lasso_dir <- file.path(
  output_dir,
  "05_Machine_Learning",
  "01_LASSO"
)

dir.create(
  lasso_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 1. Validate input alignment
stopifnot(
  nrow(expr_combat) == 12373,
  ncol(expr_combat) == 99,
  identical(
    colnames(expr_combat),
    as.character(sample_info_pca$Sample)
  )
)

# 2. Prepare predictors: samples x genes
x <- t(expr_combat)

# Remove genes with zero variance or non-finite values
gene_variance <- apply(x, 2, var)
keep_genes <- is.finite(gene_variance) & gene_variance > 0
x <- x[, keep_genes, drop = FALSE]

# Binary outcome: Control = 0, AUD = 1
y <- ifelse(sample_info_pca$Group == "AUD", 1, 0)

stopifnot(
  nrow(x) == length(y),
  all(is.finite(x)),
  length(unique(y)) == 2
)

cat("Samples:", nrow(x), "\n")
cat("Genes entering LASSO:", ncol(x), "\n")
cat("Control samples:", sum(y == 0), "\n")
cat("AUD samples:", sum(y == 1), "\n")

# 3. Set reproducible seed and run 10-fold cross-validation
set.seed(123)

cv_lasso <- cv.glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance",
  standardize = TRUE
)

# 4. Save cross-validation plot
png(
  filename = file.path(lasso_dir, "LASSO_10Fold_CV.png"),
  width = 2400,
  height = 1800,
  res = 300
)
plot(cv_lasso)
dev.off()

# 5. Extract selected genes at lambda.min
coef_min <- coef(cv_lasso, s = "lambda.min")
coef_matrix <- as.matrix(coef_min)

lasso_table <- data.frame(
  Gene_Symbol = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix[, 1]),
  stringsAsFactors = FALSE
)

lasso_selected <- subset(
  lasso_table,
  Gene_Symbol != "(Intercept)" & Coefficient != 0
)

lasso_selected <- lasso_selected[
  order(abs(lasso_selected$Coefficient), decreasing = TRUE),
]

# 6. Also extract genes at the more conservative lambda.1se
coef_1se <- as.matrix(coef(cv_lasso, s = "lambda.1se"))

lasso_1se <- data.frame(
  Gene_Symbol = rownames(coef_1se),
  Coefficient = as.numeric(coef_1se[, 1]),
  stringsAsFactors = FALSE
)

lasso_1se <- subset(
  lasso_1se,
  Gene_Symbol != "(Intercept)" & Coefficient != 0
)

lasso_1se <- lasso_1se[
  order(abs(lasso_1se$Coefficient), decreasing = TRUE),
]

# 7. Save results
write.csv(
  lasso_selected,
  file.path(lasso_dir, "LASSO_Selected_Genes_LambdaMin.csv"),
  row.names = FALSE
)

write.csv(
  lasso_1se,
  file.path(lasso_dir, "LASSO_Selected_Genes_Lambda1SE.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Lambda_Min = cv_lasso$lambda.min,
    Lambda_1SE = cv_lasso$lambda.1se,
    Genes_LambdaMin = nrow(lasso_selected),
    Genes_Lambda1SE = nrow(lasso_1se)
  ),
  file.path(lasso_dir, "LASSO_Summary.csv"),
  row.names = FALSE
)

cat("\nLASSO completed.\n")
cat("Lambda.min:", cv_lasso$lambda.min, "\n")
cat("Genes selected at lambda.min:", nrow(lasso_selected), "\n")
cat("Lambda.1se:", cv_lasso$lambda.1se, "\n")
cat("Genes selected at lambda.1se:", nrow(lasso_1se), "\n")
cat("Results saved in:", lasso_dir, "\n")

cat("\nSelected genes at lambda.min:\n")
print(lasso_selected)


# ============================================================
# STEP 9: RANDOM FOREST FEATURE SELECTION
# Candidate genes: LASSO lambda.min
# ============================================================

if (!requireNamespace("randomForest", quietly = TRUE)) {
  stop("Install the randomForest package before continuing.")
}

library(randomForest)

rf_dir <- file.path(
  output_dir,
  "05_Machine_Learning",
  "02_Random_Forest"
)

dir.create(
  rf_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 1. Use the 41 LASSO-selected genes
rf_genes <- lasso_selected$Gene_Symbol

stopifnot(
  length(rf_genes) > 1,
  all(rf_genes %in% rownames(expr_combat)),
  identical(
    colnames(expr_combat),
    as.character(sample_info_pca$Sample)
  )
)

# Samples x genes
rf_x <- as.data.frame(
  t(expr_combat[rf_genes, , drop = FALSE])
)

# Outcome factor
rf_y <- factor(
  sample_info_pca$Group,
  levels = c("Control", "AUD")
)

rf_data <- data.frame(
  Group = rf_y,
  rf_x,
  check.names = FALSE
)

stopifnot(
  nrow(rf_data) == 99,
  nlevels(rf_data$Group) == 2
)

# 2. Fit Random Forest with 500 trees
set.seed(123)

rf_model <- randomForest(
  Group ~ .,
  data = rf_data,
  ntree = 500,
  importance = TRUE,
  proximity = FALSE
)

# 3. Save variable importance
importance_matrix <- importance(rf_model)

rf_importance <- data.frame(
  Gene_Symbol = rownames(importance_matrix),
  importance_matrix,
  row.names = NULL,
  check.names = FALSE
)

# Rank by MeanDecreaseGini, when available
if ("MeanDecreaseGini" %in% names(rf_importance)) {
  rf_importance <- rf_importance[
    order(rf_importance$MeanDecreaseGini, decreasing = TRUE),
  ]
}

write.csv(
  rf_importance,
  file.path(rf_dir, "Random_Forest_Gene_Importance.csv"),
  row.names = FALSE
)

# 4. Plot top 20 genes
importance_col <- if (
  "MeanDecreaseGini" %in% names(rf_importance)
) {
  "MeanDecreaseGini"
} else {
  colnames(importance_matrix)[1]
}

top20_rf <- head(rf_importance, 20)

p_rf <- ggplot(
  top20_rf,
  aes(
    x = reorder(Gene_Symbol, .data[[importance_col]]),
    y = .data[[importance_col]]
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Random Forest: Top 20 Candidate Genes",
    subtitle = "Input genes selected by LASSO at lambda.min",
    x = "Gene",
    y = importance_col
  ) +
  theme_bw(base_size = 12)

print(p_rf)

ggsave(
  file.path(rf_dir, "Random_Forest_Top20_Importance.png"),
  plot = p_rf,
  width = 9,
  height = 7,
  dpi = 300
)

# 5. Print model summary
cat("\nRandom Forest completed.\n")
cat("Trees:", rf_model$ntree, "\n")
cat("Candidate genes:", length(rf_genes), "\n")
cat("Out-of-bag error rate:",
    round(100 * rf_model$err.rate[rf_model$ntree, "OOB"], 2),
    "%\n")

cat("\nTop 20 genes by importance:\n")
print(top20_rf)

cat("\nResults saved in:", rf_dir, "\n")

# ============================================================
# RANDOM FOREST ON ALL 12,373 COMMON GENES
# ============================================================

library(randomForest)

rf_full_dir <- file.path(
  output_dir,
  "05_Machine_Learning",
  "03_Random_Forest_All_Genes"
)

dir.create(rf_full_dir, recursive = TRUE, showWarnings = FALSE)

# Samples x genes
rf_full_x <- as.data.frame(t(expr_combat), check.names = FALSE)

# Outcome
rf_full_y <- factor(
  sample_info_pca$Group,
  levels = c("Control", "AUD")
)

stopifnot(
  nrow(rf_full_x) == 99,
  ncol(rf_full_x) == 12373,
  identical(rownames(rf_full_x),
            as.character(sample_info_pca$Sample)),
  all(is.finite(as.matrix(rf_full_x)))
)

# Fit Random Forest
set.seed(123)

rf_full_model <- randomForest(
  x = rf_full_x,
  y = rf_full_y,
  ntree = 500,
  importance = TRUE
)

# Extract and rank gene importance
rf_full_imp <- as.data.frame(
  importance(rf_full_model)
)

rf_full_imp$Gene_Symbol <- rownames(rf_full_imp)

rf_full_imp <- rf_full_imp[
  order(rf_full_imp$MeanDecreaseGini, decreasing = TRUE),
]

rownames(rf_full_imp) <- NULL

# Save complete importance table
write.csv(
  rf_full_imp,
  file.path(rf_full_dir, "Random_Forest_All_12373_Genes.csv"),
  row.names = FALSE
)

# Select top 20 genes
top20_full_rf <- head(rf_full_imp, 20)

write.csv(
  top20_full_rf,
  file.path(rf_full_dir, "Random_Forest_Top20_All_Genes.csv"),
  row.names = FALSE
)

# Plot top 20
p_rf_full <- ggplot(
  top20_full_rf,
  aes(
    x = reorder(Gene_Symbol, MeanDecreaseGini),
    y = MeanDecreaseGini
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Random Forest: Top 20 Genes",
    subtitle = "Feature selection from 12,373 common genes",
    x = "Gene",
    y = "Mean Decrease in Gini"
  ) +
  theme_bw(base_size = 12)

print(p_rf_full)

ggsave(
  file.path(rf_full_dir, "Random_Forest_Top20_All_Genes.png"),
  p_rf_full,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Genes ranked:", nrow(rf_full_imp), "\n")
cat("Out-of-bag error (%):",
    round(100 * rf_full_model$err.rate[
      rf_full_model$ntree, "OOB"
    ], 2), "\n")

cat("\nTop 20 genes:\n")
print(top20_full_rf)

cat("\nResults saved in:", rf_full_dir, "\n")

# ============================================================
# STEP 10: INTERSECT LASSO AND RANDOM FOREST
# ============================================================

# LASSO genes selected at lambda.min
lasso_genes <- unique(lasso_selected$Gene_Symbol)

# Top 20 genes from Random Forest trained on all genes
rf_top20_genes <- head(rf_full_imp$Gene_Symbol, 20)

# Find overlapping genes
overlap_genes <- intersect(
  lasso_genes,
  rf_top20_genes
)

cat("LASSO genes at lambda.min:", length(lasso_genes), "\n")
cat("Random Forest top-20 genes:", length(rf_top20_genes), "\n")
cat("Overlapping candidate genes:", length(overlap_genes), "\n\n")

print(overlap_genes)

# Save overlap results in a separate folder
overlap_dir <- file.path(
  output_dir,
  "05_Machine_Learning",
  "04_LASSO_RF_Intersection"
)

dir.create(
  overlap_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Combine importance and LASSO coefficients
overlap_table <- merge(
  lasso_selected[, c("Gene_Symbol", "Coefficient")],
  rf_full_imp[, c("Gene_Symbol", "MeanDecreaseGini")],
  by = "Gene_Symbol"
)

overlap_table <- overlap_table[
  overlap_table$Gene_Symbol %in% overlap_genes,
]

overlap_table <- overlap_table[
  order(overlap_table$MeanDecreaseGini, decreasing = TRUE),
]

write.csv(
  overlap_table,
  file.path(overlap_dir, "LASSO_RF_Overlapping_Candidates.csv"),
  row.names = FALSE
)

cat("\nOverlap table saved in:", overlap_dir, "\n")

# ============================================================
# INTERSECTION: 41 LASSO GENES vs ALL 12,373 RF-RANKED GENES
# ============================================================

# All 41 genes selected by LASSO at lambda.min
lasso_genes <- unique(lasso_selected$Gene_Symbol)

# All genes ranked by Random Forest
rf_all_genes <- unique(rf_full_imp$Gene_Symbol)

# Intersect the complete lists
overlap_all_rf <- intersect(
  lasso_genes,
  rf_all_genes
)

cat("LASSO genes:", length(lasso_genes), "\n")
cat("Genes ranked by Random Forest:", length(rf_all_genes), "\n")
cat("Genes common to both methods:", length(overlap_all_rf), "\n\n")

# Extract RF importance and LASSO coefficient for overlapping genes
overlap_table_all <- merge(
  lasso_selected[, c("Gene_Symbol", "Coefficient")],
  rf_full_imp[, c("Gene_Symbol", "MeanDecreaseGini")],
  by = "Gene_Symbol"
)

# Sort by RF importance
overlap_table_all <- overlap_table_all[
  order(overlap_table_all$MeanDecreaseGini, decreasing = TRUE),
]

print(overlap_table_all)

# Save complete intersection in a separate folder
overlap_dir <- file.path(
  output_dir,
  "05_Machine_Learning",
  "04_LASSO_RF_Intersection"
)

dir.create(
  overlap_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  overlap_table_all,
  file.path(
    overlap_dir,
    "LASSO_41_vs_RandomForest_AllGenes_Intersection.csv"
  ),
  row.names = FALSE
)

cat("\nFull intersection saved in:", overlap_dir, "\n")

# ============================================================
# STEP 11: INTEGRATE LASSO, RANDOM FOREST AND DGE RESULTS
# ============================================================

# Start with the 41 LASSO-selected genes
candidate_table <- merge(
  lasso_selected[, c("Gene_Symbol", "Coefficient")],
  rf_full_imp[, c("Gene_Symbol", "MeanDecreaseGini")],
  by = "Gene_Symbol",
  all.x = TRUE
)

# Add differential-expression statistics
dge_for_merge <- dge_results[
  ,
  c("Gene_Symbol", "logFC", "P.Value", "adj.P.Val")
]

candidate_table <- merge(
  candidate_table,
  dge_for_merge,
  by = "Gene_Symbol",
  all.x = TRUE
)

# Add expression direction and significance labels
candidate_table$DGE_Direction <- ifelse(
  candidate_table$logFC > 0,
  "Higher in AUD",
  ifelse(candidate_table$logFC < 0, "Lower in AUD", "No change")
)

candidate_table$Nominal_Significant <- (
  candidate_table$P.Value < 0.05
)

candidate_table$FDR_Significant <- (
  candidate_table$adj.P.Val < 0.05
)

# Rank by RF importance
candidate_table <- candidate_table[
  order(candidate_table$MeanDecreaseGini, decreasing = TRUE),
]

# Save integrated candidate table
candidate_dir <- file.path(
  output_dir,
  "05_Machine_Learning",
  "05_Integrated_Candidates"
)

dir.create(
  candidate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  candidate_table,
  file.path(candidate_dir, "LASSO_RF_DGE_Integrated_Candidates.csv"),
  row.names = FALSE
)

cat("Total integrated candidates:", nrow(candidate_table), "\n")
cat("Nominally significant DGE genes:",
    sum(candidate_table$Nominal_Significant, na.rm = TRUE), "\n")
cat("FDR-significant DGE genes:",
    sum(candidate_table$FDR_Significant, na.rm = TRUE), "\n\n")

print(candidate_table)

# Read the complete differential-expression results
all_genes <- read.csv(
  "/Users/abarnageethanjali/Documents/set/combine datasets/affymetrix datasets /AUD_Fresh_Analysis_Results/04_Differential_Gene_Expression/All_12373_Genes.csv",
  check.names = FALSE
)

# 1. FDR-significant genes
fdr_genes <- subset(
  all_genes,
  adj.P.Val < 0.05
)

# 2. Nominally significant genes
nominal_genes <- subset(
  all_genes,
  P.Value < 0.05
)

# 3. Nominal significance + fold-change threshold
paper_genes <- subset(
  all_genes,
  P.Value < 0.05 & abs(logFC) > 1
)

# Sort each table by significance
fdr_genes <- fdr_genes[
  order(fdr_genes$adj.P.Val), ]

nominal_genes <- nominal_genes[
  order(nominal_genes$P.Value), ]

paper_genes <- paper_genes[
  order(paper_genes$P.Value), ]

# Save the results
write.csv(
  fdr_genes,
  "FDR_Significant_Genes.csv",
  row.names = FALSE
)

write.csv(
  nominal_genes,
  "Nominally_Significant_Genes.csv",
  row.names = FALSE
)

write.csv(
  paper_genes,
  "Nominal_DEGs_PValue_FC.csv",
  row.names = FALSE
)

# Print gene counts
cat("Total genes:", nrow(all_genes), "\n")
cat("FDR-significant:", nrow(fdr_genes), "\n")
cat("Nominally significant:", nrow(nominal_genes), "\n")
cat("Nominal DEGs with |logFC| > 1:",
    nrow(paper_genes), "\n")



