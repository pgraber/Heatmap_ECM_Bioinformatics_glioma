# Set working directory
setwd(here())

# Load libraries
library(readxl)
library(gplots)
library(ComplexHeatmap)
library(fastcluster)
library(circlize)
library(here)

# Define file paths
excel_file <- here("data/heatmap/ZERO_ECMmarkers BT_patient_metadata.xlsx")
counts_file <- here("data/heatmap/Zero_ECM_Markers_BT_counts.txt")

# Read sheet names
sheet_nam <- excel_sheets(excel_file)

# Read Gliomas data
data_ECM_Gliomas <- read_excel(excel_file, sheet = "Gliomas")
data_ECM_Gliomas <- data_ECM_Gliomas[-c(nrow(data_ECM_Gliomas), nrow(data_ECM_Gliomas) - 1),]
metadata <- data_ECM_Gliomas[, 1:4]

# Extract TPM data
data_tpm <- data_ECM_Gliomas[, c(1, 5:ncol(data_ECM_Gliomas))]
rownames(data_tpm) <- data_tpm$Patient.ID

# Read average data
data_ECM_avg <- read_excel(excel_file, sheet = "Gliomas Rank")
data_tpm_avg <- as.data.frame(t(data_tpm[, -1]))
colnames(data_tpm_avg) <- data_tpm$Patient.ID
data_tpm_avg$avg_exp <- rowMeans(data_tpm_avg)

# Order data based on average expression
data_tpm_avg_order <- data_tpm_avg[order(data_tpm_avg$avg_exp, decreasing = TRUE),]
x <- nrow(data_tpm_avg_order) - 29
data_tpm_avg_order_top <- data_tpm_avg_order[c(1:30, x:nrow(data_tpm_avg_order)),]

# Extract relevant columns for TPM based on ordering
index <- which(colnames(data_tpm) %in% rownames(data_tpm_avg_order_top))
data_tpm_top <- data_tpm[, c(1, index)]

# Transpose data
data_tpm_t <- t(data_tpm_top[, -1])
colnames(data_tpm_t) <- data_tpm$Patient.ID

# Handle zero values
data_tpm_t[data_tpm_t == 0] <- 0.001

# Create log-transformed TPM matrix
logTPM <- log(as.matrix(data_tpm_t))

# Scale the matrix before using Complexheatmap
log_scale <- scale(logTPM)

# Define color ramps and clustering function
col_fun1 <- colorRamp2(c(-1, 0, 1), c("darkorange", "white", "deeppink"))
col_fun2 <- colorRamp2(c(1, 5, 10), c("seagreen1", "royalblue", "orangered"))
fh <- function(x) fastcluster::hclust(dist(x))

# Define SubtypeCol colors
SubtypeCol <- c("Glioma other" = "goldenrod", "HGG" = "darkgreen", "DMG" = "lightskyblue1")

# Match column names for metadata
m <- match(colnames(log_scale), metadata$Patient.ID)
metadata_order <- metadata[m,]
metadata_order_1 <- as.data.frame(metadata_order[, 2])

# Define HeatmapAnnotation
col_at <- HeatmapAnnotation(df = metadata_order_1, col = list(Diagnosis = SubtypeCol))

# Set font size
font_size <- 8

# Create PDF output for clustered heatmap
pdf("./output/heatmap/heatmap_gliomas_top_30_cluster_scaled.pdf", height = 11, width = 10)
Heatmap(as.matrix(log_scale), name = "z-scores", show_row_names = TRUE, show_column_names = FALSE,
        top_annotation = col_at, row_names_gp = gpar(fontsize = font_size), cluster_rows = fh,
        cluster_columns = fh)
dev.off()

# Save ordered data to file
write.table(data_tpm_avg_order, "./output/heatmap/ZERO_ECM_glioma_fil_genelist.txt", sep = "\t", quote = FALSE)

# Filtering low counts
data_counts <- read.delim(counts_file, sep = "\t", stringsAsFactors = FALSE)
library(edgeR)

# Keep rows with counts greater than threshold
keep <- rowSums(cpm(data_counts[, 3:ncol(data_counts)]) > 1) >= (0.5 * ncol(data_counts[, 3:ncol(data_counts)]) + 1)
counts <- data_counts[keep, ]
counts_left <- data_counts[!keep, ]

# Read Gliomas data again
data_ECM_Gliomas <- read_excel(excel_file, sheet = "Gliomas")
data_ECM_Gliomas <- data_ECM_Gliomas[-c(nrow(data_ECM_Gliomas), nrow(data_ECM_Gliomas) - 1), ]

# Extract relevant columns based on filtered counts
index <- which(colnames(data_ECM_Gliomas) %in% counts$gene_id)
data_ECM_Gliomas.1 <- data_ECM_Gliomas[, c(1:4, index)]

# Write filtered gene list to a file
data_tpm_out <- data_ECM_Gliomas[, -index]
writeLines(colnames(data_tpm_out)[-c(1:4)], "./output/heatmap/filtered_geneList.txt")

metadata <- data_ECM_Gliomas[, 1:4]

# Extract TPM data again
data_tpm <- data_ECM_Gliomas[, c(1, 5:ncol(data_ECM_Gliomas))]
rownames(data_tpm) <- data_tpm$Patient.ID

# Process average data again
data_tpm_avg <- data_tpm
data_tpm_avg <- data_tpm_avg[, -1]
data_tpm_avg <- as.data.frame(t(data_tpm_avg))
colnames(data_tpm_avg) <- data_tpm$Patient.ID
data_tpm_avg$avg_exp <- rowMeans(data_tpm_avg)

data_tpm_avg_order <- data_tpm_avg[order(data_tpm_avg$avg_exp, decreasing = TRUE), ]
x <- nrow(data_tpm_avg_order) - 29
data_tpm_avg_order_top <- data_tpm_avg_order[c(1:30, x:nrow(data_tpm_avg_order)), ]

# Extract relevant columns for TPM based on ordering
index <- which(colnames(data_tpm) %in% rownames(data_tpm_avg_order_top))
data_tpm_top <- data_tpm[, c(1, index)]

data_tpm_t <- t(data_tpm_top[, -1])
colnames(data_tpm_t) <- data_tpm$Patient.ID

mat <- data_tpm_t
mat[mat == 0] <- 0.001
tpm <- as.matrix(mat)

logTPM <- log(tpm)

# Complex heatmap again

## Scale the matrix before using Complexheatmap
log_scale <- scale(logTPM)

pdf("./output/heatmap/heatmap_gliomas_top_30_cluster_filtered_scaled.pdf", height = 11, width = 10)
Heatmap(as.matrix(log_scale), name = "z-scores", show_row_names = TRUE, show_column_names = FALSE,
        top_annotation = col_at, row_names_gp = gpar(fontsize = font_size), cluster_rows = fh,
        cluster_columns = fh)
dev.off()

# Save ordered data to file
write.table(data_tpm_avg_order, "./output/heatmap/ZERO_ECM_glioma_fil_genelist.txt", sep = "\t", quote = FALSE)



