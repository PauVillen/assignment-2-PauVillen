setwd("/Users/pauvillen14/Desktop/BIOINFO/DMI/assignment-2-PauVillen")

required_packages <- c("readxl", "Rtsne", "umap", "dplyr", "RColorBrewer", "ggplot2", "BiocManager")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(missing_packages) > 0) install.packages(missing_packages)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE)) BiocManager::install("circlize")

library(readxl)
library(Rtsne)
library(umap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Data loading and preprocessing
raw <- read_excel("table_s2.xlsx", sheet = "Proteomics_proteins_training", col_names = FALSE)

ms_ids <- as.character(unlist(raw[2, 3:ncol(raw)]))
ms_ids <- ms_ids[!is.na(ms_ids)]
n_samples <- length(ms_ids)

protein_ids_gene <- as.character(unlist(raw[3:nrow(raw), 2]))
protein_ids_uniprot <- as.character(unlist(raw[3:nrow(raw), 1]))
protein_ids <- ifelse(is.na(protein_ids_gene) | protein_ids_gene == "NA", protein_ids_uniprot, protein_ids_gene)

protein_matrix <- raw[3:nrow(raw), 3:(2 + n_samples)]
protein_matrix <- as.data.frame(lapply(protein_matrix, as.numeric))

protein_matrix_t <- t(protein_matrix)
rownames(protein_matrix_t) <- ms_ids
colnames(protein_matrix_t) <- protein_ids

clinical <- read_excel("table_s1.xlsx", sheet = "Clinical_information")
severity_map <- c("0" = "Healthy", "1" = "Non-COVID-19", "2" = "Non-severe", "3" = "Severe")
clinical$Severity <- severity_map[as.character(clinical$`Group d`)]

ms_to_severity <- setNames(clinical$Severity, clinical$`MS ID b`)
ms_to_patient <- setNames(clinical$`Patient ID a`, clinical$`MS ID b`)

replicate_ids <- c("F1_133N", "F3_133N", "F5_133N", "F6_133N")
keep_samples <- rownames(protein_matrix_t)[!rownames(protein_matrix_t) %in% replicate_ids]
keep_samples <- keep_samples[keep_samples %in% names(ms_to_severity)]
keep_samples <- keep_samples[!is.na(ms_to_severity[keep_samples])]

protein_matrix_t <- protein_matrix_t[keep_samples, ]

col_means <- colMeans(protein_matrix_t, na.rm = TRUE)
for (i in 1:ncol(protein_matrix_t)) {
  nas <- is.na(protein_matrix_t[, i])
  if (any(nas)) protein_matrix_t[nas, i] <- col_means[i]
}

col_var <- apply(protein_matrix_t, 2, var, na.rm = TRUE)
col_sd <- apply(protein_matrix_t, 2, sd, na.rm = TRUE)
protein_matrix_clean <- protein_matrix_t[, col_var > 1e-10 & col_sd > 1e-10 & !is.na(col_var)]
protein_matrix_clean <- protein_matrix_clean[, apply(protein_matrix_clean, 2, function(x) length(unique(x)) > 1)]
protein_matrix_clean <- protein_matrix_clean[, grepl("^[A-Z][A-Z0-9]", colnames(protein_matrix_clean))]

sample_severity <- ms_to_severity[rownames(protein_matrix_clean)]
severity_colors <- c("Healthy" = "#2196F3", "Non-COVID-19" = "#9C27B0", "Non-severe" = "#4CAF50", "Severe" = "#F44336")

# PCA
pca_result <- prcomp(protein_matrix_clean, scale. = TRUE, center = TRUE)
var_explained <- summary(pca_result)$importance[2, ] * 100

pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x),
  Severity = sample_severity[rownames(pca_result$x)]
)
pca_data <- pca_data[!is.na(pca_data$Severity), ]

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Severity)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = severity_colors) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#fafafa", color = NA),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "PCA of Proteomic Data",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Disease Severity"
  )
print(p_pca)

loadings <- as.data.frame(pca_result$rotation)

pc1_order <- order(abs(loadings[, 1]), decreasing = TRUE)[1:20]
top20_pc1 <- data.frame(Protein = rownames(loadings)[pc1_order], Loading = loadings[pc1_order, 1])
top20_pc1 <- top20_pc1[order(top20_pc1$Loading), ]
top20_pc1$Protein <- factor(top20_pc1$Protein, levels = top20_pc1$Protein)
top20_pc1$Direction <- ifelse(top20_pc1$Loading > 0, "Positive", "Negative")

p_pc1 <- ggplot(top20_pc1, aes(x = Protein, y = Loading, fill = Direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#E74C3C", "Negative" = "#3498DB"), guide = "none") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9, color = "#333333"),
    plot.title = element_text(face = "bold", size = 13),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#fafafa", color = NA)
  ) +
  labs(title = "Top 20 Proteins Contributing to PC1", x = NULL, y = "PC1 Loading")
print(p_pc1)
ggsave("pc1_loadings.png", p_pc1, width = 8, height = 6, dpi = 150)

pc2_order <- order(abs(loadings[, 2]), decreasing = TRUE)[1:20]
top20_pc2 <- data.frame(Protein = rownames(loadings)[pc2_order], Loading = loadings[pc2_order, 2])
top20_pc2 <- top20_pc2[order(top20_pc2$Loading), ]
top20_pc2$Protein <- factor(top20_pc2$Protein, levels = top20_pc2$Protein)
top20_pc2$Direction <- ifelse(top20_pc2$Loading > 0, "Positive", "Negative")

p_pc2 <- ggplot(top20_pc2, aes(x = Protein, y = Loading, fill = Direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#E74C3C", "Negative" = "#3498DB"), guide = "none") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9, color = "#333333"),
    plot.title = element_text(face = "bold", size = 13),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#fafafa", color = NA)
  ) +
  labs(title = "Top 20 Proteins Contributing to PC2", x = NULL, y = "PC2 Loading")
print(p_pc2)
ggsave("pc2_loadings.png", p_pc2, width = 8, height = 6, dpi = 150)

# t-SNE
set.seed(42)
n_samples_tsne <- nrow(protein_matrix_clean)
max_perp <- floor((n_samples_tsne - 1) / 3) - 1
perp_high <- min(20, max_perp)
perp_low <- min(10, max_perp)

tsne_high <- Rtsne(protein_matrix_clean, perplexity = perp_high, verbose = FALSE, max_iter = 1000, pca = FALSE)
tsne_low <- Rtsne(protein_matrix_clean, perplexity = perp_low, verbose = FALSE, max_iter = 1000, pca = FALSE)

tsne_df_high <- data.frame(tSNE1 = tsne_high$Y[, 1], tSNE2 = tsne_high$Y[, 2], Sample = rownames(protein_matrix_clean))
tsne_df_high$Severity <- sample_severity[tsne_df_high$Sample]
tsne_df_high <- tsne_df_high[!is.na(tsne_df_high$Severity), ]

tsne_df_low <- data.frame(tSNE1 = tsne_low$Y[, 1], tSNE2 = tsne_low$Y[, 2], Sample = rownames(protein_matrix_clean))
tsne_df_low$Severity <- sample_severity[tsne_df_low$Sample]
tsne_df_low <- tsne_df_low[!is.na(tsne_df_low$Severity), ]

make_tsne_plot <- function(df, perp_val) {
  ggplot(df, aes(x = tSNE1, y = tSNE2, color = Severity)) +
    geom_point(size = 4, alpha = 0.85) +
    scale_color_manual(values = severity_colors) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#fafafa", color = NA),
      legend.title = element_text(face = "bold")
    ) +
    labs(
      title = paste0("t-SNE (perplexity = ", perp_val, ")"),
      x = "t-SNE 1", y = "t-SNE 2", color = "Disease Severity"
    )
}

p_tsne_high <- make_tsne_plot(tsne_df_high, perp_high)
p_tsne_low <- make_tsne_plot(tsne_df_low, perp_low)
print(p_tsne_high)
print(p_tsne_low)

# UMAP
set.seed(42)
umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1
umap_result <- umap(protein_matrix_clean, config = umap_config)

umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], Sample = rownames(protein_matrix_clean))
umap_df$Severity <- sample_severity[umap_df$Sample]
umap_df <- umap_df[!is.na(umap_df$Severity), ]

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Severity)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = severity_colors) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#fafafa", color = NA),
    legend.title = element_text(face = "bold")
  ) +
  labs(title = "UMAP of Proteomic Data", x = "UMAP 1", y = "UMAP 2", color = "Disease Severity")
print(p_umap)

# Heatmap
table_s6 <- read_excel("1-s2.0-S0092867420306279-mmc6.xlsx", sheet = "Prot Severe vs Healthy", skip = 1)
colnames(table_s6) <- c("Proteins", "Gene_Symbol", "fd", "P_value", "P_value_adj")

dep_proteins <- as.character(table_s6$Proteins)
dep_proteins <- dep_proteins[!is.na(dep_proteins)]
gene_symbols <- as.character(table_s6$Gene_Symbol)
matrix_cols <- colnames(protein_matrix_clean)

normalize_name <- function(x) toupper(trimws(x))

dep_proteins_in_matrix <- c()
for (i in seq_along(dep_proteins)) {
  genes <- normalize_name(trimws(strsplit(gene_symbols[i], ",")[[1]]))
  matched <- matrix_cols[sapply(matrix_cols, function(col) {
    col_genes <- normalize_name(trimws(strsplit(col, ",")[[1]]))
    any(genes %in% col_genes)
  })]
  if (length(matched) > 0) dep_proteins_in_matrix <- c(dep_proteins_in_matrix, matched[1])
}
dep_proteins_in_matrix <- unique(dep_proteins_in_matrix)

if (length(dep_proteins_in_matrix) > 0) {
  heatmap_matrix <- t(protein_matrix_clean[, dep_proteins_in_matrix, drop = FALSE])
  
  severe_samples <- names(sample_severity)[sample_severity == "Severe"]
  healthy_samples <- names(sample_severity)[sample_severity == "Healthy"]
  patients_for_heatmap <- c(severe_samples, healthy_samples)
  patients_available <- intersect(patients_for_heatmap, colnames(heatmap_matrix))
  heatmap_matrix <- heatmap_matrix[, patients_available, drop = FALSE]
  
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
  row_var_hm <- apply(heatmap_matrix_scaled, 1, var, na.rm = TRUE)
  heatmap_matrix_scaled <- heatmap_matrix_scaled[!is.na(row_var_hm) & row_var_hm > 1e-10, ]
  
  severity_col <- sample_severity[patients_available]
  sex_col <- as.character(clinical$`Sex g`[match(ms_to_patient[patients_available], clinical$`Patient ID a`)])
  age_col <- clinical$`Age (year)`[match(ms_to_patient[patients_available], clinical$`Patient ID a`)]
  age_group <- ifelse(is.na(age_col), NA, ifelse(age_col < 40, "<40", ifelse(age_col <= 60, "40-60", ">60")))
  
  annotation_df <- data.frame(
    Group = severity_col,
    Gender = sex_col,
    Age_Group = age_group,
    row.names = patients_available
  )
  annotation_df <- annotation_df[!is.na(annotation_df$Group), ]
  heatmap_matrix_scaled <- heatmap_matrix_scaled[, rownames(annotation_df), drop = FALSE]
  
  ha <- HeatmapAnnotation(
    df = annotation_df,
    col = list(
      Group = c("Healthy" = "#2196F3", "Severe" = "#F44336"),
      Gender = c("0" = "#9C27B0", "1" = "#FF9800"),
      Age_Group = c("<40" = "#C8E6C9", "40-60" = "#388E3C", ">60" = "#1B5E20")
    ),
    annotation_name_side = "left",
    annotation_height = unit(3, "mm")
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  
  ht <- Heatmap(
    heatmap_matrix_scaled,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 7),
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    column_title = "Differentially Expressed Proteins: Severe vs Healthy",
    column_title_gp = gpar(fontsize = 13, fontface = "bold")
  )
  
  png("heatmap_dep.png", width = 2400, height = 3200, res = 200)
  draw(ht)
  dev.off()
}

# DEP overlap with top 50 PC1 and PC2
top50_pc1 <- rownames(loadings)[order(abs(loadings[, 1]), decreasing = TRUE)[1:50]]
top50_pc2 <- rownames(loadings)[order(abs(loadings[, 2]), decreasing = TRUE)[1:50]]

overlap_pc1 <- dep_proteins_in_matrix[sapply(dep_proteins_in_matrix, function(p) {
  any(normalize_name(trimws(strsplit(p, ",")[[1]])) %in% normalize_name(trimws(strsplit(paste(top50_pc1, collapse = ","), ",")[[1]])))
})]

overlap_pc2 <- dep_proteins_in_matrix[sapply(dep_proteins_in_matrix, function(p) {
  any(normalize_name(trimws(strsplit(p, ",")[[1]])) %in% normalize_name(trimws(strsplit(paste(top50_pc2, collapse = ","), ",")[[1]])))
})]

cat("DEPs in top 50 PC1 contributors:", length(overlap_pc1), "\n")
if (length(overlap_pc1) > 0) cat(overlap_pc1, sep = "\n")
cat("DEPs in top 50 PC2 contributors:", length(overlap_pc2), "\n")
if (length(overlap_pc2) > 0) cat(overlap_pc2, sep = "\n")