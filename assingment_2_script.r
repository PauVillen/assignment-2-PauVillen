# Assignment 2
setwd("/Users/pauvillen14/Desktop/BIOINFO/DMI/assignment-2-PauVillen")

# Import all the necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
library(cluster)   
library(factoextra) 
library(ggpubr)
library(dendextend)

# -------------------------------------------------------------------------
# EXERCISE 1 - K-means clustering
# -------------------------------------------------------------------------

# --- PREPARE DATA (Table s2) ---
# Read table setting correct headers
patient_id_row <- read_excel("table_s2.xlsx", sheet = "Proteomics_proteins_training", 
                             skip = 1, n_max = 1, col_names = FALSE)
patient_names <- as.character(unlist(patient_id_row))
df2_raw <- read_excel("table_s2.xlsx", sheet = "Proteomics_proteins_training", 
                     skip = 2, col_names = FALSE)
colnames(df2_raw) <- c("Proteins", "Gene_Symbol", patient_names)

# Remove patients identified as outliers 
outliers_to_remove <- c("F2_129N", "F2_130N", "F1_130C", "F4_130C", "F1_130N", 
                        "F2_130C", "F2_127C")

# Keep all the data frame except the outliers
cols_to_keep <- setdiff(colnames(df2_raw), outliers_to_remove)
df2_raw <- df2_raw[, cols_to_keep]

# Filter missing values (delete proteins with >= 20% missing values)
df2_numeric <- df2_raw %>% 
  select(-Gene_Symbol) %>% 
  column_to_rownames("Proteins")
df2_numeric <- as.data.frame(lapply(df2_numeric, as.numeric))
rownames(df2_numeric) <- df2_raw$Proteins
missing_rate <- rowMeans(is.na(df2_numeric))
df2_filtered <- df2_numeric[missing_rate < 0.20, ]

# Impute missing values as in the paper: Replace all NA with 0
df2_imputed <- df2_filtered
df2_imputed[is.na(df2_imputed)] <- 0

# Transpose the data frame (we want Patient IDs as rows and Proteins as columns)
df2_transposed <- t(df2_imputed)

# Scale variables before clustering to have mean 0 and sd 1 so that they are comparable
df2_scaled <- scale(df2_transposed)

# --- K-MEANS CLUSTERING ---
# Choosing the optimal number of clusters with Elbow Method
fviz_nbclust(df2_scaled, kmeans, method = "wss") +
  labs(title = "Elbow Method (WSS)")

# Choosing the optimal number of clusters with Silhouette Method
fviz_nbclust(df2_scaled, kmeans, method = "silhouette") +
  labs(title = "Silhouette Method")

# Set the optimal k = 2 based on previous methods and perform K-means clustering
set.seed(123)
k_optimal <- 2
km <- kmeans(df2_scaled, centers = k_optimal, nstart = 25)

# Visualize the clusters
fviz_cluster(km, data = df2_scaled,
             geom = "text",             
             ellipse.type = "convex",    
             ggtheme = theme_minimal(),
             main = "K-Means Clustering") + 
    guides(color = guide_legend(override.aes = list(label = "")))


# --- PREPARE DATA (Table s1) ---
# Load data: treating "/", "NA", and empty strings as missing values
df1 <- read_excel("table_s1.xlsx", sheet = 2, na = c("", "NA", "/"))

# Remove almost empty column
df1 <- df1 %>% select(-`MSRep ID c`)

# Clean column names (replace spaces with underscores)
names(df1) <- gsub(" ", "_", names(df1))
names(df1) <- gsub("\\(", "", names(df1))
names(df1) <- gsub("\\)", "", names(df1))

# Ensure numeric columns
df1$BMI_h <- as.numeric(df1$BMI_h)
df1$`CRP_i,_mg/L` <- as.numeric(df1$`CRP_i,_mg/L`)

# Date conversion & time difference calculation
df1 <- df1 %>%
  mutate(
    # Convert columns to Date format (assuming YYYY-MM-DD)
    Onset_date_f = as.Date(Onset_date_f),
    Admission_date = as.Date(Admission_date),
    Date_of_progression_to_severe_state = as.Date(Date_of_progression_to_severe_state),
    
    # Calculate differences in days
    time_onset_to_admission = as.numeric(Admission_date - Onset_date_f),
    time_admission_to_severe = as.numeric(Date_of_progression_to_severe_state - Admission_date)
  )


# Create group labels
df1 <- df1 %>%
  mutate(
    Group = case_when(
      Group_d == 0 ~ "healthy",
      Group_d == 1 ~ "non_covid",
      Group_d == 2 ~ "non_severe",
      Group_d == 3 ~ "severe"
    )
  )

# Create Sex labels (Factor)
df1 <- df1 %>%
  mutate(Sex_g = factor(Sex_g, 
                        levels = c(0, 1),            
                        labels = c("Female", "Male"))) 

# --- VISUALIZE PREVIOUS CLUSTERS WITH CLINICAL INFO ---
df1$name <- as.character(df1$`MS_ID_b`)
p <- fviz_cluster(km, data = df2_scaled, geom = "point", ellipse.type = "convex")
plot_data <- p$data

# Match IDs
final_plot_data <- left_join(plot_data, df1, by = "name")

# Check for missing clinical info 
missing_clinical <- final_plot_data %>% filter(is.na(Group))
if(nrow(missing_clinical) > 0) {
  message("Warning: The following patients have no Clinical Group info:")
  print(missing_clinical$name)
}

# Set clean names
final_plot_data$Group <- factor(final_plot_data$Group,
                                levels = c("healthy", "non_covid", "non_severe", "severe"),
                                labels = c("Healthy", "Non-COVID-19", "COVID-19 Non-severe", "COVID-19 Severe"))

# Plot
ggplot(final_plot_data, aes(x = x, y = y)) +
  
  # Draw the Cluster Areas 
  stat_chull(aes(group = cluster), 
             geom = "polygon", 
             fill = "grey",   
             color = "black", 
             alpha = 0.1) +   
  
  # Draw the Points (Colored by Clinical Group)
  geom_point(aes(color = Group, shape = cluster), size = 3) +
  
  # Formatting
  scale_color_brewer(palette = "Set1", na.value = "grey50") + 

  scale_shape_manual(values = c(16, 17)) +
  
  labs(title = "K-Means Clusters Colored by Clinical Group",
       subtitle = "Shapes: Cluster | Colors: Clinical Group",
       x = "Dim1 (PCA)", y = "Dim2 (PCA)",
       color = "Clinical Group", 
       shape = "Cluster") +      
  
  guides(fill = "none") +
  
  theme_minimal()


# -------------------------------------------------------------------------
# EXERCISE 1 - Hierarchical clustering
# -------------------------------------------------------------------------

# Dissimilarity matrix
group_colors <- c(
  "healthy"    = "blue",  # Blue
  "non_covid"  = "green",  # Teal
  "non_severe" = "yellow",  # Yellow
  "severe"     = "red"   # Red
)
df1$color <- group_colors[as.character(df1$Group)]
d_mat <- dist(df2_scaled, method = "euclidean")

# Silhouette Method
fviz_nbclust(df2_scaled, FUN = hcut, method = "silhouette") +
  labs(title = "Silhouette Method")

# Perform the clustering
hcl <- hclust(d_mat, method = "complete")

# Prepare the plot including clinical information
dend <- as.dendrogram(hcl)
dend_labels <- labels(dend)
ordered_colors <- df1$color[match(dend_labels, df1$MS_ID_b)]
labels_colors(dend) <- ordered_colors
plot(dend, 
     main = "Dendrogram Colored by Clinical Group",
     ylab = "Height (Euclidean Distance)",
     cex = 0.6)

legend("topright", 
       legend = names(group_colors), 
       fill = group_colors, 
       bty = "n", 
       cex = 0.8)

rect.hclust(hcl, k = 2, border = "red")
hc_clusters <- cutree(hcl, k = 2)
print(table(hc_clusters))



