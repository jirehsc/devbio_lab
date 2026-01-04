# Combined figure: column dendrogram (top) + heatmap (middle) + PCA biplot (bottom)
# - PCA is computed on specimen-level data using correlation (i.e. z-scored traits).
# - Heatmap shows treatment means (rows = combined treatments, columns = traits).
# - Traits (columns) are clustered using correlation distance (1 - cor).
# - Treatments (rows) are clustered and treatment clusters are used to group specimens in the PCA
#   biplot (ellipses around specimen points based on their treatment cluster).
#
# Expected input: metadata2.csv in working directory (the dataset you provided).
#
# Required packages:
# install.packages(c("tidyverse","ggdendro","ggrepel","patchwork","scales","grid","ggalt"))
#
# Notes:
# - If a trait/column is constant across rows it will be removed from PCA and (if needed) handled for the heatmap.
# - You can change `k_clusters` to alter the number of treatment clusters used for PCA grouping.
# - The script tries to align the top column dendrogram with the heatmap columns.
#
# Usage:
# 1. Save this file and run in RStudio: source("combined_correlation_heatmap_pca.R")
# 2. Output file: combined_correlation_heatmap_pca.png
#
# Author: generated for user

library(tidyverse)
library(ggdendro)
library(ggrepel)
library(patchwork)
library(scales)
library(grid)
  # for geom_encircle if needed

set.seed(42)

# ---------------------------
# User parameters
# ---------------------------
csv_file <- "metadata2.csv"
group_cols <- c("ethanol_pre.treatment", "saline_treatment")  # used to make treatment labels
k_clusters <- 3            # number of treatment clusters to define (used to draw ellipses in PCA)
label_top_loadings <- 18   # how many variable loadings to label on PCA
output_file <- "combined_correlation_heatmap_pca.png"
width <- 14
height <- 12
dpi <- 300
# ---------------------------

# ---------------------------
# 1) Read data
# ---------------------------
df <- readr::read_csv(csv_file, show_col_types = FALSE)

# create combined treatment label (E{ethanol}_S{saline})
if (!all(group_cols %in% colnames(df))) {
  stop("Grouping columns not found: ", paste(setdiff(group_cols, colnames(df)), collapse = ", "))
}
df <- df %>% mutate(treatment = paste0("E", ethanol_pre.treatment, "_S", saline_treatment))

# ---------------------------
# 2) Select trait columns (numeric), exclude IDs / grouping columns
# ---------------------------
exclude_cols <- c("specimen_id", "species", "variety_name", "variety_type", "block", "health_status", group_cols, "treatment")
trait_df <- df %>%
  select(-any_of(exclude_cols)) %>%
  select(where(is.numeric)) %>%
  select(where(~ !all(is.na(.))))   # drop columns that are all NA

if (ncol(trait_df) < 2) stop("Need at least two numeric trait columns for PCA/heatmap. Check your data and exclude_cols.")

trait_cols <- colnames(trait_df)

# Remove specimens with all-NA traits (if any)
keep_rows <- rowSums(is.na(trait_df)) < ncol(trait_df)
if (any(!keep_rows)) {
  message("Removing ", sum(!keep_rows), " specimens with all NA trait values.")
  df <- df[keep_rows, ]
  trait_df <- trait_df[keep_rows, ]
}

# ---------------------------
# 3) Treatment means for heatmap
# ---------------------------
df_means <- df %>%
  group_by(treatment) %>%
  summarise(across(all_of(trait_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

mat_means <- df_means %>% select(all_of(trait_cols)) %>% as.matrix()
rownames(mat_means) <- df_means$treatment

# ---------------------------
# 4) Column clustering based on correlation (1 - cor)
#    This clusters traits (columns) by correlation across treatments
# ---------------------------
# If any column is constant across treatments, cor will be NA; drop constant columns for clustering but keep track
col_sds_means <- apply(mat_means, 2, sd, na.rm = TRUE)
const_cols_means <- names(col_sds_means)[is.na(col_sds_means) | col_sds_means <= 1e-12]
if (length(const_cols_means) > 0) {
  message("Constant trait columns across treatments (removed for trait clustering): ", paste(const_cols_means, collapse = ", "))
}
mat_for_col_clust <- mat_means[, setdiff(colnames(mat_means), const_cols_means), drop = FALSE]

# compute correlation among traits (columns)
cor_cols <- tryCatch(cor(mat_for_col_clust, use = "pairwise.complete.obs"), error = function(e) stop("Column correlation failed: ", e$message))
# convert to distance (1 - cor) and cluster
dist_cols <- as.dist(1 - cor_cols)
hc_cols <- hclust(dist_cols, method = "complete")
# column order including const columns appended at end (or at start): we'll place const_cols at end
col_order <- c(hc_cols$labels[hc_cols$order], const_cols_means)

# ---------------------------
# 5) Row clustering (treatments) based on treatment means (Euclidean)
# ---------------------------
# For rows, we can cluster using scaled values (z-score columns)
mat_means_scaled <- scale(mat_means, center = TRUE, scale = TRUE)
mat_means_scaled[is.na(mat_means_scaled)] <- 0
dist_rows <- dist(mat_means_scaled, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "complete")
row_order <- hc_rows$labels[hc_rows$order]

# Optionally cut treatment dendrogram into k clusters for grouping
treatment_cluster <- cutree(hc_rows, k = k_clusters)
treatment_cluster_df <- tibble(treatment = names(treatment_cluster), cluster = as.factor(treatment_cluster))

# ---------------------------
# 6) Prepare heatmap data (long) and plotting objects
# ---------------------------
heatmap_df <- as.data.frame(mat_means_scaled) %>%
  rownames_to_column("treatment") %>%
  pivot_longer(-treatment, names_to = "trait", values_to = "value") %>%
  mutate(
    trait = factor(trait, levels = col_order),
    treatment = factor(treatment, levels = rev(row_order))  # so top row is first in row_order
  )

heatmap_plot <- ggplot(heatmap_df, aes(x = trait, y = treatment, fill = value)) +
  geom_tile(color = "grey30") +
  scale_fill_gradient2(low = "#1a9850", mid = "black", high = "#d73027", midpoint = 0,
                       name = "z-score") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),    # hide here; we'll show labels below dendrogram
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 2, r = 6, b = 0, l = 6)
  )

# ---------------------------
# 7) Build column dendrogram plot (traits) to be put on TOP
#    Align x positions to 1..ntraits
# ---------------------------
# For dendro_data we must pass the hc for the set used; if const_cols_means present,
# create combined labels so the dendrogram covers only variable columns; we'll draw labels for all columns below or above.
ddata_cols <- dendro_data(hc_cols, type = "rectangle")
seg_cols <- segment(ddata_cols)
labs_cols <- label(ddata_cols)

# Build a top dendrogram plotting branches upward (default orientation).
# We'll not put trait labels here (we'll enable heatmap x-axis labels) to avoid double printing.
dendro_top <- ggplot() +
  geom_segment(data = seg_cols, aes(x = x, y = y, xend = xend, yend = yend), size = 0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, length(col_order) + 0.5)) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 6, b = 0, l = 6))

# ---------------------------
# 8) Prepare PCA (specimen-level) on correlation (i.e. z-scored specimen data)
# ---------------------------
# Prepare specimen-level trait matrix and remove constant columns
spec_trait_df <- df %>% select(all_of(trait_cols))
sds_spec <- apply(spec_trait_df, 2, sd, na.rm = TRUE)
zero_var_spec <- names(sds_spec)[is.na(sds_spec) | sds_spec <= 1e-12]
if (length(zero_var_spec) > 0) {
  message("Dropping constant/near-constant trait columns for specimen PCA: ", paste(zero_var_spec, collapse = ", "))
}
spec_pca_traits <- setdiff(colnames(spec_trait_df), zero_var_spec)
mat_spec <- df %>% select(all_of(spec_pca_traits)) %>% as.matrix()
rownames(mat_spec) <- df$specimen_id %||% seq_len(nrow(mat_spec))  # specimen IDs if available

# z-score specimen matrix (PCA on correlation matrix)
mat_spec_scaled <- scale(mat_spec, center = TRUE, scale = TRUE)
mat_spec_scaled[is.na(mat_spec_scaled)] <- 0

# Run PCA
pca_spec <- prcomp(mat_spec_scaled, center = FALSE, scale. = FALSE)

# Prepare specimen scores and attach treatment & cluster membership
scores_spec <- as.data.frame(pca_spec$x[, 1:2]) %>%
  rownames_to_column("specimen_row") %>%
  # join specimen metadata (treatment etc.)
  left_join(df %>% select(specimen_row = row_number(), specimen_id, treatment), by = "specimen_row")

# If specimen_id exists, use it for labels; otherwise use row number
if ("specimen_id" %in% colnames(df)) {
  scores_spec <- scores_spec %>% mutate(specimen_label = df$specimen_id)
} else {
  scores_spec <- scores_spec %>% mutate(specimen_label = as.character(specimen_row))
}

# assign each specimen the cluster of its treatment (clusters computed from treatment means)
scores_spec <- scores_spec %>% left_join(treatment_cluster_df, by = "treatment")

# Determine which clusters have >=3 specimens for stat_ellipse
cluster_counts <- scores_spec %>% count(cluster) %>% rename(n = n)
clusters_with_ellipse <- cluster_counts %>% filter(n >= 3) %>% pull(cluster)
message("Clusters with >=3 specimens (ellipses will be drawn): ", paste(clusters_with_ellipse, collapse = ", "))

# Loadings (variables) for PC1/PC2. Use same scaled variables as pca input
loadings <- as.data.frame(pca_spec$rotation[, 1:2], check.names = FALSE) %>% rownames_to_column("variable")
colnames(loadings)[2:3] <- c("PC1", "PC2")

# scale loadings so arrows fit into scores range
mult <- min(
  (max(scores_spec$PC1) - min(scores_spec$PC1)) / (max(loadings$PC1) - min(loadings$PC1) + 1e-9),
  (max(scores_spec$PC2) - min(scores_spec$PC2)) / (max(loadings$PC2) - min(loadings$PC2) + 1e-9)
) * 0.6
loadings <- loadings %>% mutate(PC1s = PC1 * mult, PC2s = PC2 * mult, len = sqrt(PC1^2 + PC2^2))
loadings_label <- loadings %>% arrange(desc(len)) %>% slice(1:min(label_top_loadings, nrow(loadings)))

# Percent variance labels
pct_var <- function(i) round(100 * (pca_spec$sdev[i]^2 / sum(pca_spec$sdev^2)), 2)
xlab <- paste0("PC1 (", pct_var(1), "%)")
ylab <- paste0("PC2 (", pct_var(2), "%)")

# PCA plot: specimens colored by treatment cluster
pca_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(data = scores_spec, aes(x = PC1, y = PC2, color = cluster), size = 2.8) +
  # draw ellipses only for clusters with >=3 specimens (stat_ellipse uses covariance)
  {if (length(clusters_with_ellipse) > 0)
    stat_ellipse(data = filter(scores_spec, cluster %in% clusters_with_ellipse),
                 aes(x = PC1, y = PC2, group = cluster, color = cluster),
                 level = 0.68, linetype = 2, size = 0.8)
  } +
  # labels for specimen points (optional); we keep short and readable
  geom_text_repel(data = scores_spec, aes(x = PC1, y = PC2, label = specimen_label),
                  size = 3, max.overlaps = 30, box.padding = 0.35) +
  # variable arrows
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
               arrow = arrow(length = unit(0.02, "npc")), color = "blue", alpha = 0.6) +
  geom_text_repel(data = loadings_label, aes(x = PC1s, y = PC2s, label = variable),
                  color = "black", size = 3, max.overlaps = 50, box.padding = 0.4) +
  theme_minimal(base_size = 12) +
  labs(x = xlab, y = ylab, color = paste0("treatment cluster (k=", k_clusters, ")")) +
  theme(plot.margin = margin(t = 6, r = 8, b = 6, l = 8))

# ---------------------------
# 9) Create a full column label strip (so we can show trait names under top dendrogram)
#    This draws the column names (rotated) in alignment with the dendrogram/heatmap
# ---------------------------
col_labels <- tibble(trait = col_order, x = seq_along(col_order)) %>%
  mutate(label = trait)

col_label_plot <- ggplot(col_labels, aes(x = x, y = 0, label = label)) +
  geom_text(aes(x = x, y = 0, label = label), angle = 90, hjust = 1, size = 3) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, length(col_order) + 0.5)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 6, b = 0, l = 6))

# We will stack plots: dendro_top / col_label_plot / heatmap_plot / pca_plot
# But dendro_top already uses same x scale (0.5 .. n+0.5). col_label_plot aligns labels.

# ---------------------------
# 10) Combine and save
# ---------------------------
combined <- dendro_top / col_label_plot / heatmap_plot / pca_plot +
  plot_layout(heights = c(0.8, 0.5, 4.5, 4.2))

print(combined)
ggsave(output_file, combined, width = width, height = height, dpi = dpi)
message("Saved combined plot to: ", output_file)