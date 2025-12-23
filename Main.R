
#Loading data

md <- read.csv("metadata2.csv", header = TRUE)
names(md)

summary(metadata2)

#install/loading packages

library(tidyverse)
library(dplyr)
library(dbplyr)
library(ggplot2)
library(pheatmap)
library(car)
library(emmeans)
library(lme4)
library(lmerTest)
library(janitor)

# Data Pipelining and Cleaning
################################

library(ggpubr)

summary(metadata2)

#njchnfcjdnhc

###### DATA CLEANING ######

#import dataset

md <- read.csv("metadata2.csv", header = TRUE)
View(md)

#load library
library(readxl)
library(tidyverse)
library(janitor)

#observe/ inspect data
head(md)
tail(md)
str(md)
summary(md)
###note: NA in columns: LFW,LDW, TW, RWC, WD 
###note: dead plants in 3 rows (specimenid: CX2, CY2, EY2,)

#data types and errors
# conversion to numeric to factor by column "block"
md$block <- as.factor(md$block)
#checking
class(md$block)
levels(md$block)

library(dplyr)
library(janitor)

clean <- md %>%
  clean_names() %>%  
  mutate(
    across(
      c(leaf_area, leaf_number, shoot_dry_weight, root_dry_weight, 
        shoot_fresh_weight, root_fresh_weight, leaf_fresh_weight, 
        leaf_dry_weight, turgid_weight, leaf_rwc, water_deficit, 
        plant_height, stem_diameter, shoot_length, root_length, 
        root_to_shoot_ratio, block_mortality_rate, block_survival_rate),
      ~ as.numeric(.x)
    ),
    specimen_id = as.character(specimen_id),
    block = as.factor(block)
  )

View(clean)

# Remove the dead 
clean <- clean %>%
  filter(
    !(coalesce(leaf_area, 0) == 0 &
        coalesce(leaf_number, 0) == 0 &
        coalesce(shoot_dry_weight, 0) == 0 &
        coalesce(root_dry_weight, 0) == 0)
  )

#verify
class(clean$block)
class(clean$ethanol_pre_treatment)
print(paste("Columns:", ncol(clean)))
print(paste("Rows:", nrow(clean)))
names(clean)
head(clean)
view(clean)


# Data Groupings ----------------------------------------------------------

md_groups <- clean %>%
  mutate(
    # convert to factors with labels
    saline_treatment = factor(
      saline_treatment,
      levels = c(0, 50, 100),
      labels = c("Control", "S1", "S2")
    ),
    ethanol_pre_treatment = factor(
      ethanol_pre_treatment,
      levels = c(0, 20),
      labels = c("None", "Ethanol")
    ),
    
    # create group labels
    group = case_when(
      saline_treatment == "Control" & ethanol_pre_treatment == "None"     ~ "CX",
      saline_treatment == "Control" & ethanol_pre_treatment == "Ethanol"  ~ "EX",
      saline_treatment == "S1"      & ethanol_pre_treatment == "None"     ~ "CY",
      saline_treatment == "S1"      & ethanol_pre_treatment == "Ethanol"  ~ "EY",
      saline_treatment == "S2"      & ethanol_pre_treatment == "None"     ~ "CZ",
      saline_treatment == "S2"      & ethanol_pre_treatment == "Ethanol"  ~ "EZ",
      TRUE ~ NA_character_
    )
  )

table(md_groups$group, useNA = "ifany")

#reference levels
md_groups <- md_groups %>%
  mutate(
    saline_treatment = relevel(saline_treatment, ref = "Control"),
    ethanol_pre_treatment = relevel(ethanol_pre_treatment, ref = "None")
  )

unique(clean$saline_treatment)
unique(clean$ethanol_pre_treatment)

clean <- clean %>% 
  mutate(
    saline_treatment = factor(saline_treatment,
                              levels = c("Control", "S1", "S2")),
    ethanol_pre_treatment = factor(ethanol_pre_treatment,
                                   levels = c("None", "Ethanol"))
  )

table(md_groups$group)

#data-splitting

group_dfs <- md_groups %>%
  filter(!is.na(group)) %>%
  split(.$group)

#assignment of dfs

CX <- group_dfs$CX   
EX <- group_dfs$EX   
CY <- group_dfs$CY   
EY <- group_dfs$EY   
CZ <- group_dfs$CZ   
EZ <- group_dfs$EZ   

#Detecting outlier using IQR-------------------------

##IQR= Q3-Q1
calculate_outliers <- function(data, trait) {
  x <- data[[trait]]
  
  Q1  <- quantile(x, 0.25, na.rm = TRUE)
  Q3  <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  
  cat("\n========================\n")
  cat("Trait:", trait, "\n")
  cat("Q1:", Q1, " Q3:", Q3, " IQR:", IQR_val, "\n")
  cat("Lower bound:", lower_bound, " Upper bound:", upper_bound, "\n")
  
  # treat NA as not an outlier
  out_raw <- x < lower_bound | x > upper_bound
  outliers_logical <- ifelse(is.na(out_raw), FALSE, out_raw)
  
  
  outlier_values <- x[outliers_logical]
  
  cat("Number of outliers:", length(outlier_values), "\n")
  if (length(outlier_values) > 0) {
    cat("Outlier values:", outlier_values, "\n")
  }
  
  return(outliers_logical)
}
trait_cols <- trimws(trait_cols)
# List of all traits
trait_cols <- c(
  "leaf_area", "leaf_number", "shoot_dry_weight", "root_dry_weight", 
  "shoot_fresh_weight", "root_fresh_weight", "leaf_fresh_weight", 
  "leaf_dry_weight", "turgid_weight", "leaf_rwc", "water_deficit", 
  "plant_height", "stem_diameter", "shoot_length", "root_length", 
  "root_to_shoot_ratio"
)

# Apply per trait and see outliers each trait
trait_outliers <- lapply(trait_cols, function(tr) calculate_outliers(clean, tr))

# FLAG outliers inside df (if viewed, additional coloumns, true=outliers, false=normal)
clean_outliers <- clean
for (tr in trait_cols) {
  out_name <- paste0(tr, "_outlier")
  clean_outliers[[out_name]] <- trait_outliers[[tr]]
}

clean_outliers <- clean %>%
  mutate(across(all_of(trait_cols), ~ {
    x <- .x[!is.na(.x)]
    if(length(x) < 4) return(FALSE)  # Skip if too few data
    Q1 <- quantile(x, 0.25)
    Q3 <- quantile(x, 0.75)
    IQR_val <- Q3 - Q1
    ( .x < (Q1 - 1.5*IQR_val) | .x > (Q3 + 1.5*IQR_val) ) & !is.na(.x)
  }, .names = "{.col}_outlier"))

# Sanitize ALL column names once
colnames(clean) <- trimws(colnames(clean))
# visualize trait
boxplot(clean$root_fresh_weight)


# Normality Test per treatment
shapiro.test(CX$leaf_area)
shapiro.test(EX$leaf_area)
shapiro.test(CY$leaf_area)
shapiro.test(EY$leaf_area)
shapiro.test(CZ$leaf_area)
shapiro.test(EZ$leaf_area)

shapiro.test(CX$leaf_number)
shapiro.test(EX$leaf_number)
shapiro.test(CY$leaf_number)
shapiro.test(EY$leaf_number)
shapiro.test(CZ$leaf_number)
shapiro.test(EZ$leaf_number)

shapiro.test(CX$shoot_dry_weight)
shapiro.test(EX$shoot_dry_weight)
shapiro.test(CY$shoot_dry_weight)
shapiro.test(EY$shoot_dry_weight)
shapiro.test(CZ$shoot_dry_weight)
shapiro.test(EZ$shoot_dry_weight)

shapiro.test(CX$root_dry_weight)
shapiro.test(EX$root_dry_weight)
shapiro.test(CY$root_dry_weight)
shapiro.test(EY$root_dry_weight)
shapiro.test(CZ$root_dry_weight)
shapiro.test(EZ$root_dry_weight)

shapiro.test(CX$shoot_fresh_weight)
shapiro.test(EX$shoot_fresh_weight)
shapiro.test(CY$shoot_fresh_weight)
shapiro.test(EY$shoot_fresh_weight)
shapiro.test(CZ$shoot_fresh_weight)
shapiro.test(EZ$shoot_fresh_weight)

shapiro.test(CX$root_fresh_weight)
shapiro.test(EX$root_fresh_weight)
shapiro.test(CY$root_fresh_weight)
shapiro.test(EY$root_fresh_weight)
shapiro.test(CZ$root_fresh_weight)
shapiro.test(EZ$root_fresh_weight)

shapiro.test(CX$plant_height)
shapiro.test(EX$plant_height)
shapiro.test(CY$plant_height)
shapiro.test(EY$plant_height)
shapiro.test(CZ$plant_height)
shapiro.test(EZ$plant_height)

shapiro.test(CX$stem_diameter)
shapiro.test(EX$stem_diameter)
shapiro.test(CY$stem_diameter)
shapiro.test(EY$stem_diameter)
shapiro.test(CZ$stem_diameter)
shapiro.test(EZ$stem_diameter)

shapiro.test(CX$shoot_length)
shapiro.test(EX$shoot_length)
shapiro.test(CY$shoot_length)
shapiro.test(EY$shoot_length)
shapiro.test(CZ$shoot_length)
shapiro.test(EZ$shoot_length)

shapiro.test(CX$root_length)
shapiro.test(EX$root_length)
shapiro.test(CY$root_length)
shapiro.test(EY$root_length)
shapiro.test(CZ$root_length)
shapiro.test(EZ$root_length)

shapiro.test(CX$root_to_shoot_ratio)
shapiro.test(EX$root_to_shoot_ratio)
shapiro.test(CY$root_to_shoot_ratio)
shapiro.test(EY$root_to_shoot_ratio)
shapiro.test(CZ$root_to_shoot_ratio)
shapiro.test(EZ$root_to_shoot_ratio)

shapiro.test(CX$health_status)
shapiro.test(EX$health_status)
shapiro.test(CY$health_status)
shapiro.test(EY$health_status)
shapiro.test(CZ$health_status)
#shapiro.test(EZ$health_status)
# ----------------------------
# Correlational analysis + PCA / Heatmap / Dendrogram (REPLACEMENT BLOCK)
# (Paste this over your existing correlational / PCA section)
# ----------------------------

# Ensure plotting packages are available (do not attach new packages here to avoid masking)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Please install ggrepel")
if (!requireNamespace("ggdendro", quietly = TRUE)) stop("Please install ggdendro")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install patchwork")

# Confirm required objects
if (!exists("clean") || !exists("md_groups")) stop("Objects `clean` and/or `md_groups` not found. Run earlier sections first.")

# Define trait names (trimmed)
trait_cols <- c(
  "leaf_area","leaf_number","shoot_dry_weight","root_dry_weight",
  "shoot_fresh_weight","root_fresh_weight","leaf_fresh_weight",
  "leaf_dry_weight","turgid_weight","leaf_rwc","water_deficit",
  "plant_height","stem_diameter","shoot_length","root_length",
  "root_to_shoot_ratio"
)
trait_cols <- trimws(trait_cols)

# Keep only traits actually present in `clean`
present_traits <- intersect(trait_cols, colnames(clean))
if (length(present_traits) == 0) stop("No trait columns found in `clean` for correlational/PCA analysis.")

# ---- Correlation matrix (base-R indexing) ----
cor_results <- stats::cor(clean[, present_traits, drop = FALSE], use = "pairwise.complete.obs")
message("Pairwise correlation (rounded 3):")
print(round(cor_results, 3))

# ---- PCA + heatmap settings ----
normalize_before_mean <- TRUE
normalization_method <- "zscore"   # options: "zscore","minmax","log1p"
force_normalize_if_not <- TRUE

# helper: normalization check (works on a data.frame or matrix)
check_normalized <- function(df_numeric, method = "zscore", tol_mean = 1e-6, tol_sd = 1e-6, tol_range = 1e-6) {
  df_numeric <- as.data.frame(df_numeric)
  res <- data.frame(
    trait = colnames(df_numeric),
    mean = sapply(df_numeric, function(x) mean(x, na.rm = TRUE)),
    sd   = sapply(df_numeric, function(x) sd(x, na.rm = TRUE)),
    min  = sapply(df_numeric, function(x) min(x, na.rm = TRUE)),
    max  = sapply(df_numeric, function(x) max(x, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  if (method == "zscore") {
    res$is_normal <- (abs(res$mean) <= tol_mean) & (abs(res$sd - 1) <= tol_sd)
    return(list(per_trait = res, all_normal = all(res$is_normal, na.rm = TRUE), method = method))
  } else if (method == "minmax") {
    res$is_normal <- (abs(res$min - 0) <= tol_range) & (abs(res$max - 1) <= tol_range)
    return(list(per_trait = res, all_normal = all(res$is_normal, na.rm = TRUE), method = method))
  } else if (method == "log1p") {
    res$is_nonneg <- res$min >= 0
    warning("log1p check only verifies non-negativity of raw data.")
    return(list(per_trait = res, all_normal = all(res$is_nonneg, na.rm = TRUE), method = method))
  } else stop("Unsupported normalization method.")
}

# helper: apply normalization (returns modified data.frame)
normalize_df <- function(df, trait_cols, method = "zscore") {
  out <- df
  if (method == "zscore") {
    for (c in trait_cols) {
      v <- df[[c]]
      m <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
      if (is.na(s) || s <= .Machine$double.eps) out[[c]] <- ifelse(is.na(v), NA, 0) else out[[c]] <- (v - m)/s
    }
  } else if (method == "minmax") {
    for (c in trait_cols) {
      v <- df[[c]]; mn <- min(v, na.rm = TRUE); mx <- max(v, na.rm = TRUE); rng <- mx - mn
      if (is.na(rng) || rng <= .Machine$double.eps) out[[c]] <- ifelse(is.na(v), NA, 0) else out[[c]] <- (v - mn)/rng
    }
  } else if (method == "log1p") {
    for (c in trait_cols) {
      v <- df[[c]]; mn <- min(v, na.rm = TRUE)
      if (!is.na(mn) && mn < 0) { shift <- abs(mn) + 1e-8; out[[c]] <- log1p(v + shift) } else out[[c]] <- log1p(v)
    }
  } else stop("Unsupported normalization method.")
  out
}

# Build numeric trait candidates safely (base-R)
trait_candidates <- clean[, vapply(clean, is.numeric, logical(1)), drop = FALSE]
# drop all-NA columns
trait_candidates <- trait_candidates[, vapply(trait_candidates, function(x) !all(is.na(x)), logical(1)), drop = FALSE]
trait_cols_pca <- intersect(colnames(trait_candidates), present_traits)
if (length(trait_cols_pca) == 0) stop("No numeric trait columns available for PCA/heatmap.")

# Optional normalization (specimen-level)
if (normalize_before_mean) {
  check_res <- check_normalized(clean[, trait_cols_pca, drop = FALSE], method = normalization_method)
  if (!isTRUE(check_res$all_normal) && force_normalize_if_not) {
    message("Applying specimen-level normalization (", normalization_method, ")")
    clean <- normalize_df(clean, trait_cols_pca, method = normalization_method)
    trait_candidates <- clean[, trait_cols_pca, drop = FALSE]
  } else {
    message("Specimen-level traits appear normalized or no action requested.")
  }
}

# Remove grouping columns if present in trait list (safety)
group_cols <- c("ethanol_pre_treatment", "saline_treatment", "group", "block")
trait_cols_pca <- setdiff(trait_cols_pca, intersect(trait_cols_pca, group_cols))
if (length(trait_cols_pca) == 0) stop("No trait columns left for PCA after removing grouping columns.")

# Aggregate to treatment means using base aggregate (robust, no tidyselect)
agg_by_eth_sal <- aggregate(
  x = clean[, trait_cols_pca, drop = FALSE],
  by = list(ethanol = as.character(md_groups$ethanol_pre_treatment), saline = as.character(md_groups$saline_treatment)),
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Build df_means with proper names
df_means <- as.data.frame(agg_by_eth_sal, stringsAsFactors = FALSE)
# columns: ethanol, saline, then trait columns
names(df_means)[1:2] <- c("ethanol_pre_treatment", "saline_treatment")

# Convert numeric-like labels back if possible
# Map ethanol factor labels "None"/"Ethanol" -> numeric codes, same for saline
df_means$ethanol_val <- ifelse(tolower(df_means$ethanol_pre_treatment) %in% c("none","none "), 0,
                               ifelse(tolower(df_means$ethanol_pre_treatment) %in% c("ethanol","ethanol "), 20, as.numeric(df_means$ethanol_pre_treatment)))
df_means$saline_val <- ifelse(tolower(df_means$saline_treatment) %in% c("control","control "), 0,
                              ifelse(tolower(df_means$saline_treatment) %in% c("s1","s1 "), 50,
                                     ifelse(tolower(df_means$saline_treatment) %in% c("s2","s2 "), 100, as.numeric(df_means$saline_treatment))))

df_means$treatment <- paste0("E", df_means$ethanol_val, "_S", df_means$saline_val)
df_means$treatment_label <- df_means$treatment

# optional mapping to friendly labels
treatment_map <- c("E0_S0"="Control","E0_S50"="S1","E0_S100"="S2","E20_S0"="Eth","E20_S50"="S1 + Eth","E20_S100"="S2 + Eth")
df_means$treatment_label <- ifelse(df_means$treatment %in% names(treatment_map), treatment_map[df_means$treatment], df_means$treatment_label)

# Build matrix for heatmap/PCA (rows = treatments)
mat <- as.matrix(df_means[, trait_cols_pca, drop = FALSE])
rownames(mat) <- df_means$treatment_label

# Guard: need at least one row
if (nrow(mat) < 1) stop("No treatment rows in df_means; aggregation failed.")

# Identify non-constant traits for PCA
sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
nonconst_traits <- names(sds)[!(is.na(sds) | sds <= 1e-12)]

# If fewer than 2 non-constant traits, draw only heatmap
if (length(nonconst_traits) < 2) {
  warning("Not enough non-constant traits for PCA. Drawing heatmap only.")
  mat_scaled_heat <- scale(mat, center = TRUE, scale = TRUE)
  mat_scaled_heat[is.na(mat_scaled_heat)] <- 0

  # Build heatmap long data.frame using base-R (avoid pivot_longer)
  heatmap_df <- do.call(rbind, lapply(seq_len(nrow(mat_scaled_heat)), function(i) {
    data.frame(treatment = rownames(mat_scaled_heat)[i],
               trait = colnames(mat_scaled_heat),
               value = as.numeric(mat_scaled_heat[i, ]),
               stringsAsFactors = FALSE)
  }))

  heatmap_df$trait <- factor(heatmap_df$trait, levels = colnames(mat_scaled_heat))
  heatmap_df$treatment <- factor(heatmap_df$treatment, levels = rev(rownames(mat_scaled_heat)))

  heatmap_plot <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = trait, y = treatment, fill = value)) +
    ggplot2::geom_tile(color = "grey30") +
    ggplot2::scale_fill_gradient2(low = "#00B200", mid = "black", high = "#D7191C", midpoint = 0, name = "z-score") +
    ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title = ggplot2::element_blank())

  print(heatmap_plot)
} else {
  # Full pipeline: heatmap + dendrogram + PCA biplot
  mat_scaled_heat <- scale(mat, center = TRUE, scale = TRUE); mat_scaled_heat[is.na(mat_scaled_heat)] <- 0
  mat_pca <- mat[, nonconst_traits, drop = FALSE]
  mat_scaled_pca <- scale(mat_pca, center = TRUE, scale = TRUE); mat_scaled_pca[is.na(mat_scaled_pca)] <- 0

  # Column clustering (for heatmap)
  hc_cols <- stats::hclust(stats::dist(t(mat_scaled_heat), method = "euclidean"), method = "complete")
  cluster_col_order <- hc_cols$labels[hc_cols$order]

  # Build heatmap long df (base-R)
  col_order <- cluster_col_order
  col_order <- intersect(col_order, colnames(mat_scaled_heat))
  if (length(col_order) != ncol(mat_scaled_heat)) col_order <- c(col_order, setdiff(colnames(mat_scaled_heat), col_order))

  heatmap_df <- do.call(rbind, lapply(seq_len(nrow(mat_scaled_heat)), function(i) {
    data.frame(treatment = rownames(mat_scaled_heat)[i],
               trait = colnames(mat_scaled_heat),
               value = as.numeric(mat_scaled_heat[i, ]),
               stringsAsFactors = FALSE)
  }))
  heatmap_df$trait <- factor(heatmap_df$trait, levels = col_order)
  heatmap_df$treatment <- factor(heatmap_df$treatment, levels = rev(rownames(mat_scaled_heat)))

  heatmap_plot <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = trait, y = treatment, fill = value)) +
    ggplot2::geom_tile(color = "grey30") +
    ggplot2::scale_fill_gradient2(low = "#00B200", mid = "black", high = "#D7191C", midpoint = 0,
                                  limits = c(min(heatmap_df$value, na.rm = TRUE), max(heatmap_df$value, na.rm = TRUE)),
                                  name = "z-score") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
                   axis.ticks.x = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())

  # Dendrogram (top)
  ddata <- tryCatch(ggdendro::dendro_data(hc_cols, type = "rectangle"), error = function(e) NULL)
  if (!is.null(ddata)) {
    segments_df <- ddata$segments
    dendro_plot_top <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), size = 0.6) +
      ggplot2::scale_x_continuous(expand = c(0,0), limits = c(0.5, length(col_order) + 0.5)) +
      ggplot2::theme_void()
  } else {
    dendro_plot_top <- NULL
    warning("Dendrogram construction failed; continuing without dendrogram.")
  }

  # PCA (biplot)
  pca <- stats::prcomp(mat_scaled_pca, center = FALSE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, 1:2], check.names = FALSE); scores <- tibble::rownames_to_column(scores, "treatment_label")
  colnames(scores)[2:3] <- c("PC1","PC2")
  # Attach grouping info by base-R merge if available
  if ("treatment_label" %in% colnames(df_means)) {
    join_cols <- intersect(c("treatment_label","ethanol_pre_treatment","saline_treatment"), colnames(df_means))
    if (length(join_cols) > 1) scores <- merge(scores, df_means[, join_cols, drop = FALSE], by = "treatment_label", all.x = TRUE, sort = FALSE)
  }

  loadings <- as.data.frame(pca$rotation[, 1:2], check.names = FALSE); loadings <- tibble::rownames_to_column(loadings, "trait")
  colnames(loadings)[2:3] <- c("PC1","PC2")
  mult <- min((max(scores$PC1)-min(scores$PC1))/(max(loadings$PC1)-min(loadings$PC1)+1e-9),
              (max(scores$PC2)-min(scores$PC2))/(max(loadings$PC2)-min(loadings$PC2)+1e-9)) * 0.6
  loadings$PC1s <- loadings$PC1 * mult; loadings$PC2s <- loadings$PC2 * mult
  loadings$len <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  loadings_label <- loadings[order(-loadings$len), , drop = FALSE][1:min(18, nrow(loadings)), ]

  pct_var <- function(i) round(100*(pca$sdev[i]^2/sum(pca$sdev^2)), 2)
  xlab <- paste0("PC1 (", pct_var(1), "%)"); ylab <- paste0("PC2 (", pct_var(2), "%)")

  pca_plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, color = "grey70") +
    ggplot2::geom_point(data = scores, ggplot2::aes(PC1, PC2), color = "red", size = 3) +
    ggrepel::geom_text_repel(data = scores, ggplot2::aes(PC1, PC2, label = treatment_label), color = "red", size = 3.2, max.overlaps = Inf) +
    ggplot2::geom_segment(data = loadings, ggplot2::aes(x = 0, y = 0, xend = PC1s, yend = PC2s), arrow = grid::arrow(length = unit(0.02, "npc")), color = "blue", alpha = 0.6) +
    ggrepel::geom_text_repel(data = loadings_label, ggplot2::aes(x = PC1s, y = PC2s, label = trait), size = 3) +
    ggplot2::stat_ellipse(data = scores, ggplot2::aes(x = PC1, y = PC2), linetype = 2, color = "darkred", level = 0.68, inherit.aes = FALSE) +
    ggplot2::theme_minimal() + ggplot2::labs(x = xlab, y = ylab) + ggplot2::theme(plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 8))

  # Combine and print (patchwork)
  if (!is.null(dendro_plot_top)) {
    combined <- (dendro_plot_top / heatmap_plot) / pca_plot + patchwork::plot_layout(heights = c(0.7, 4.5, 4.0))
  } else {
    combined <- (heatmap_plot) / pca_plot + patchwork::plot_layout(heights = c(4.5, 4.0))
  }
  print(combined)
}

shapiro.test(CX$block_mortality_rate)
shapiro.test(EX$block_mortality_rate)
shapiro.test(CY$block_mortality_rate)
shapiro.test(EY$block_mortality_rate)
shapiro.test(CZ$block_mortality_rate)
shapiro.test(EZ$block_mortality_rate)

shapiro.test(CX$block_survival_rate)
shapiro.test(EX$block_survival_rate)
shapiro.test(CY$block_survival_rate)
shapiro.test(EY$block_survival_rate)
shapiro.test(CZ$block_survival_rate)
shapiro.test(EZ$block_survival_rate)


#list of traits that cannot be tested due to limited sample size: leaf fresh weight, lead dw, leaf rwc, water deficit, turgid weight
#shapiro.test(CX$leaf_fresh_weight)  #cant be testesd
#shapiro.test(EX$leaf_rwc)
shapiro.test(CY$turgid_weight)
shapiro.test(EY$leaf_dry_weight)
shapiro.test(CZ$water_deficit)




################summary

traits <- c("leaf_area", "leaf_number", "shoot_dry_weight", "root_dry_weight", 
            "shoot_fresh_weight", "root_fresh_weight", "plant_height", 
            "stem_diameter", "shoot_length", "root_length", "root_to_shoot_ratio",
            "health_status")

normal_counts <- c(6, 2, 6, 5, 3, 1, 6, 4, 5, 5, 4, 3)  # p>0.05 per trait

summary_table <- data.frame(
  Trait = traits,
  Normal = ifelse(normal_counts >= 4, "NORMAL", "NON-NORMAL"),
  Prop_Normal = round(normal_counts/6, 2)
)
print(summary_table)


# Homogeneity testing
fligner.test(list(CX$leaf_area, EX$leaf_area, CY$leaf_area, EY$leaf_area, CZ$leaf_area, EZ$leaf_area))
fligner.test(list(CX$leaf_number, EX$leaf_number, CY$leaf_number, EY$leaf_number, CZ$leaf_number, EZ$leaf_number))
fligner.test(list(CX$shoot_dry_weight, EX$shoot_dry_weight, CY$shoot_dry_weight, EY$shoot_dry_weight, CZ$shoot_dry_weight, EZ$shoot_dry_weight))
fligner.test(list(CX$root_dry_weight, EX$root_dry_weight, CY$root_dry_weight, EY$root_dry_weight, CZ$root_dry_weight, EZ$root_dry_weight))
fligner.test(list(CX$shoot_fresh_weight, EX$shoot_fresh_weight, CY$shoot_fresh_weight, EY$shoot_fresh_weight, CZ$shoot_fresh_weight, EZ$shoot_fresh_weight))
fligner.test(list(CX$root_fresh_weight, EX$root_fresh_weight, CY$root_fresh_weight, EY$root_fresh_weight, CZ$root_fresh_weight, EZ$root_fresh_weight))
fligner.test(list(CX$plant_height, EX$plant_height, CY$plant_height, EY$plant_height, CZ$plant_height, EZ$plant_height))
fligner.test(list(CX$stem_diameter, EX$stem_diameter, CY$stem_diameter, EY$stem_diameter, CZ$stem_diameter, EZ$stem_diameter))
fligner.test(list(CX$shoot_length, EX$shoot_length, CY$shoot_length, EY$shoot_length, CZ$shoot_length, EZ$shoot_length))
fligner.test(list(CX$root_length, EX$root_length, CY$root_length, EY$root_length, CZ$root_length, EZ$root_length))
fligner.test(list(CX$root_to_shoot_ratio, EX$root_to_shoot_ratio, CY$root_to_shoot_ratio, EY$root_to_shoot_ratio, CZ$root_to_shoot_ratio, EZ$root_to_shoot_ratio))
fligner.test(list(CX$health_status, EX$health_status, CY$health_status, EY$health_status, CZ$health_status, EZ$health_status))
fligner.test(list(CX$block_mortality_rate, EX$block_mortality_rate, CY$block_mortality_rate, EY$block_mortality_rate, CZ$block_mortality_rate, EZ$block_mortality_rate))
fligner.test(list(CX$block_survival_rate, EX$block_survival_rate, CY$block_survival_rate, EY$block_survival_rate, CZ$block_survival_rate, EZ$block_survival_rate))


leveneTest(leaf_area ~ group, data = md_groups)
leveneTest(leaf_number ~ group, data = md_groups)
leveneTest(shoot_dry_weight ~ group, data = md_groups)
leveneTest(root_dry_weight ~ group, data = md_groups)
leveneTest(shoot_fresh_weight ~ group, data = md_groups)
leveneTest(root_fresh_weight ~ group, data = md_groups)
leveneTest(plant_height ~ group, data = md_groups)
leveneTest(stem_diameter ~ group, data = md_groups)
leveneTest(shoot_length ~ group, data = md_groups)
leveneTest(root_length ~ group, data = md_groups)
leveneTest(root_to_shoot_ratio ~ group, data = md_groups)
leveneTest(health_status ~ group, data = md_groups)
leveneTest(block_mortality_rate ~ group, data = md_groups)
leveneTest(block_survival_rate ~ group, data = md_groups)

# Descriptive analysis of leaf area
min(CX$leaf_area)
min(EX$leaf_area)
min(CY$leaf_area)
min(EY$leaf_area)
min(CZ$leaf_area)
min(EZ$leaf_area)

max(CX$leaf_area)
max(EX$leaf_area)
max(CY$leaf_area)
max(EY$leaf_area)
max(CZ$leaf_area)
max(EZ$leaf_area)

mean(CX$leaf_area)
mean(EX$leaf_area)
mean(CY$leaf_area)
mean(EY$leaf_area)
mean(CZ$leaf_area)
mean(EZ$leaf_area)

median(CX$leaf_area)
median(EX$leaf_area)
median(CY$leaf_area)
median(EY$leaf_area)
median(CZ$leaf_area)
median(EZ$leaf_area)

sd(CX$leaf_area)
sd(EX$leaf_area)
sd(CY$leaf_area)
sd(EY$leaf_area)
sd(CZ$leaf_area)
sd(EZ$leaf_area)

var(CX$leaf_area)
var(EX$leaf_area)
var(CY$leaf_area)
var(EY$leaf_area)
var(CZ$leaf_area)
var(EZ$leaf_area)

# Descriptive analysis of leaf number
min(CX$leaf_area)
min(EX$leaf_area)
min(CY$leaf_area)
min(EY$leaf_area)
min(CZ$leaf_area)
min(EZ$leaf_area)

max(CX$leaf_number)
max(EX$leaf_number)
max(CY$leaf_number)
max(EY$leaf_number)
max(CZ$leaf_number)
max(EZ$leaf_number)

mean(CX$leaf_number)
mean(EX$leaf_number)
mean(CY$leaf_number)
mean(EY$leaf_number)
mean(CZ$leaf_number)
mean(EZ$leaf_number)

median(CX$leaf_number)
median(EX$leaf_number)
median(CY$leaf_number)
median(EY$leaf_number)
median(CZ$leaf_number)
median(EZ$leaf_number)

sd(CX$leaf_number)
sd(EX$leaf_number)
sd(CY$leaf_number)
sd(EY$leaf_number)
sd(CZ$leaf_number)
sd(EZ$leaf_number)

var(CX$leaf_number)
var(EX$leaf_number)
var(CY$leaf_number)
var(EY$leaf_number)
var(CZ$leaf_number)
var(EZ$leaf_number)

# Descriptive analysis of shoot dry weight

min(CX$shoot_dry_weight)
min(EX$shoot_dry_weight)
min(CY$shoot_dry_weight)
min(EY$shoot_dry_weight)
min(CZ$shoot_dry_weight)
min(EZ$shoot_dry_weight)

max(CX$shoot_dry_weight)
max(EX$shoot_dry_weight)
max(CY$shoot_dry_weight)
max(EY$shoot_dry_weight)
max(CZ$shoot_dry_weight)
max(EZ$shoot_dry_weight)

mean(CX$shoot_dry_weight)
mean(EX$shoot_dry_weight)
mean(CY$shoot_dry_weight)
mean(EY$shoot_dry_weight)
mean(CZ$shoot_dry_weight)
mean(EZ$shoot_dry_weight)

median(CX$shoot_dry_weight)
median(EX$shoot_dry_weight)
median(CY$shoot_dry_weight)
median(EY$shoot_dry_weight)
median(CZ$shoot_dry_weight)
median(EZ$shoot_dry_weight)

sd(CX$shoot_dry_weight)
sd(EX$shoot_dry_weight)
sd(CY$shoot_dry_weight)
sd(EY$shoot_dry_weight)
sd(CZ$shoot_dry_weight)
sd(EZ$shoot_dry_weight)

var(CX$shoot_dry_weight)
var(EX$shoot_dry_weight)
var(CY$shoot_dry_weight)
var(EY$shoot_dry_weight)
var(CZ$shoot_dry_weight)
var(EZ$shoot_dry_weight)

# Descriptive analysis of root dry weight

min(CX$root_dry_weight)
min(EX$root_dry_weight)
min(CY$root_dry_weight)
min(EY$root_dry_weight)
min(CZ$root_dry_weight)
min(EZ$root_dry_weight)

max(CX$root_dry_weight)
max(EX$root_dry_weight)
max(CY$root_dry_weight)
max(EY$root_dry_weight)
max(CZ$root_dry_weight)
max(EZ$root_dry_weight)

mean(CX$root_dry_weight)
mean(EX$root_dry_weight)
mean(CY$root_dry_weight)
mean(EY$root_dry_weight)
mean(CZ$root_dry_weight)
mean(EZ$root_d_weight)

median(CX$root_dry_weight)
median(EX$root_dry_weight)
median(CY$root_dry_weight)
median(EY$root_dry_weight)
median(CZ$root_dry_weight)
median(EZ$root_dry_weight)

sd(CX$root_dry_weight)
sd(EX$root_dry_weight)
sd(CY$root_dry_weight)
sd(EY$root_dry_weight)
sd(CZ$root_dry_weight)
sd(EZ$root_dry_weight)

var(CX$root_dry_weight)
var(EX$root_dry_weight)
var(CY$root_dry_weight)
var(EY$root_dry_weight)
var(CZ$root_dry_weight)
var(EZ$root_dry_weight)

# Descriptive analysis of shoot fresh weight

min(CX$shoot_fresh_weight)
min(EX$shoot_fresh_weight)
min(CY$shoot_fresh_weight)
min(EY$shoot_fresh_weight)
min(CZ$shoot_fresh_weight)
min(EZ$shoot_fresh_weight)

max(CX$shoot_fresh_weight)
max(EX$shoot_fresh_weight)
max(CY$shoot_fresh_weight)
max(EY$shoot_fresh_weight)
max(CZ$shoot_fresh_weight)
max(EZ$shoot_fresh_weight)

mean(CX$shoot_fresh_weight)
mean(EX$shoot_fresh_weight)
mean(CY$shoot_fresh_weight)
mean(EY$shoot_fresh_weight)
mean(CZ$shoot_fresh_weight)
mean(EZ$shoot_fresh_weight)

median(CX$shoot_fresh_weight)
median(EX$shoot_fresh_weight)
median(CY$shoot_fresh_weight)
median(EY$shoot_fresh_weight)
median(CZ$shoot_fresh_weight)
median(EZ$shoot_fresh_weight)

sd(CX$shoot_fresh_weight)
sd(EX$shoot_fresh_weight)
sd(CY$shoot_fresh_weight)
sd(EY$shoot_fresh_weight)
sd(CZ$shoot_fresh_weight)
sd(EZ$shoot_fresh_weight)

var(CX$shoot_fresh_weight)
var(EX$shoot_fresh_weight)
var(CY$shoot_fresh_weight)
var(EY$shoot_fresh_weight)
var(CZ$shoot_fresh_weight)
var(EZ$shoot_fresh_weight)

# Descriptive analysis of root fresh weight

min(CX$root_fresh_weight)
min(EX$root_fresh_weight)
min(CY$root_fresh_weight)
min(EY$root_fresh_weight)
min(CZ$root_fresh_weight)
min(EZ$root_fresh_weight)

max(CX$root_fresh_weight)
max(EX$root_fresh_weight)
max(CY$root_fresh_weight)
max(EY$root_fresh_weight)
max(CZ$root_fresh_weight)
max(EZ$root_fresh_weight)

mean(CX$root_fresh_weight)
mean(EX$root_fresh_weight)
mean(CY$root_fresh_weight)
mean(EY$root_fresh_weight)
mean(CZ$root_fresh_weight)
mean(EZ$root_fresh_weight)

median(CX$root_fresh_weight)
median(EX$root_fresh_weight)
median(CY$root_fresh_weight)
median(EY$root_fresh_weight)
median(CZ$root_fresh_weight)
median(EZ$root_fresh_weight)

sd(CX$root_fresh_weight)
sd(EX$root_fresh_weight)
sd(CY$root_fresh_weight)
sd(EY$root_fresh_weight)
sd(CZ$root_fresh_weight)
sd(EZ$root_fresh_weight)

var(CX$root_fresh_weight)
var(EX$root_fresh_weight)
var(CY$root_fresh_weight)
var(EY$root_fresh_weight)
var(CZ$root_fresh_weight)
var(EZ$root_fresh_weight)

# Descriptive analysis of plant height

min(CX$plant_height)
min(EX$plant_height)
min(CY$plant_height)
min(EY$plant_height)
min(CZ$plant_height)
min(EZ$plant_height)

max(CX$plant_height)
max(EX$plant_height)
max(CY$plant_height)
max(EY$plant_height)
max(CZ$plant_height)
max(EZ$plant_height)

mean(CX$plant_height)
mean(EX$plant_height)
mean(CY$plant_height)
mean(EY$plant_height)
mean(CZ$plant_height)
mean(EZ$plant_height)

median(CX$plant_height)
median(EX$plant_height)
median(CY$plant_height)
median(EY$plant_height)
median(CZ$plant_height)
median(EZ$plant_height)

sd(CX$plant_height)
sd(EX$plant_height)
sd(CY$plant_height)
sd(EY$plant_height)
sd(CZ$plant_height)
sd(EZ$plant_height)

var(CX$plant_height)
var(EX$plant_height)
var(CY$plant_height)
var(EY$plant_height)
var(CZ$plant_height)
var(EZ$plant_height)

# Descriptive analysis of stem diameter

min(CX$stem_diameter)
min(EX$stem_diameter)
min(CY$stem_diameter)
min(EY$stem_diameter)
min(CZ$stem_diameter)
min(EZ$stem_diameter)

max(CX$stem_diameter)
max(EX$stem_diameter)
max(CY$stem_diameter)
max(EY$stem_diameter)
max(CZ$stem_diameter)
max(EZ$stem_diameter)

mean(CX$stem_diameter)
mean(EX$stem_diameter)
mean(CY$stem_diameter)
mean(EY$stem_diameter)
mean(CZ$stem_diameter)
mean(EZ$stem_diameter)

median(CX$stem_diameter)
median(EX$stem_diameter)
median(CY$stem_diameter)
median(EY$stem_diameter)
median(CZ$stem_diameter)
median(EZ$stem_diameter)

sd(CX$stem_diameter)
sd(EX$stem_diameter)
sd(CY$stem_diameter)
sd(EY$stem_diameter)
sd(CZ$stem_diameter)
sd(EZ$stem_diameter)

var(CX$stem_diameter)
var(EX$stem_diameter)
var(CY$stem_diameter)
var(EY$stem_diameter)
var(CZ$stem_diameter)
var(EZ$stem_diameter)

# Descriptive analysis of shoot length

min(CX$shoot_length)
min(EX$shoot_length)
min(CY$shoot_length)
min(EY$shoot_length)
min(CZ$shoot_length)
min(EZ$shoot_length)

max(CX$shoot_length)
max(EX$shoot_length)
max(CY$shoot_length)
max(EY$shoot_length)
max(CZ$shoot_length)
max(EZ$shoot_length)

mean(CX$shoot_length)
mean(EX$shoot_length)
mean(CY$shoot_length)
mean(EY$shoot_length)
mean(CZ$shoot_length)
mean(EZ$shoot_length)

median(CX$shoot_length)
median(EX$shoot_length)
median(CY$shoot_length)
median(EY$shoot_length)
median(CZ$shoot_length)
median(EZ$shoot_length)

sd(CX$shoot_length)
sd(EX$shoot_length)
sd(CY$shoot_length)
sd(EY$shoot_length)
sd(CZ$shoot_length)
sd(EZ$shoot_length)

var(CX$shoot_length)
var(EX$shoot_length)
var(CY$shoot_length)
var(EY$shoot_length)
var(CZ$shoot_length)
var(EZ$shoot_length)

# Descriptive analysis of root length

min(CX$root_length)
min(EX$root_length)
min(CY$root_length)
min(EY$root_length)
min(CZ$root_length)
min(EZ$root_length)

max(CX$root_length)
max(EX$root_length)
max(CY$root_length)
max(EY$root_length)
max(CZ$root_length)
max(EZ$root_length)

mean(CX$root_length)
mean(EX$root_length)
mean(CY$root_length)
mean(EY$root_length)
mean(CZ$root_length)
mean(EZ$root_length)

median(CX$root_length)
median(EX$root_length)
median(CY$root_length)
median(EY$root_length)
median(CZ$root_length)
median(EZ$root_length)

sd(CX$root_length)
sd(EX$root_length)
sd(CY$root_length)
sd(EY$root_length)
sd(CZ$root_length)
sd(EZ$root_length)

var(CX$root_length)
var(EX$root_length)
var(CY$root_length)
var(EY$root_length)
var(CZ$root_length)
var(EZ$root_length)

# Descriptive analysis of root to shoot ratio

min(CX$root_to_shoot_ratio)
min(EX$root_to_shoot_ratio)
min(CY$root_to_shoot_ratio)
min(EY$root_to_shoot_ratio)
min(CZ$root_to_shoot_ratio)
min(EZ$root_to_shoot_ratio)

max(CX$root_to_shoot_ratio)
max(EX$root_to_shoot_ratio)
max(CY$root_to_shoot_ratio)
max(EY$root_to_shoot_ratio)
max(CZ$root_to_shoot_ratio)
max(EZ$root_to_shoot_ratio)

mean(CX$root_to_shoot_ratio)
mean(EX$root_to_shoot_ratio)
mean(CY$root_to_shoot_ratio)
mean(EY$root_to_shoot_ratio)
mean(CZ$root_to_shoot_ratio)
mean(EZ$root_to_shoot_ratio)

median(CX$root_to_shoot_ratio)
median(EX$root_to_shoot_ratio)
median(CY$root_to_shoot_ratio)
median(EY$root_to_shoot_ratio)
median(CZ$root_to_shoot_ratio)
median(EZ$root_to_shoot_ratio)

sd(CX$root_to_shoot_ratio)
sd(EX$root_to_shoot_ratio)
sd(CY$root_to_shoot_ratio)
sd(EY$root_to_shoot_ratio)
sd(CZ$root_to_shoot_ratio)
sd(EZ$root_to_shoot_ratio)

var(CX$root_to_shoot_ratio)
var(EX$root_to_shoot_ratio)
var(CY$root_to_shoot_ratio)
var(EY$root_to_shoot_ratio)
var(CZ$root_to_shoot_ratio)
var(EZ$root_to_shoot_ratio)

# Descriptive analysis per treatment

summary(md)
summary(CX)
summary(EX)
summary(CY)
summary(EY)
summary(CZ)
summary(EZ)

# Boxplot per treatment

boxplot(leaf_area ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Leaf Area",
        main = "Leaf Area by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(leaf_number ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Leaf Number",
        main = "Leaf Number by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(shoot_dry_weight ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Shoot Dry Weight",
        main = "Shoot Dry Weight by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(root_dry_weight ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Root Dry Weight",
        main = "Root Dry Weight by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(shoot_fresh_weight ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Shoot Fresh Weight",
        main = "Shoot Fresh Weight by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(root_fresh_weight ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Root Fresh Weight",
        main = "Root Fresh Weight by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(root_fresh_weight ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Root Fresh Weight",
        main = "Root Fresh Weight by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(plant_height ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Plant Height",
        main = "Plant Height by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(stem_diameter ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Stem Diameter",
        main = "Stem Diameter by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(shoot_length ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Shoot Length",
        main = "Shoot Length by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(root_length ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Root Length",
        main = "Root Length by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

boxplot(root_to_shoot_ratio ~ group,
        data = md_groups,
        xlab = "Treatment Group",
        ylab = "Root to Shoot Ratio",
        main = "Root to Shoot Ratio by Treatment",
        col = c("lightblue", "orange", "lightgreen", "pink", "purple", "yellow"))

# Bar plot by traits
library(dplyr)
library(ggplot2)

# Summarize data of leaf area (mean ± SE)
leaf_area_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(leaf_area, na.rm = TRUE),
    se   = sd(leaf_area, na.rm = TRUE) / sqrt(n())
  )

# Bar plot for leaf area
ggplot(leaf_area_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Leaf Area by Treatment",
    x = "Treatment Group",
    y = "Leaf Area"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of leaf number (mean ± SE)
leaf_number_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(leaf_number, na.rm = TRUE),
    se   = sd(leaf_number, na.rm = TRUE) / sqrt(n())
  )

# Barplot for leaf number
ggplot(leaf_number_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Leaf Number by Treatment",
    x = "Treatment Group",
    y = "Leaf Number"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of shoot dry weight (mean ± SE)
shoot_dry_weight_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(shoot_dry_weight, na.rm = TRUE),
    se   = sd(shoot_dry_weight, na.rm = TRUE) / sqrt(n())
  )

# Barplot for shoot dry weight
ggplot(shoot_dry_weight_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Shoot Dry Weight by Treatment",
    x = "Treatment Group",
    y = "Shoot Dry Weight"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of root dry weight (mean ± SE)
root_dry_weight_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(root_dry_weight, na.rm = TRUE),
    se   = sd(root_dry_weight, na.rm = TRUE) / sqrt(n())
  )

# Barplot for root dry weight
ggplot(root_dry_weight_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Root Dry Weight by Treatment",
    x = "Treatment Group",
    y = "Root Dry Weight"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of shoot fresh weight (mean ± SE)
shoot_fresh_weight_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(shoot_fresh_weight, na.rm = TRUE),
    se   = sd(shoot_fresh_weight, na.rm = TRUE) / sqrt(n())
  )

# Barplot for shoot fresh weight
ggplot(shoot_fresh_weight_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Shoot Fresh Weight by Treatment",
    x = "Treatment Group",
    y = "Shoot Fresh Weight"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of root fresh weight (mean ± SE)
root_fresh_weight_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(root_fresh_weight, na.rm = TRUE),
    se   = sd(root_fresh_weight, na.rm = TRUE) / sqrt(n())
  )

# Barplot for root fresh weight
ggplot(root_fresh_weight_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Root Fresh Weight by Treatment",
    x = "Treatment Group",
    y = "Root Fresh Weight"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of plant height (mean ± SE)
plant_height_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(plant_height, na.rm = TRUE),
    se   = sd(plant_height, na.rm = TRUE) / sqrt(n())
  )

# Barplot for plant height
ggplot(plant_height_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Plant Height by Treatment",
    x = "Treatment Group",
    y = "Plant Height"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of stem diameter (mean ± SE)
stem_diameter_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(stem_diameter, na.rm = TRUE),
    se   = sd(stem_diameter, na.rm = TRUE) / sqrt(n())
  )

# Barplot for stem diameter 
ggplot(stem_diameter_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Stem Diameter by Treatment",
    x = "Treatment Group",
    y = "Stem Diameter"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of shoot length (mean ± SE)
shoot_length_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(shoot_length, na.rm = TRUE),
    se   = sd(shoot_length, na.rm = TRUE) / sqrt(n())
  )

# Barplot for shoot length
ggplot(shoot_length_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Shoot Length by Treatment",
    x = "Treatment Group",
    y = "Shoot Length"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of root length (mean ± SE)
root_length_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(root_length, na.rm = TRUE),
    se   = sd(root_length, na.rm = TRUE) / sqrt(n())
  )

# Barplot for root length
ggplot(root_length_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Root Length by Treatment",
    x = "Treatment Group",
    y = "Root Length"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )

# Summarize data of root to shoot ratio (mean ± SE)
root_to_shoot_ratio_summary <- md_groups %>%
  group_by(group) %>%
  summarise(
    mean = mean(root_to_shoot_ratio, na.rm = TRUE),
    se   = sd(root_to_shoot_ratio, na.rm = TRUE) / sqrt(n())
  )

# Barplot for root to shoot ratio
ggplot(root_to_shoot_ratio_summary, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.7) +
  labs(
    title = "Mean Root to Shoot Ratio by Treatment",
    x = "Treatment Group",
    y = "Root to Shoot Ratio"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )




# =================== INFERENTIAL STATISTICS ==========================

library(tidyverse)
library(car)
library(emmeans)
library(MASS)
library(janitor)


md_groups$group <- factor(md_groups$group, levels = c("CX","EX","CY","EY","CZ","EZ"))
md_groups$block <- factor(md_groups$block)

# 1. Meeting assumptions; Transformation of Data using log and box cox (Non-normal traits: Leaf_number, sfw, rfw)

  #Transform shoot fresh weight to normal using log
md_groups$log_shoot_fresh_weight <- log(md_groups$shoot_fresh_weight + 0.01)

  #Transform root fresh weight and leaf number to normal using box cox with plot
# Leaf number 
model_leaf <- lm((leaf_number + 0.5) ~ group + block, data = md_groups)
bc_leaf <- boxcox(model_leaf); lambda_leaf <- bc_leaf$x[which.max(bc_leaf$y)]
md_groups$bc_leaf_number <- ((md_groups$leaf_number + 0.5)^lambda_leaf - 1)/lambda_leaf

#qqplot
par(mfrow=c(1,2)); 
qqnorm(residuals(lm(leaf_number~group+block,md_groups)),main="Leaf Orig");qqline(residuals(lm(leaf_number~group+block,md_groups)),col="red")
qqnorm(residuals(lm(bc_leaf_number~group+block,md_groups)),main=paste("Leaf λ=",round(lambda_leaf,2)));qqline(residuals(lm(bc_leaf_number~group+block,md_groups)),col="blue")
par(mfrow=c(1,1))

# Root fresh weight 
model_rfw <- lm(root_fresh_weight ~ group + block, data = md_groups)
bc_rfw <- boxcox(model_rfw); lambda_rfw <- bc_rfw$x[which.max(bc_rfw$y)]
md_groups$bc_root_fresh_weight <- (md_groups$root_fresh_weight^lambda_rfw - 1)/lambda_rfw  

#qqplot 
par(mfrow=c(1,2))
qqnorm(residuals(lm(root_fresh_weight~group+block,md_groups)), main="RFW - Original"); qqline(residuals(lm(root_fresh_weight~group+block,md_groups)), col="red")
qqnorm(residuals(lm(bc_root_fresh_weight~group+block,md_groups)), main=paste("RFW Box-Cox λ=",round(lambda_rfw,2))); qqline(residuals(lm(bc_root_fresh_weight~group+block,md_groups)), col="blue")
par(mfrow=c(1,1))



#2.check normality again with transformed data: summary table for all traits showing "NORMAL" or NON-NORMAL"

all_traits <- c("leaf_area", "bc_leaf_number", "shoot_dry_weight", "root_dry_weight", 
                "log_shoot_fresh_weight", "bc_root_fresh_weight", "plant_height", 
                "stem_diameter", "shoot_length", "root_length", "root_to_shoot_ratio")

normality_table <- data.frame(
  Trait = all_traits,
  Normal_Groups = sapply(all_traits, function(trait) {
    count <- 0
    for(g in c("CX","EX","CY","EY","CZ","EZ")) {
      x <- na.omit(md_groups[[trait]][md_groups$group == g])
      if(length(x) >= 3 && shapiro.test(x)$p.value > 0.05) count <- count + 1
    }
    count
  }),
  Prop = sapply(all_traits, function(trait) paste0(
    sum(sapply(c("CX","EX","CY","EY","CZ","EZ"), function(g){
      x <- na.omit(md_groups[[trait]][md_groups$group == g]); length(x)>=3 && shapiro.test(x)$p.value>0.05
    })), "/6")),
  Status = sapply(all_traits, function(trait) {
    n <- sum(sapply(c("CX","EX","CY","EY","CZ","EZ"), function(g){
      x <- na.omit(md_groups[[trait]][md_groups$group == g]); length(x)>=3 && shapiro.test(x)$p.value>0.05
    }))
    ifelse(n>=4, "NORMAL", "NON-NORMAL")
  })
); print(normality_table)


#make dataframe "norm_traits"

norm_traits <- md_groups[, c("leaf_area", "bc_leaf_number", "log_shoot_fresh_weight", 
                             "bc_root_fresh_weight", "shoot_dry_weight", "root_dry_weight", 
                             "plant_height", "stem_diameter", "shoot_length", "root_length", 
                             "root_to_shoot_ratio")]

norm_traits <- norm_traits[complete.cases(norm_traits), ]
norm_traits[, c("bc_leaf_number", "log_shoot_fresh_weight", "bc_root_fresh_weight")] <- 
  round(norm_traits[, c("bc_leaf_number", "log_shoot_fresh_weight", "bc_root_fresh_weight")], 4)

cat("\n=== NORM_TRAITS READY (NA-FREE, NO group/block) ===\n")
cat("Dimensions:", nrow(norm_traits), "rows x", ncol(norm_traits), "traits\n")
print(head(norm_traits))

# DEFINE traits_final
traits_final <- all_traits



# 3. RCBD one-way anova 

#fit model
rcbd_models <- lapply(traits_final, function(trait) {
  formula <- as.formula(paste(trait, "~ group + block"))
  lm(formula, data = md_groups)
})

#apply to all traits
rcbd_anova_results <- lapply(rcbd_models, anova)
names(rcbd_anova_results) <- traits_final

#view result to each trait
cat("RCBD ANOVA RESULTS PER TRAIT:\n")
for(i in seq_along(rcbd_anova_results)) {
  cat("\n=== ", traits_final[i], " ===\n")
  print(rcbd_anova_results[[i]])
}



# 4. Post hoc test for significant traits (p<0.05 ONLY). 

cat("\n=== POST-HOC TUKEYHSD (Leaf Area & Leaf Number Only) ===\n")

# Leaf Area
model_leaf_area <- aov(leaf_area ~ group + block, data = md_groups)
cat("Leaf Area ANOVA p =", format.pval(anova(model_leaf_area)$"Pr(>F)"[1]), "\n")
if(anova(model_leaf_area)$"Pr(>F)"[1] < 0.05) {
  cat("Leaf Area TukeyHSD:\n")
  print(TukeyHSD(model_leaf_area)$group[, "p adj"])
}

cat("\n")

# Leaf Number (transformed)
model_leaf_num <- aov(bc_leaf_number ~ group + block, data = md_groups)
cat("Leaf Number ANOVA p =", format.pval(anova(model_leaf_num)$"Pr(>F)"[1]), "\n")
if(anova(model_leaf_num)$"Pr(>F)"[1] < 0.05) {
  cat("Leaf Number TukeyHSD:\n")
  print(TukeyHSD(model_leaf_num)$group[, "p adj"])
}

#============================================================

# Correlational Analysis --------------------------------------------------

library(tidyverse)
library(janitor)
library(ggdendro)
library(ggrepel)
library(patchwork)
library(scales)
library(grid)

set.seed(42)

# ----------------------------
# 0) Read data
# ----------------------------
md <- read.csv("metadata2.csv", header = TRUE, stringsAsFactors = FALSE)

# ----------------------------
# 1) Basic cleaning (snake_case column names, numeric coercion)
# ----------------------------
clean <- md %>%
  janitor::clean_names() %>%
  mutate(
    across(
      c(leaf_area, leaf_number, shoot_dry_weight, root_dry_weight,
        shoot_fresh_weight, root_fresh_weight, leaf_fresh_weight,
        leaf_dry_weight, turgid_weight, leaf_rwc, water_deficit,
        plant_height, stem_diameter, shoot_length, root_length,
        root_to_shoot_ratio, block_mortality_rate, block_survival_rate),
      ~ as.numeric(.x)
    ),
    specimen_id = as.character(specimen_id),
    block = as.factor(block)
  )

# Remove rows that appear to be dead (all zeros in main measured cols)
clean <- clean %>%
  filter(
    !(coalesce(leaf_area, 0) == 0 &
        coalesce(leaf_number, 0) == 0 &
        coalesce(shoot_dry_weight, 0) == 0 &
        coalesce(root_dry_weight, 0) == 0)
  )

# ----------------------------
# 2) Group/factor creation used elsewhere
# ----------------------------
md_groups <- clean %>%
  mutate(
    saline_treatment = factor(saline_treatment, levels = c(0, 50, 100), labels = c("Control", "S1", "S2")),
    ethanol_pre_treatment = factor(ethanol_pre_treatment, levels = c(0, 20), labels = c("None", "Ethanol")),
    group = case_when(
      saline_treatment == "Control" & ethanol_pre_treatment == "None"     ~ "CX",
      saline_treatment == "Control" & ethanol_pre_treatment == "Ethanol"  ~ "EX",
      saline_treatment == "S1"      & ethanol_pre_treatment == "None"     ~ "CY",
      saline_treatment == "S1"      & ethanol_pre_treatment == "Ethanol"  ~ "EY",
      saline_treatment == "S2"      & ethanol_pre_treatment == "None"     ~ "CZ",
      saline_treatment == "S2"      & ethanol_pre_treatment == "Ethanol"  ~ "EZ",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    saline_treatment = relevel(saline_treatment, ref = "Control"),
    ethanol_pre_treatment = relevel(ethanol_pre_treatment, ref = "None")
  )

# split per-group if needed later
group_dfs <- md_groups %>% filter(!is.na(group)) %>% split(.$group)
# e.g. CX <- group_dfs$CX etc.

# ----------------------------
# 3) Outlier detection (IQR) - robust trait list (no trailing spaces)
# ----------------------------
trait_cols <- c(
  "leaf_area", "leaf_number", "shoot_dry_weight", "root_dry_weight",
  "shoot_fresh_weight", "root_fresh_weight", "leaf_fresh_weight",
  "leaf_dry_weight", "turgid_weight", "leaf_rwc", "water_deficit",
  "plant_height", "stem_diameter", "shoot_length", "root_length",
  "root_to_shoot_ratio"
)

# keep only traits present in data
present_traits <- intersect(trait_cols, names(clean))
missing_traits <- setdiff(trait_cols, present_traits)
if (length(missing_traits) > 0) warning("Missing trait columns (skipped): ", paste(missing_traits, collapse = ", "))
if (length(present_traits) == 0) stop("No trait columns found in data.")

calculate_outliers <- function(data, trait) {
  x <- data[[trait]]
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR_val
  upper <- Q3 + 1.5 * IQR_val
  out_raw <- x < lower | x > upper
  ifelse(is.na(out_raw), FALSE, out_raw)
}

trait_outliers <- lapply(present_traits, function(tr) calculate_outliers(clean, tr))
names(trait_outliers) <- present_traits

clean_outliers <- clean
for (tr in names(trait_outliers)) {
  clean_outliers[[paste0(tr, "_outlier")]] <- trait_outliers[[tr]]
}

# ----------------------------
# Correlational analysis placeholder
# ----------------------------
# (User's correlational analysis should be placed here.)
# Example:
# cor_results <- cor(clean %>% select(all_of(present_traits)), use = "pairwise.complete.obs")
# print(cor_results)

# ----------------------------
# PCA + Heatmap pipeline (adapted from PCA2.2.R)
# ----------------------------

# User settings
normalize_before_mean <- TRUE       # normalize per-sample before aggregating
normalization_method <- "zscore"    # "zscore", "minmax", or "log1p"
force_normalize_if_not <- TRUE

# helper functions for normalization/checking
check_normalized <- function(df_numeric, method = "zscore", tol_mean = 1e-6, tol_sd = 1e-6, tol_range = 1e-6) {
  df_numeric <- as.data.frame(df_numeric)
  res <- tibble::tibble(
    trait = colnames(df_numeric),
    mean = sapply(df_numeric, function(x) mean(x, na.rm = TRUE)),
    sd = sapply(df_numeric, function(x) sd(x, na.rm = TRUE)),
    min = sapply(df_numeric, function(x) min(x, na.rm = TRUE)),
    max = sapply(df_numeric, function(x) max(x, na.rm = TRUE))
  )
  if (method == "zscore") {
    res <- res %>% mutate(is_normal = (abs(mean) <= tol_mean) & (abs(sd - 1) <= tol_sd))
    return(list(per_trait = res, all_normal = all(res$is_normal, na.rm = TRUE), method = "zscore"))
  } else if (method == "minmax") {
    res <- res %>% mutate(is_normal = (abs(min - 0) <= tol_range) & (abs(max - 1) <= tol_range))
    return(list(per_trait = res, all_normal = all(res$is_normal, na.rm = TRUE), method = "minmax"))
  } else if (method == "log1p") {
    res <- res %>% mutate(is_nonnegative = min >= 0)
    warning("log1p check: only non-negativity verified.")
    return(list(per_trait = res, all_normal = all(res$is_nonnegative, na.rm = TRUE), method = "log1p"))
  } else stop("Unsupported method.")
}

normalize_df <- function(df, trait_cols, method = "zscore") {
  df_out <- df
  if (method == "zscore") {
    for (c in trait_cols) {
      v <- df[[c]]
      m <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
      if (is.na(s) || s <= .Machine$double.eps) df_out[[c]] <- ifelse(is.na(v), NA, 0) else df_out[[c]] <- (v - m) / s
    }
  } else if (method == "minmax") {
    for (c in trait_cols) {
      v <- df[[c]]; mn <- min(v, na.rm = TRUE); mx <- max(v, na.rm = TRUE); rng <- mx - mn
      if (is.na(rng) || rng <= .Machine$double.eps) df_out[[c]] <- ifelse(is.na(v), NA, 0) else df_out[[c]] <- (v - mn) / rng
    }
  } else if (method == "log1p") {
    for (c in trait_cols) {
      v <- df[[c]]; mn <- min(v, na.rm = TRUE)
      if (!is.na(mn) && mn < 0) { shift <- abs(mn) + 1e-8; df_out[[c]] <- log1p(v + shift) } else df_out[[c]] <- log1p(v)
    }
  } else stop("Unsupported normalization method.")
  df_out
}

# trait columns selected for PCA/heatmap (numeric, non-all-NA)
exclude_cols_pca <- c("specimen_id", "species", "variety_name", "variety_type",
                      "block", "health_status", "block_mortality_rate", "block_survival_rate",
                      "group", "ethanol_pre_treatment", "saline_treatment")
trait_candidates <- clean %>% select(-any_of(exclude_cols_pca)) %>% select(where(is.numeric))
trait_candidates <- trait_candidates %>% select(where(~ !all(is.na(.))))
trait_cols_pca <- colnames(trait_candidates)
if (length(trait_cols_pca) == 0) stop("No numeric trait columns found for PCA/heatmap.")

# Optional: normalize per-sample before aggregating to treatment means
if (normalize_before_mean) {
  check_res <- check_normalized(clean %>% select(all_of(trait_cols_pca)), method = normalization_method)
  if (!check_res$all_normal && force_normalize_if_not) {
    clean <- normalize_df(clean, trait_cols_pca, method = normalization_method)
    trait_candidates <- clean %>% select(all_of(trait_cols_pca))
  }
}

# Robust aggregation: use md_groups factors (ensures proper grouping)
df_means <- md_groups %>%
  group_by(ethanol_pre_treatment, saline_treatment) %>%
  summarise(across(all_of(trait_cols_pca), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(
    ethanol_val = case_when(as.character(ethanol_pre_treatment) %in% c("None", "none") ~ 0,
                            as.character(ethanol_pre_treatment) %in% c("Ethanol", "ethanol") ~ 20,
                            TRUE ~ as.numeric(as.character(ethanol_pre_treatment))),
    saline_val = case_when(as.character(saline_treatment) %in% c("Control", "control") ~ 0,
                           as.character(saline_treatment) %in% c("S1", "s1") ~ 50,
                           as.character(saline_treatment) %in% c("S2", "s2") ~ 100,
                           TRUE ~ as.numeric(as.character(saline_treatment)))
  ) %>%
  mutate(
    treatment = paste0("E", ethanol_val, "_S", saline_val),
    treatment_label = treatment
  ) %>%
  relocate(treatment, treatment_label)

# Map to nicer labels if desired
treatment_map <- c(
  "E0_S0"   = "Control",
  "E0_S50"  = "S1",
  "E0_S100" = "S2",
  "E20_S0"  = "Eth",
  "E20_S50" = "S1 + Eth",
  "E20_S100"= "S2 + Eth"
)
df_means <- df_means %>% mutate(treatment_label = ifelse(treatment %in% names(treatment_map), treatment_map[treatment], treatment_label))

# Build matrix (rows = treatments, cols = traits) and set rownames
mat <- df_means %>% select(all_of(trait_cols_pca)) %>% as.matrix()
rownames(mat) <- df_means$treatment_label

# Check number of non-constant traits
sds <- apply(mat, 2, sd, na.rm = TRUE)
zero_var <- names(sds)[is.na(sds) | (sds <= 1e-12)]
nonconst_traits <- setdiff(colnames(mat), zero_var)

if (length(nonconst_traits) < 2) {
  warning("Not enough non-constant traits for PCA (need >= 2). Skipping PCA; heatmap will still be drawn.")
  # create heatmap-scaled matrix and draw heatmap only
  mat_scaled_heat <- scale(mat, center = TRUE, scale = TRUE)
  mat_scaled_heat[is.na(mat_scaled_heat)] <- 0
  
  # column clustering (may still work with 1 col)
  hc_cols <- try(hclust(dist(t(mat_scaled_heat), method = "euclidean"), method = "complete"), silent = TRUE)
  
  # prepare heatmap df
  col_order <- colnames(mat_scaled_heat)
  heatmap_df <- as.data.frame(mat_scaled_heat) %>%
    rownames_to_column("treatment") %>%
    pivot_longer(-treatment, names_to = "trait", values_to = "value") %>%
    mutate(trait = factor(trait, levels = col_order), treatment = factor(treatment, levels = rev(rownames(mat_scaled_heat))))
  
  heatmap_plot <- ggplot(heatmap_df, aes(x = trait, y = treatment, fill = value)) +
    geom_tile(color = "grey30") +
    scale_fill_gradient2(low = "#00B200", mid = "black", high = "#D7191C", midpoint = 0, name = "z-score") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title = element_blank())
  
  print(heatmap_plot)
  
} else {
  # proceed with full pipeline (heatmap + column dendrogram + PCA biplot)
  mat_scaled_heat <- scale(mat, center = TRUE, scale = TRUE)
  mat_scaled_heat[is.na(mat_scaled_heat)] <- 0
  
  mat_pca <- mat[, nonconst_traits, drop = FALSE]
  mat_scaled_pca <- scale(mat_pca, center = TRUE, scale = TRUE)
  mat_scaled_pca[is.na(mat_scaled_pca)] <- 0
  
  # clustering for columns (heatmap)
  hc_cols <- hclust(dist(t(mat_scaled_heat), method = "euclidean"), method = "complete")
  cluster_col_order <- hc_cols$labels[hc_cols$order]
  
  # desired labels (optional mapping by position, adjust if your trait order differs)
  desired_trait_labels <- c("LFW", "TW", "LDW", "RWC", "WD", "RFW", "PH", "RL",
                            "LN", "LA", "RDW", "RSR", "SL", "SFW", "SDW", "SD")
  
  if (length(desired_trait_labels) == ncol(mat_scaled_heat)) {
    colnames(mat_scaled_heat) <- desired_trait_labels
    pca_pos <- match(colnames(mat_pca), colnames(mat))
    colnames(mat_scaled_pca) <- desired_trait_labels[pca_pos]
    col_order <- desired_trait_labels
  } else {
    col_order <- cluster_col_order
    col_order <- intersect(col_order, colnames(mat_scaled_heat))
    if (length(col_order) != ncol(mat_scaled_heat)) {
      missing_cols <- setdiff(colnames(mat_scaled_heat), col_order)
      col_order <- c(col_order, missing_cols)
    }
  }
  cor(clean %>% select(all_of(present_traits)),
      use = "pairwise.complete.obs")
  X_complete <- clean %>%
    dplyr::select(dplyr::all_of(present_traits)) %>%
    tidyr::drop_na()
  
  cor_results <- cor(X_complete, use = "complete.obs")
  # heatmap df and plot
  heatmap_df <- as.data.frame(mat_scaled_heat) %>%
    rownames_to_column("treatment") %>%
    pivot_longer(-treatment, names_to = "trait", values_to = "value") %>%
    mutate(trait = factor(trait, levels = col_order), treatment = factor(treatment, levels = rev(rownames(mat_scaled_heat))))
  
  heatmap_plot <- ggplot(heatmap_df, aes(x = trait, y = treatment, fill = value)) +
    geom_tile(color = "grey30") +
    scale_fill_gradient2(low = "#00B200", mid = "black", high = "#D7191C", midpoint = 0,
                         limits = c(min(heatmap_df$value, na.rm = TRUE), max(heatmap_df$value, na.rm = TRUE)),
                         name = "z-score") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), axis.title = element_blank(), panel.grid = element_blank())
  
  # dendrogram (top)
  ddata <- dendro_data(hc_cols, type = "rectangle")
  segments_df <- segment(ddata)
  dendro_plot_top <- ggplot() +
    geom_segment(data = segments_df, aes(x = x, y = y, xend = xend, yend = yend), size = 0.6) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, length(col_order) + 0.6)) +
    theme_void()
  
  # PCA
  pca <- prcomp(mat_scaled_pca, center = FALSE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, 1:2], check.names = FALSE) %>% rownames_to_column("treatment_label")
  colnames(scores)[2:3] <- c("PC1", "PC2")
  # attach grouping info for plotting (if needed)
  scores <- scores %>% left_join(df_means %>% select(treatment_label, ethanol_pre_treatment, saline_treatment), by = "treatment_label")
  
  loadings <- as.data.frame(pca$rotation[, 1:2], check.names = FALSE) %>% rownames_to_column("trait")
  colnames(loadings)[2:3] <- c("PC1", "PC2")
  mult <- min(
    (max(scores$PC1) - min(scores$PC1)) / (max(loadings$PC1) - min(loadings$PC1) + 1e-9),
    (max(scores$PC2) - min(scores$PC2)) / (max(loadings$PC2) - min(loadings$PC2) + 1e-9)
  ) * 0.6
  loadings <- loadings %>% mutate(PC1s = PC1 * mult, PC2s = PC2 * mult, len = sqrt(PC1^2 + PC2^2))
  loadings_label <- loadings %>% arrange(desc(len)) %>% slice(1:min(18, nrow(loadings)))
  
  pct_var <- function(i) round(100 * (pca$sdev[i]^2 / sum(pca$sdev^2)), 2)
  xlab <- paste0("PC1 (", pct_var(1), "%)")
  ylab <- paste0("PC2 (", pct_var(2), "%)")
  
  pca_plot <- ggplot() +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_point(data = scores, aes(PC1, PC2), color = "red", size = 3) +
    geom_text_repel(data = scores, aes(PC1, PC2, label = treatment_label), color = "red", size = 3.2, max.overlaps = Inf) +
    geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1s, yend = PC2s), arrow = arrow(length = unit(0.02, "npc")), color = "blue", alpha = 0.6) +
    geom_text_repel(data = loadings_label, aes(x = PC1s, y = PC2s, label = trait), size = 3) +
    stat_ellipse(data = scores, aes(x = PC1, y = PC2), linetype = 2, color = "darkred", level = 0.68, inherit.aes = FALSE) +
    theme_minimal() + labs(x = xlab, y = ylab) +
    theme_set(theme_bw()) + #To make the gray background white
    theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #Eliminate the grids
  
  
  combined <- (dendro_plot_top / heatmap_plot) / pca_plot + plot_layout(heights = c(0.7, 4.5, 4.0))
  print(combined)
  # ggsave("combined_heatmap_pca_top_dendro.png", combined, width = 14, height = 12, dpi = 300)
}

# End of Main.R
