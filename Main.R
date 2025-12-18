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

summary(metadata)

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
  

# Verify/recheck
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

# List of all traits
trait_cols <- c(
  "leaf_area", "leaf_number", "shoot_dry_weight", "root_dry_weight", 
  "shoot_fresh_weight", "root_fresh_weight", "leaf_fresh_weight", 
  "leaf_dry_weight", "turgid_weight", "leaf_rwc", "water_deficit", 
  "plant_height", "stem_diameter", "shoot_length", "root_length", 
  "root_to_shoot_ratio "
)

# Apply per trait and see outliers each trait
trait_outliers <- lapply(trait_cols, function(tr) calculate_outliers(clean, tr))
names(trait_outliers) <- trait_cols

# FLAG outliers inside df (if viewed, additional coloumns, true=outliers, false=normal)
clean_outliers <- clean
for (tr in trait_cols) {
  out_name <- paste0(tr, "_outlier")
  clean_outliers[[out_name]] <- trait_outliers[[tr]]
}


# visualize trait
boxplot(clean$root_fresh_weight)

