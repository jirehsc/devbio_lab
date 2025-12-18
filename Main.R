
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
  "root_to_shoot_ratio"
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
shapiro.test(EZ$health_status)

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


