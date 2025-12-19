
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

# Homogeneity test per treatment

fligner.test(list(CX$leaf_area, EX$leaf_area, CY$leaf_area, EY$leaf_area, CZ$leaf_area, EZ$leaf_area))
fligner.test(list(CX$leaf_number, EX$leaf_number, CY$leaf_number, EY$leaf_number, CZ$leaf_number, EZ$leaf_number))
fligner.test(list(CX$shoot_dry_weight, EX$shoot_dry_weight, CY$shoot_dry_weight, EY$shoot_fresh_weight, CZ$shoot_dry_weight, EZ$shoot_dry_weight))
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

# Descriptive analysis per treatment

# Minimum data per treatment
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

