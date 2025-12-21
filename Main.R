
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

clean_outliers <- clean %>%
  mutate(across(all_of(trait_cols), ~ {
    x <- .x[!is.na(.x)]
    if(length(x) < 4) return(FALSE)  # Skip if too few data
    Q1 <- quantile(x, 0.25)
    Q3 <- quantile(x, 0.75)
    IQR_val <- Q3 - Q1
    ( .x < (Q1 - 1.5*IQR_val) | .x > (Q3 + 1.5*IQR_val) ) & !is.na(.x)
  }, .names = "{.col}_outlier"))


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


#list of traits that cannot be tested due to limited sample size: leaf fresh weight, lead dw, leaf rwc, water deficit, turgid weight
shapiro.test(CX$leaf_fresh_weight)  #cant be testesd
shapiro.test(EX$leaf_rwc)
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




