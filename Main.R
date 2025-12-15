#Loading data

md <- read.csv("metadata.csv", header = TRUE)
names(md)

summary(metadata)

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
library(ggp
library(janitor)

# Data Pipelining and Cleaning
################################

library(ggpubr)

summary(metadata)

#njchnfcjdnhc

###### DATA CLEANING ######

#import dataset

data <- read.csv("metadata.csv", header = TRUE)
View(data)

#load library
library(readxl)
library(tidyverse)
library(janitor)

#observe/ inspect data
head(data)
tail(data)
str(data)
summary(data)
###note: NA in columns: LFW,LDW, TW, RWC, WD 
###note: dead plants in 3 rows (specimenid: CX2, CY2, EY2,)

#data types and errors

clean_data <- rw %>%
  clean_names() %>%  
  mutate(
    # Fix data types
    across(c(leaf_area, leaf_number, shoot_dry_weight, root_dry_weight, 
             shoot_fresh_weight, root_fresh_weight, leaf_fresh_weight, 
             leaf_dry_weight, turgid_weight, leaf_rwc, water_deficit, 
             plant_height, stem_diameter, shoot_length, root_length, 
             root_to_shoot_ratio, health_status, block_mortality_rate, 
             block_survival_rate), 
           ~ case_when(
             .x == "class" ~ NA_real_,
             TRUE ~ as.numeric(as.character(.x))
           )),
    
    specimen_id = as.character(specimen_id)
    
    # make block as factor
    block = as.factor(block)
  ) %>%
  
  view(clean_data)
  
  # Remove the dead 
  filter(!(leaf_area == 0 & leaf_number == 0 & shoot_dry_weight == 0 & root_dry_weight == 0)) %>%
  distinct(specimen_id, .keep_all = TRUE)  %>% 
  
  #NAs to 0
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Verify/recheck
class(clean_data$block)
class(clean_data$ethanol_pre_treatment)
print(paste("Columns:", ncol(clean_data)))
print(paste("Rows:", nrow(clean_data)))
names(clean_data)  
head(clean_data)
view(clean_data)
