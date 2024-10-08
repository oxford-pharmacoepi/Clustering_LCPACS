library(dplyr)
library(here)
library(pbatR)
library(lubridate)
source(here("Study/Functions.R"))

dir_data     <- "D:/Projects/GeneticDeterminants_LC_PACS/"
dir_results  <- "D:/Projects/GeneticDeterminants_LC_PACS/Results/"

# Load datasets -----
hes <- loadHesData(dir_data)
ukb <- loadUKBData(dir_data)

# Create the cohorts -----
n_symptoms <- 1
washout_period <- 30

source(here("Study/1-CreateCohorts.R"))
