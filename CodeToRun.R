#################################################################################
# Clustering individuals with and without PACS/LC from the UKBB COVID-19 survey #
#################################################################################

# Libraries needed for this study
library(dplyr)
library(here)
library(pbatR)
library(lubridate)
library(readr)
library(patchwork)
library(ggplot2)
library(mice)
library(flextable)
library(stringr)
library(ftExtra)
library(officer)
library(ggplot2)
library(poLCA)
library(IsingFit)
library(igraph)
library(tidyr)
source(here("Study/Functions.R"))

dir_data     <- "D:/Projects/Clustering_LCPACS/"
dir_results  <- "D:/Projects/Clustering_LCPACS/Results"

# Load datasets -----
hes <- loadHesData()
ukb <- loadUKBData(dir_data)

# Create the cohorts and do clustering
n_symptoms <- 1
washout_period <- 30

dir.create(dir_results, recursive = TRUE)
dir.create(paste0(dir_results,"/Tables"), recursive = TRUE)
dir.create(paste0(dir_results,"/GWAS"), recursive = TRUE)

source(here::here("Study/LoadCovariates.R"))
source(here::here("Study/DataCleaning.R"))
source(here("Study/1-CreateCohorts.R"))
source(here("Study/2-Clustering.R"))
