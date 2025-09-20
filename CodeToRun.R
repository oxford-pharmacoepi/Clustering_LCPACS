#################################################################################
# Clustering individuals with PCC from the UKBB COVID-19 survey and doing GWAS #
#################################################################################
# Libraries needed for this study
library(dplyr)
library(here)
library(pbatR)
library(lubridate)
library(readr)
library(ggplot2)
library(patchwork)
library(mice)
library(flextable)
library(stringr)
library(ftExtra)
library(officer)
library(ggplot2)
library(poLCA)
library(IsingFit)
library(igraph)
library(kableExtra)
library(qqman)
library(VariantAnnotation)
library(GenomicRanges)
library(tidyr)
library(forcats)
library(wordcloud)
library(RColorBrewer)
library(data.table)
library(gwasrapidd)
source(here("Study/Functions.R"))

dir_data     <- "/Users/spet5356/Documents/GitHub/Clustering_LCPACS/UKBiobank"
dir_results  <- "/Users/spet5356/Documents/GitHub/Clustering_LCPACS/Results"
dir_clust_results  <- "/Users/spet5356/Documents/GitHub/Clustering_LCPACS/Results/Clustering"
dir_gwas_results  <- "/Users/spet5356/Documents/GitHub/Clustering_LCPACS/Results/GWAS"

# Load datasets -----
hes <- loadHesData()
ukb <- loadUKBData(dir_data)

# Create the cohorts and do clustering
n_symptoms <- 1
washout_period <- 30

dir.create(dir_clust_results, recursive = TRUE)
dir.create(dir_gwas_results, recursive = TRUE)
dir.create(paste0(dir_clust_results,"/Tables"), recursive = TRUE) # clustering characteristics
dir.create(paste0(dir_clust_results,"/GWAS"), recursive = TRUE) # input files for UKBB RAP
dir.create(paste0(dir_gwas_results,"/Tables"), recursive = TRUE) # characteristics of GWAS cohorts
dir.create(paste0(dir_gwas_results,"/GWAS"), recursive = TRUE) # MP, QQP and FUMA
dir.create(paste0(dir_gwas_results,"/Traits"), recursive = TRUE) # post FUMA, reported trait analysis

## PART 1: Use UKBB data to study clustering of PCC and prepare GWAS inputs
# Load covariates and data for the study
source(here::here("Study/LoadCovariates.R"))
source(here::here("Study/DataCleaning.R"))
# Create PCC and control cohorts
source(here("Study/1-CreateCohorts.R"))
# Perform clustering and characterisation. Needs UKBB characteristics data in dir_clust_results
# Results for this analysis can be inspected in the Shiny
source(here("Study/2-Clustering.R"))
# Create final GWAS cohorts and analyse demographics
source(here("Study/3-GWASFiles.R"))

## PART 2: After GWAS, check all results
# Post-process results. 
# Needs GWAS summary statistics inside dir_gwas_results.
# Needs also "GCF_000001405.40.gz" VCF file to add rsIDs to GRCh38 chr/pos coordinates
source(here("Study/4-PostProcessing.R"))
# Checks FUMA output and analyses traits. 
# Needs each GWAS FUMA output in its own folder (snps, leadSNPs, genes, annot) in dir_gwas_results. 
# Needs also GWAS catalog information in a csv for traits.
source(here("Study/5-PostFUMA.R"))
