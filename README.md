# Clustering post-acute COVID-19 patients

This repository contains code for a study on clustering individuals with post-COVID condition (PCC) from the UK Biobank COVID-19 survey and performing genome-wide association studies (GWAS). The code is structured into several parts, each of which corresponds to a different stage of the project.

Important:\
- The R scripts in this repository are not meant to be run one after the other without interruption.\
- Some steps depend on external processing (e.g., analyses run on the UK Biobank Research Analysis Platform, FUMA web outputs).\
- Clinical and methodological decisions are involved at certain stages and should be carefully reviewed before proceeding.

------------------------------------------------------------------------

## Repository structure

### Top-level scripts

-   CodeToRun.R – main driver script showing the overall workflow.
-   Study/Functions.R – helper functions used across analyses.
-   Study/LoadCovariates.R – loads baseline covariates (e.g., age, sex, etc.).
-   Study/DataCleaning.R – data preprocessing for PCC and controls.
-   Study/1-CreateCohorts.R – creates PCC and control cohorts.
-   Study/2-Clustering.R – runs latent class analysis (LCA) clustering of symptoms and characterisation of clusters.
-   Study/3-GWASFiles.R – prepares input files for GWAS (to be uploaded and run on UK Biobank RAP).
-   Study/4-PostProcessing.R – processes GWAS results after UKBB RAP runs, adds rsIDs, generates summary statistics plots.
-   Study/5-PostFUMA.R – processes and interprets FUMA outputs, annotates traits using GWAS catalog information.

## Required data and inputs

To run this project, you will need:

-   UK Biobank data (via approved application):
    -   COVID-19 survey responses
    -   Linked HES (Hospital Episode Statistics) data
    -   Genotype data (for GWAS on RAP)
-   Reference files:
    -   VCF file for rsID mapping: `GCF_000001405.40.gz` (GRCh38 reference).
-   GWAS results:
    -   Summary statistics generated on the UK Biobank Research Analysis Platform.
-   FUMA outputs:
    -   For each GWAS, FUMA folders (`snps/`, `leadSNPs/`, `genes/`, `annot/`) exported into the results directory.
-   GWAS Catalog traits CSV (for trait enrichment analysis).

------------------------------------------------------------------------

## Workflow Overview

1.  Prepare cohorts and clustering (Parts 1–2)
    -   Run covariate loading, cleaning, cohort creation, and clustering scripts.
    -   Inspect cluster characteristics in generated tables/figures.
2.  Prepare GWAS cohorts
    -   Use `3-GWASFiles.R` to generate input datasets.
    -   Upload to UKBB RAP for GWAS analysis (not done in R here).
3.  Post-processing GWAS results (Part 2)
    -   Download summary statistics from RAP.
    -   Run `4-PostProcessing.R` to annotate, QC, plot results and prepare data for FUMA analysis (not done in R here).
4.  FUMA analyses
    -   Upload GWAS results to [FUMA](https://fuma.ctglab.nl/).
    -   Export FUMA outputs and place them in `Results/GWAS/`.
    -   Run `5-PostFUMA.R` for trait enrichment and interpretation.

------------------------------------------------------------------------

## Notes

-   This repository contains code only.
-   No UKBB data, results, or sensitive outputs are included.
-   Users must have their own UK Biobank approval to access and analyse the data.
-   Clinical and methodological decisions (e.g., PCC definitions, cluster number selection, GWAS model choices) should be carefully considered before replicating analyses.
