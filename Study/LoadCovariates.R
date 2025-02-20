# Baseline characteristics ----
bd <- as_tibble(read.table(paste0("UKBiobank/baselineCharacteristics.tab"), header = TRUE, sep = "\t"))

bd |>
  dplyr::select("eid" = "f.eid",
         "sex" = "f.31.0.0",
         "year_of_birth" = "f.34.0.0",
         "ethnic_background" = "f.21000.0.0",
         "body_mass_index" = "f.21001.0.0",
         "smoking_status"  = "f.20116.0.0",
         "index_of_multiple_deprivation_england" = "f.26410.0.0",
         "index_of_multiple_deprivation_scotland" = "f.26427.0.0",
         "index_of_multiple_deprivation_wales" = "f.26426.0.0") |>
  mutate("index_of_multiple_deprivation" = index_of_multiple_deprivation_england) |>
  mutate("index_of_multiple_deprivation" = if_else(is.na(index_of_multiple_deprivation),
                                                   index_of_multiple_deprivation_scotland,
                                                   index_of_multiple_deprivation)) |>
  mutate("index_of_multiple_deprivation" = if_else(is.na(index_of_multiple_deprivation),
                                                   index_of_multiple_deprivation_wales,
                                                   index_of_multiple_deprivation)) |>
  dplyr::select(-starts_with("index_of_multiple_deprivation_")) |>
  write.csv(paste0(dir_results, "/LoadCovariates_baselineCharacteristics.csv"))

# Biomarkers ----
bd |>
  dplyr::select("eid" = "f.eid",
         "alanine_aminotransferase" = "f.30620.0.0",
         "albumin" = "f.30600.0.0",
         "alkaline_phosphatase" = "f.30610.0.0",
         "apolipoprotein_a" = "f.30630.0.0",
         "apolipoprotein_b" = "f.30640.0.0",
         "aspartate_aminotransferase" = "f.30650.0.0",
         "c_reactive_protein" = "f.30710.0.0",
         "calcium" = "f.30680.0.0",
         "cholesterol" = "f.30690.0.0",
         "creatinine" = "f.30700.0.0",
         "cystatin_c" = "f.30720.0.0",
         "direct_bilirubin" = "f.30660.0.0",
         "gamma_glutamyltransferase" = "f.30730.0.0",
         "glucose" = "f.30740.0.0",
         "hbauc" = "f.30750.0.0",
         "hdl_cholesterol" = "f.30760.0.0",
         "igf" = "f.30770.0.0",
         "ldl_direct" = "f.30780.0.0",
         "lipoprotein_a" = "f.30790.0.0",
         "oestradiol" = "f.30800.0.0",
         "phosphate" = "f.30810.0.0",
         "rheumatoid_factor" = "f.30820.0.0",
         "shbg" = "f.30830.0.0",
         "testosterone" = "f.30850.0.0",
         "total_bilirubin" = "f.30840.0.0",
         "total_protein" = "f.30860.0.0",
         "triglycerides" = "f.30870.0.0",
         "urate" = "f.30880.0.0",
         "urea" = "f.30670.0.0",
         "vitamin_d" = "f.30890.0.0") |>
  write.csv(paste0(dir_results,"/LoadCovariates_biomarkers.csv"))

# Comorbidities ----
hes_data <- loadHesData() # load hes data
codes <- as_tibble(read.csv(paste0("UKBiobank/Read_code/codes.csv"))) # Read codes

comorbidities <- codes$phenotype |> unique()

cohorts <- tibble("eid" = as.integer(),
                  "episode_date" = as.Date(x = integer(0), origin = "2000-01-01"),
                  "value" = as.character())

for(ci in comorbidities){
  codes_i <- codes |>
    filter(phenotype == ci)

  phenotype_i <- codes_i |> pull(phenotype) |> unique()
  codes_i     <- codes_i |> pull(icd10code)

  cohorts <- cohorts |>
    union_all(
      hes_data |>
        filter(diag_icd10 %in% codes_i) |>
        group_by(eid) |>
        dplyr::select("eid","episode_date") |>
        filter(episode_date == min(episode_date)) |>
        distinct() |>
        ungroup() |>
        mutate(value = phenotype_i)
    )
}

cohorts |>
  write.csv(paste0(dir_results,"/LoadCovariates_commorbidities.csv"))

rm(list = c("bd", "codes", "cohorts", "hes_data", "ci", "codes_i", "comorbidities",
            "phenotype_i"))
