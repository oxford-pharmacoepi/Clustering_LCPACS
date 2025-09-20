# Load data ----
health_questionnaire <- loadHealthAndWellBeingQuestionnaire()
covid19_result <- loadCovid19Result()

# Create initial cohort ----
cohort <- covid19_result %>%
  dplyr::select(eid, specdate, result) %>%
  recordAttrition() %>%
  filter(result == 1) %>%
  dplyr::select(-"result") %>%
  distinct() %>%
  recordAttrition("Restrict to participants with a positive COVID-19 test result") %>%
  inner_join(
    health_questionnaire %>%
      filter(!is.na(questionnaire_started)),
    by = "eid"
  ) %>%
  recordAttrition("Restrict to participants that answered the health questionnaire") %>%
  filter(specdate < questionnaire_started) %>%
  recordAttrition("Restrict to participants with positive test before answering questionnaire")

cohort_latest <- cohort %>%
  group_by(eid) %>%
  filter(specdate == max(specdate)) %>%
  ungroup() %>%
  restoreAttrition(cohort) %>%
  recordAttrition("Restrict to the latest record.")

cohort_filtered <- cohort_latest %>%
  filter(
    ((specdate + washout_period) < questionnaire_started) &
      (specdate > (questionnaire_started - 365))
  ) %>%
  recordAttrition(paste0("Restrict to participants with infections within 1 year and ", washout_period, " days"))

# Identify symptom columns
symptom_cols <- gsub("symptom_", "", colnames(cohort_filtered %>% dplyr::select(starts_with("symptom"))))

# Process symptom variables
for (symptom in symptom_cols) {
  var <- paste0("symptom_", symptom)
  var_len <- paste0("length_", symptom)
  
  cohort_filtered <- cohort_filtered %>%
    mutate(
      !!sym(var) := if_else(!!sym(var) < 0, NA, !!sym(var)),
      !!sym(var) := if_else(
        questionnaire_started - !!sym(var_len) < specdate &
          !is.na(!!sym(var_len)) & !is.na(!!sym(var)),
        0, !!sym(var)
      ),
      !!sym(var_len) := if_else(!!sym(var) == 0, NA, !!sym(var_len))
    )
}

cohort_aggregated <- cohort_filtered %>%
  rowwise() %>%
  mutate(symptom = sum(c_across(starts_with("symptom")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(length = do.call(pmin, c(dplyr::select(., starts_with("length")), na.rm = TRUE))) %>%
  dplyr::select(-starts_with("length_"))

longCovid_cases <- cohort_aggregated %>%
  restoreAttrition(cohort_filtered) %>%
  filter((questionnaire_started - length) > specdate) %>%
  recordAttrition("Restrict to participants with at least one non-persistent symptom") %>%
  filter(symptom >= n_symptoms) %>%
  recordAttrition("Restrict to participants with multiple symptoms")

longCovid_controls <- cohort_aggregated %>%
  restoreAttrition(cohort_filtered) %>%
  filter(symptom < n_symptoms & symptom >= 0) %>%
  recordAttrition("Restrict to participants without multiple symptoms")

# longCovid_cohort <- longCovid_cases %>%
#   dplyr::select(eid, specdate, questionnaire_started) %>%
#   mutate(state = 1) %>%
#   union_all(
#     longCovid_controls %>%
#       dplyr::select(eid, specdate, questionnaire_started) %>%
#       mutate(state = 0)
#   )

set.seed(11)
data_clustering <- longCovid_cases %>%
  dplyr::select(-c(questionnaire_started, symptom, length, specdate))

data_clustering[, 2:ncol(data_clustering)] <- data_clustering[, 2:ncol(data_clustering)] + 1  # Adjust for poLCA
  
x_vars <- c("1")
names_symptoms <- colnames(longCovid_cases %>%
                             dplyr::select(-c(eid, specdate, questionnaire_started, symptom, length)))
cols <- paste0("cbind(", paste(names_symptoms, collapse = ","), ")")
f <- with(data_clustering, as.formula(sprintf("%s ~ %s", cols, paste(x_vars, collapse = " + "))))

lc <- poLCA(f, data_clustering, nclass=3, maxiter=5000, graphs = FALSE,
            tol=1e-5, na.rm=FALSE,
            nrep=10, verbose=TRUE, calc.se=FALSE)

data_clustering <- data_clustering %>%
  mutate(cluster_assignment = lc$predclass) %>%
  dplyr::left_join(
    longCovid_cases %>% 
      dplyr::select(c("eid", "specdate")),
    by = "eid"
  ) %>%
  compute()

# Add all GWAS variables to cohorts of interest
# SAVE all cohorts in a UKBB GWAS format

# Add GWAS variables
data_clustering <- data_clustering |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  dplyr::select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with the same sex and genetic sex recorded") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")

longCovid_controls <- longCovid_controls |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  dplyr::select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with the same sex and genetic sex recorded") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate") |>
  dplyr::select(-c("questionnaire_started", "symptom", "length"))

data_clustering_1vs23 <- as.phe(data_clustering |> 
                                  dplyr::mutate(
                                    state = dplyr::if_else(
                                      cluster_assignment == 1, 1, 0
                                    )
                                  ) |>
                                  dplyr::select("IID" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                  mutate("FID" = IID) |> 
                            relocate("FID"), "FID", "IID")

data_clustering_2vs13 <- as.phe(data_clustering |> 
                                  dplyr::mutate(
                                    state = dplyr::if_else(
                                      cluster_assignment == 2, 1, 0
                                    )
                                  ) |>
                                  dplyr::select("IID" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                  mutate("FID" = IID) |> 
                                  relocate("FID"), "FID", "IID")

data_clustering_3vs12 <- as.phe(data_clustering |> 
                                  dplyr::mutate(
                                    state = dplyr::if_else(
                                      cluster_assignment == 3, 1, 0
                                    )
                                  ) |>
                                  dplyr::select("IID" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                  mutate("FID" = IID) |> 
                                  relocate("FID"), "FID", "IID")

data_clustering_attrition <- attr(data_clustering, "cohort_attrition")

data_clustering_1vsno <- as.phe(data_clustering |> 
                                  dplyr::filter(cluster_assignment == 1) |>
                                  dplyr::mutate(state = 1) |>
                                  dplyr::select(-"cluster_assignment") |>
                                  union_all(
                                    longCovid_controls %>%
                                      mutate(state = 0)
                                  ) |>
                                  dplyr::select("IID" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                  mutate("FID" = IID) |> 
                                  relocate("FID"), "FID", "IID")

data_clustering_2vsno <- as.phe(data_clustering |> 
                                  dplyr::filter(cluster_assignment == 2) |>
                                  dplyr::mutate(state = 1) |>
                                  dplyr::select(-"cluster_assignment") |>
                                  union_all(
                                    longCovid_controls %>%
                                      mutate(state = 0)
                                  ) |>
                                  dplyr::select("IID" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                  mutate("FID" = IID) |> 
                                  relocate("FID"), "FID", "IID")

data_clustering_3vsno <- as.phe(data_clustering |> 
                                  dplyr::filter(cluster_assignment == 3) |>
                                  dplyr::mutate(state = 1) |>
                                  dplyr::select(-"cluster_assignment") |>
                                  union_all(
                                    longCovid_controls %>%
                                      mutate(state = 0)
                                  ) |>
                                  dplyr::select("IID" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                  mutate("FID" = IID) |> 
                                  relocate("FID"), "FID", "IID")

data_clustering_1vsall <- as.phe(data_clustering_1vs23 |> 
                                  union_all(
                                    longCovid_controls %>%
                                      dplyr::mutate(state = 0) |>
                                      dplyr::select("pid" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                      dplyr::mutate("id" = pid) |> 
                                      dplyr::relocate("pid")
                                  ), "pid", "id")

data_clustering_2vsall <- as.phe(data_clustering_2vs13 |> 
                                   union_all(
                                     longCovid_controls %>%
                                       dplyr::mutate(state = 0) |>
                                       dplyr::select("pid" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                       dplyr::mutate("id" = pid) |> 
                                       dplyr::relocate("pid")
                                   ), "pid", "id")

data_clustering_3vsall <- as.phe(data_clustering_3vs12 |> 
                                   union_all(
                                     longCovid_controls %>%
                                       dplyr::mutate(state = 0) |>
                                       dplyr::select("pid" = "eid", "state", "BC_age" = "age_when_infected", "BC_sex" = "sex", "BC_age2" = "age2", "BC_agesex" = "agesex", "batch", starts_with("pc")) |> 
                                       dplyr::mutate("id" = pid) |> 
                                       dplyr::relocate("pid")
                                   ), "pid", "id")

controls_attrition <- attr(longCovid_controls, "cohort_attrition")

save(data_clustering_1vs23,
     data_clustering_2vs13,
     data_clustering_3vs12,
     data_clustering_1vsno,
     data_clustering_2vsno,
     data_clustering_3vsno,
     data_clustering_1vsall,
     data_clustering_2vsall,
     data_clustering_3vsall,
     data_clustering_attrition,
     controls_attrition,
     file = paste0(dir_clust_results,"/GWAS/","cohortsData.Rdata"))

write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_1vs23.phe"), data_clustering_1vs23)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_2vs13.phe"), data_clustering_2vs13)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_3vs12.phe"), data_clustering_3vs12)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_1vsno.phe"), data_clustering_1vsno)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_2vsno.phe"), data_clustering_2vsno)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_3vsno.phe"), data_clustering_3vsno)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_1vsall.phe"), data_clustering_1vsall)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_2vsall.phe"), data_clustering_2vsall)
write.phe(paste0(dir_clust_results,"/GWAS/","data_clustering_3vsall.phe"), data_clustering_3vsall)

# Get input in format necessary for UKBB WGS GWAS pipeline
write.table(data_clustering_1vs23 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_1vs23.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vs13 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_2vs13.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vs12 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_3vs12.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_1vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_1vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_2vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_3vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_1vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_1vsall.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_2vsall.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id"), paste0(dir_clust_results,"/GWAS/","covariate_3vsall.txt"), row.names = FALSE, quote = FALSE)

write.table(data_clustering_1vs23 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_1vs23.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vs13 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_2vs13.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vs12 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_3vs12.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_1vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_1vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_2vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_3vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_1vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_1vsall.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_2vsall.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID")), paste0(dir_clust_results,"/GWAS/","sample_included_3vsall.txt"), row.names = FALSE, quote = FALSE)

write.table(data_clustering_1vs23 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_1vs23.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vs13 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_2vs13.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vs12 |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_3vs12.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_1vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_1vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_2vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vsno |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_3vsno.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_1vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_1vsall.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_2vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_2vsall.txt"), row.names = FALSE, quote = FALSE)
write.table(data_clustering_3vsall |>
              dplyr::rename("FID" = "pid",
                            "IID" = "id") |> dplyr::select(c("FID", "IID", "state")), paste0(dir_clust_results,"/GWAS/","phenotype_3vsall.txt"), row.names = FALSE, quote = FALSE)

# Characteristics of all cohorts
dir_gwas <- "/Users/spet5356/Documents/GitHub/Clustering_LCPACS/GWAS"

# Load characteristics data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_clust_results,"/GWAS","/CleanData_baselineCharacteristics.csv"))) |> dplyr::select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_clust_results,"/GWAS","/CleanData_biomarkers.csv")))  |> dplyr::select(-c("X"))
comorbidities <- as_tibble(read.csv(paste0(dir_clust_results,"/GWAS","/CleanData_comorbidities.csv")))  |> dplyr::select(-c("X"))

# List of six datasets
datasets <- list(
  "1vs23" = data_clustering_1vs23,
  "2vs13" = data_clustering_2vs13,
  "3vs12" = data_clustering_3vs12,
  "1vsno" = data_clustering_1vsno,
  "2vsno" = data_clustering_2vsno,
  "3vsno" = data_clustering_3vsno,
  "1vsall" = data_clustering_1vsall,
  "2vsall" = data_clustering_2vsall,
  "3vsall" = data_clustering_3vsall
)

# Matching cohort names for each dataset
# (controls always first, cases second)
cohort_names <- list(
  "1vs23" = c("Controls (clust2+3)", "Cases (clust1)"),
  "2vs13" = c("Controls (clust1+3)", "Cases (clust2)"),
  "3vs12" = c("Controls (clust1+2)", "Cases (clust3)"),
  "1vsno" = c("Controls (COVID)", "Cases (clust1)"),
  "2vsno" = c("Controls (COVID)", "Cases (clust2)"),
  "3vsno" = c("Controls (COVID)", "Cases (clust3)"),
  "1vsall" = c("Controls (clust2+3,COVID)", "Cases (clust1)"),
  "2vsall" = c("Controls (clust1+3,COVID)", "Cases (clust2)"),
  "3vsall" = c("Controls (clust1+2,COVID)", "Cases (clust3)")
)

# A wrapper function for whole pipeline ----
make_characterisation_docx <- function(dat, cohorts, suffix, dir_clust_results) {
  name <- suffix
  name_cohort <- cohorts
  
  # Step 1: table creation
  longCovid_clust <- tableOneStep1(dat |> 
                                     dplyr::rename("eid" = "id") |>
                                     dplyr::left_join(cohort_aggregated |>
                                                        dplyr::select(c("eid", "specdate"))), 
                                   baselineCharacteristics, biomarkers, name, name_cohort)
  
  # Step 2: cleaning + ordering
  y <- longCovid_clust |>
    mutate(order = row_number()) |>
    dplyr::select(-c("order")) |>
    mutate(`Risk factor` = if_else(`Risk factor` == "Sociodemographics factors", 
                                   "Sociodemographic factors", 
                                   `Risk factor`))
  y <- y |>
    filter(`Risk factor` == "N") |>
    rbind(
      y |>
        filter(`Risk factor` != "N")
    )
  
  # Step 3: formatting flextable
  merge <- c(1:nrow(y))[y$`Risk factor` %in% c("Sociodemographic factors", 
                                               "Comorbidities [Cases (%)]", 
                                               "Biomarkers [Mean (SD)]")]
  
  row_fill <- c(1:nrow(y))[!y$`Long COVID clustering_Controls` == " "]
  
  ft <- y |>
    flextable() |>
    span_header() |>
    bold(i = 1, part = "header") |>
    bold(i = 2, part = "header") |>
    align(j = c(2,3), align = "center", part = "all") |>
    width(j = 1, width = 5, unit = "cm") |>
    bg(bg = "#F2F2F2", i = 1, j = c(2,3), part = "header") 
  
  for(ii in merge){
    ft <- ft |>
      merge_at(i = ii, j = c(1:3)) |>
      hline(i = ii-1, border = fp_border(color = "black")) |>
      bold(i = ii)
  }
  
  # Step 4: export
  save_as_docx(ft, path = paste0(dir_clust_results, "/Tables/BaselineCharacteristics_clustering_", suffix, ".docx"))
}

# Run it for all datasets ----
for (nm in names(datasets)) {
  make_characterisation_docx(
    dat = datasets[[nm]],
    cohorts = cohort_names[[nm]],
    suffix = nm,
    dir_clust_results = dir_gwas
  )
}

# Create a table one with less data with all cohort together
charac_tables <- list()
make_characterisation <- function(dat, cohorts, suffix) {
  name <- suffix
  name_cohort <- cohorts
  longCovid_clust <- tableOneStep1(dat |> 
                                     dplyr::rename("eid" = "id") |>
                                     dplyr::left_join(cohort_aggregated |>
                                                        dplyr::select(c("eid", "specdate"))), 
                                   baselineCharacteristics, biomarkers, name, name_cohort)
  
  y <- longCovid_clust |>
    mutate(order = row_number()) |>
    dplyr::select(-c("order")) |>
    mutate(`Risk factor` = if_else(`Risk factor` == "Sociodemographics factors", 
                                   "Sociodemographic factors", 
                                   `Risk factor`))
  y <- y |>
    filter(`Risk factor` == "N") |>
    rbind(
      y |>
        filter(`Risk factor` != "N")
    )
  
  # # Step 3: formatting flextable
  # merge <- c(1:nrow(y))[y$`Risk factor` %in% c("Sociodemographic factors", 
  #                                              "Comorbidities [Cases (%)]", 
  #                                              "Biomarkers [Mean (SD)]")]
  # 
  # row_fill <- c(1:nrow(y))[!y$`Long COVID clustering_Controls` == " "]
  # 
  # ft <- y |>
  #   flextable() |>
  #   span_header() |>
  #   bold(i = 1, part = "header") |>
  #   bold(i = 2, part = "header") |>
  #   align(j = c(2,3), align = "center", part = "all") |>
  #   width(j = 1, width = 5, unit = "cm") |>
  #   bg(bg = "#F2F2F2", i = 1, j = c(2,3), part = "header") 
  # 
  # for(ii in merge){
  #   ft <- ft |>
  #     merge_at(i = ii, j = c(1:3)) |>
  #     hline(i = ii-1, border = fp_border(color = "black")) |>
  #     bold(i = ii)
  # }
  # 
  # ft
  
  y <- y |> mutate(GWAS = suffix)
  return(y)
}


# Run it for all datasets ----
for (nm in names(datasets)) {
  charac_tables[[nm]] <- make_characterisation(
    dat = datasets[[nm]],
    cohorts = cohort_names[[nm]],
    suffix = nm
  )
}

# Start with the first GWAS tibble
table_combined <- charac_tables[[1]] |> dplyr::select(-"GWAS")

# Loop through the rest and join by "Risk factor"
for (i in 2:length(charac_tables)) {
  table_combined <- full_join(
    table_combined,
    charac_tables[[i]] |> dplyr::select(-"GWAS"),
    by = "Risk factor",
    suffix = c("", paste0("_", names(charac_tables)[i]))
  )
}

ft_full <- table_combined %>%
  flextable() %>%
  bold(part = "header") %>%
  align(j = 2:ncol(table_combined), align = "center", part = "all") %>%
  autofit()

# Identify main section rows to merge and bold
merge <- which(table_combined$`Risk factor` %in% c(
  "Sociodemographic factors", 
  "Comorbidities [Cases (%)]", 
  "Biomarkers [Mean (SD)]"
))

# Identify non-empty rows for shading
row_fill <- which(!apply(table_combined[, -1], 1, function(x) all(x == " ")))

# Create flextable
ft <- table_combined %>%
  flextable() %>%
  span_header() %>%
  bold(part = "header") %>%
  align(j = 2:ncol(table_combined), align = "center", part = "all") %>%
  width(j = 1, width = 5, unit = "cm") %>%
  bg(bg = "#F2F2F2", i = 1, j = 2:ncol(table_combined), part = "header")

# Merge and bold main sections
for(ii in merge){
  ft <- ft %>%
    merge_at(i = ii, j = 1:ncol(table_combined)) %>%
    hline(i = ii-1, border = fp_border(color = "black")) %>%
    bold(i = ii)
}

# shade rows with data
ft <- ft %>%
  bg(i = row_fill, bg = "#F9F9F9")

# Preview
ft

# Export full table
save_as_docx(ft, path = paste0(dir_clust_results, "/Tables/", "Table1_all_rows.docx"))

# Slice first 15 rows
table_first15 <- table_combined %>%
  slice(1:15)

merge <- which(table_combined$`Risk factor` %in% c(
  "Sociodemographic factors"
))

row_fill <- which(!apply(table_first15[, -1], 1, function(x) all(x == " ")))

colnames(table_first15)[colnames(table_first15) == "Controls (clust1+3,COVID)"] <- "Controls\n(clust1+3,\nCOVID)"
colnames(table_first15)[colnames(table_first15) == "Controls (clust2+3,COVID)"] <- "Controls\n(clust2+3,\nCOVID)"
colnames(table_first15)[colnames(table_first15) == "Controls (clust1+2,COVID)"] <- "Controls\n(clust1+2,\nCOVID)"

# Create flextable
ft <- table_first15 %>%
  flextable() %>%
  span_header() %>%
  bold(part = "header") %>%
  align(j = 2:ncol(table_first15), align = "center", part = "all") %>%
  bg(bg = "#F2F2F2", i = 1, j = 2:ncol(table_first15), part = "header")

# Merge and bold main sections
for(ii in merge){
  ft <- ft %>%
    merge_at(i = ii, j = 1:ncol(table_first15)) %>%
    hline(i = ii-1, border = fp_border(color = "black")) %>%
    bold(i = ii)
}

# shade rows with data
ft <- ft %>%
  bg(i = row_fill, bg = "#F9F9F9")

# Preview
ft

# Export full table
save_as_docx(ft, path = paste0(dir_clust_results, "/Tables/", "Table1_15_rows.docx"))

# As LaTeX not formatted
table_combined %>%
  slice(1:15) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE,
        caption = "Table 1. First 15 baseline characteristics of all GWAS cohorts") %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  save_kable(file = paste0(dir_clust_results, "/Tables/", "Table1_15_rows.tex"))

