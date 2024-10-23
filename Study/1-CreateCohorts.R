# Load databases ----
health_questionnaire <- loadHealthQuestionnaire(ukb)
covid19_result       <- loadCovid19Result(ukb)
sequela_table        <- loadSequelaTable()
genetic_data         <- loadGeneticData(ukb)
waves_data           <- loadWaveData(ukb)

# Create LC cohort ----
x <- covid19_result |>
  select("eid","specdate","result") |>
  recordAttrition() |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  inner_join(
    health_questionnaire |>
      filter(!is.na(questionnaire_started)),
    by = "eid"
  ) |>
  recordAttrition("Restrict to participants that answered the health and well-being web-questionnaire") |>
  filter(questionnaire_started != as.Date("1999-01-01")) |>
  recordAttrition("Restrict to participants that answered yes/no to all the questions from the health and well-being web-questionnaire.") |>
  filter(specdate < questionnaire_started) |>
  recordAttrition("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction before answering the questionnaire.")

x1 <- x |> group_by(eid) |>
  filter(specdate == max(specdate)) |>
  ungroup() |>
  restoreAttrition(x) |>
  recordAttrition("Restrict to the latest record.")

x2 <- x1 |>
  filter(((specdate+washout_period) < questionnaire_started) & (specdate > (questionnaire_started-365))) |>
  recordAttrition(paste0("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction between 1 year and ",washout_period," days before answering the questionnaire."))

x3 <- x2 |>
  rowwise() |>
  mutate("symptom" = across(starts_with("symptom")) |> sum()) |>
  ungroup() |>
  restoreAttrition(x2) %>%
  mutate(length = do.call(pmax, c(select(., starts_with("length")), na.rm = TRUE))) |>
  select(-starts_with("length_")) |>
  filter((questionnaire_started-length) > specdate) |>
  recordAttrition("Restrict to people with no persistent symptoms") |>
  inner_join(
    loadGeneticData(ukb), by = "eid"
  ) |>
  filter(!is.na(pc1)) |>
  recordAttrition("Restrict to people with available genetic data") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with equal sex and genetic sex recorded")
  

longCovid_cases_cohort <- x3 |>
  filter(symptom >= n_symptoms) |>
  recordAttrition("Restrict to people that reported at least one symptom of the questionnaire")

longCovid_controls_cohort <- x3 |>
  filter(symptom < n_symptoms & symptom >= 0) |>
  recordAttrition("Restrict to people that did not report any symptom of the questionnaire")

longCovid_cohort <- longCovid_cases_cohort |>
  select("eid", "specdate", "questionnaire_started", "year_birth") |>
  mutate(state = 1) |>
  union_all(
    longCovid_controls_cohort |>
      select("eid", "specdate", "questionnaire_started", "year_birth") |>
      mutate(state = 0)
  ) |>
  mutate("age_when_infected" = year(specdate) - year_birth) |>
  select(-c("year_birth"))

rm(list = c("health_questionnaire", "longCovid_cases_cohort", "longCovid_controls_cohort", 
            "x", "x1", "x2", "x3"))

# Create PACS cohort ----
t <- hes |>
  recordAttrition() |>
  mutate(diag_icd10   = if_else(!diag_icd10 %in% sequela_table$icd10_code, NA, diag_icd10),
         episode_date = if_else(is.na(diag_icd10), NA, episode_date)) |>
  distinct() |>
  recordAttrition("Restrict to PACS records") |>
  inner_join(
    covid19_result |>
      select("eid","specdate","result"),
    by = "eid",
    relationship = "many-to-many"
  ) |>
  recordAttrition("Restrict to participants with COVID-19 linkage data") |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  mutate(episode_date = if_else(episode_date >= (specdate-365), episode_date, NA)) |>
  mutate(episode_date = if_else(episode_date <= (specdate+365), episode_date, NA)) |>
  mutate(diag_icd10 = if_else(is.na(episode_date), NA, diag_icd10)) |>
  distinct() |>
  recordAttrition("Restrict the study period between one year before and after testing positive for COVID-19.") |>
  mutate(exclusion = if_else(
    episode_date >= (specdate-365) & episode_date <= specdate & (!is.na(diag_icd10)),
    1,0
  ))

t1 <- t |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t) |>
  filter(exclusion == 0) |>
  recordAttrition("Restrict to participants that did not have any PACS event one year prior the COVID-19 infection") |>
  mutate(exclusion = if_else((!is.na(diag_icd10)) &
                               (episode_date > specdate) &
                               (episode_date <= (specdate+washout_period)), 1, 0))
t2 <- t1 |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t1) |>
  filter(exclusion == 0) |>
  recordAttrition(paste0("Restrict to participants that did not have any PACS event after ", washout_period, " days of the infection")) |>
  mutate(diagnoses = if_else(!is.na(diag_icd10), 1, 0))

t3 <- t2 |>
  group_by(eid) |>
  mutate(cases = sum(diagnoses)) |>
  ungroup() |>
  restoreAttrition(t2)

pacs_cases_cohort <- t3 |>
  filter(cases >= n_symptoms) |>
  filter(!is.na(diag_icd10)) |>
  recordAttrition("Restrict to cases")

pacs_controls_cohort <- t3 |>
  filter(cases >= 0 & cases < n_symptoms) |>
  distinct() |>
  recordAttrition("Restrict to controls")

pacs_cohort  <- pacs_cases_cohort |>
  union_all(pacs_controls_cohort) |>
  select("eid", "specdate", "state" = "cases") |>
  group_by(eid) |>
  mutate(specdate = min(specdate),
         state = if_else(state < n_symptoms, 0, 1)) |>
  distinct() |>
  ungroup() |>
  restoreAttrition(pacs_cases_cohort)

rm(list = c("t", "t1", "t2", "t3", "pacs_cases_cohort", "pacs_controls_cohort"))

# Prepare .phe files and save cohorts -----
longCovid_attrition <- attr(longCovid_cohort, "cohort_attrition")
pacs_attrition <- attr(pacs_cohort, "cohort_attrition")

save(longCovid_cohort,
     pacs_cohort,
     longCovid_attrition,
     pacs_attrition,
     file = paste0(dir_data,"Results/cohortsData.Rdata"))


longCovid_cohort <- as.phe(longCovid_cohort |> select("IID" = "eid", "state") |> mutate("FID" = IID) |> relocate("FID"), "FID", "IID")
pacs_cohort      <- as.phe(pacs_cohort |> select("IID" = "eid", "state") |> mutate("FID" = IID) |> relocate("FID"), "FID", "IID")

write.phe(paste0(dir_data,"/Results/LongCovid_cohort.phe"), longCovid_cohort)
write.phe(paste0(dir_data,"/Results/Pacs_cohort.phe"), pacs_cohort)

rm(list = c("longCovid_cohort", "pacs_cohort", "longCovid_attrition", "pacs_attrition"))

# Create LC validation cohorts ----
t1 <- waves_data |>
  recordAttrition() |>
  filter(across(ends_with("_antibody_test_result"), ~!is.na(.))) |>
  recordAttrition("Participants with no missing antibody test data within the first 6 waves") |>
  filter(across(ends_with("self-reported_date_that_antibody_test_sample_was_collected"), ~ !. %in% c(as.Date("1900-01-01"), as.Date("1999-01-01")))) |>
  recordAttrition("Participants with a valid self-reported date that antibody test sample was collected") |>
  mutate("infection" = rowSums(across(ends_with("_antibody_test_result")), na.rm = TRUE)) |> 
  filter(infection > 0) |>
  recordAttrition("Participants with at least one reported COVID-19 infection") 

for(i in 0:5){
  date_col <- paste0("w",i,"_self-reported_date_that_antibody_test_sample_was_collected")
  inf_col  <- paste0("w",i,"_antibody_test_result")
  t1[,date_col][!t1[,inf_col]]<- NA
}

t1 <- t1 |>
  rowwise() |>
  mutate("infection_date" = min(c_across(ends_with("self-reported_date_that_antibody_test_sample_was_collected")), na.rm = TRUE)) |>
  ungroup() 

# Website or post


  
  

    


  
