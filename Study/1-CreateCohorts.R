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

longCovid_cohort <- longCovid_cases %>%
  dplyr::select(eid, specdate, questionnaire_started) %>%
  mutate(state = 1) %>%
  union_all(
    longCovid_controls %>%
      dplyr::select(eid, specdate, questionnaire_started) %>%
      mutate(state = 0)
  )

# Create control cohort for clustering
# With controls now
lc_controls <- covid19_result |>
  recordAttrition() %>%
  dplyr::select("eid","specdate","result") |>
  filter(result == 0) |>
  dplyr::select(-"result") |>
  distinct() |>
    recordAttrition("Restrict to participants with a negative COVID-19 test result") |>
  inner_join(
    health_questionnaire |>
      filter(!is.na(questionnaire_started)),
    by = "eid"
  ) |>
   recordAttrition("Restrict to participants that answered the health and well-being web-questionnaire") |>
 #   filter(questionnaire_started != as.Date("1999-01-01")) |>
  #  recordAttrition("Restrict to participants that answered yes/no to all the questions from the health and well-being web-questionnaire.") |>
  filter(specdate < questionnaire_started) |>
  recordAttrition("Restrict to participants without a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction before answering the questionnaire.")

lc_controls <- lc_controls |>
  group_by(eid) |>
  filter(specdate == max(specdate)) |>
  ungroup() |>
  restoreAttrition(lc_controls) %>%
  recordAttrition("Restrict to the latest record.")

lc_controls <- lc_controls |>
  filter(((specdate+washout_period) < questionnaire_started) & (specdate > (questionnaire_started-365))) |>
   recordAttrition(paste0("Restrict to participants without a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction between 1 year and ",washout_period," days before answering the questionnaire."))

for(s in symptom_cols) {
  var = paste0("symptom_",s)
  varl = paste0("length_",s)
  lc_controls <- lc_controls |>
    dplyr::mutate(
      !!sym(var) := dplyr::if_else(!!sym(var) < 0, NA, !!sym(var)) # next do -1 and -3 differently
    ) |>
    dplyr::mutate(
      !!sym(var) := dplyr::if_else(questionnaire_started - !!sym(varl) < specdate & !is.na(!!sym(varl)) & !is.na(!!sym(var)), 0, !!sym(var)) # if persistent symptom, turn into no symptom
    ) |>
    dplyr::mutate(
      !!sym(varl) := dplyr::if_else(!!sym(var) == 0, NA, !!sym(varl)) # if persistent symptom, length to NA
    )
}

lc_controls <- lc_controls |>
  rowwise() |>
  mutate("symptom" = across(starts_with("symptom")) |> sum(na.rm = TRUE)) |>
  ungroup() %>%
   restoreAttrition(lc_controls) %>%
  mutate(length = do.call(pmin, c(dplyr::select(., starts_with("length")), na.rm = TRUE))) |>
  dplyr::select(- dplyr::starts_with("length_"))

lc_controls <- lc_controls |>
  filter((questionnaire_started-length) > specdate) |>
    recordAttrition("Restrict to people with at least one non persistent symptom") |>
  filter(symptom >= n_symptoms) |>
  recordAttrition("Restrict to people that reported at least one symptom of the questionnaire")

# Remove if they have a positive test between the negative test and the questionnaire
lc_controls_remove <- lc_controls |>
  left_join(
    covid19_result |>
      dplyr::select("eid","specdate","result") |>
      filter(result == 1) |>
      dplyr::select(-"result") |>
      distinct() |>
      rename("coviddate" = "specdate"),
    by = "eid"
  ) |>
  dplyr::select(c("eid", "specdate", "questionnaire_started", "coviddate")) |>
  dplyr::filter(coviddate > specdate & coviddate < questionnaire_started) |>
  dplyr::pull("eid") |>
  unique()

lc_controls <- lc_controls |>
  dplyr::filter(!(eid %in% lc_controls_remove)) %>%
  recordAttrition("Removed people who had a SARS-Cov-2 infection between the negative test and the questionnaire")

# Remove anyone who is considered for the cases directly
lc_controls <- lc_controls |>
  dplyr::anti_join(
    longCovid_cases |>
      dplyr::select("eid"),
    by = "eid"
  ) %>%
  recordAttrition("Removed people who are in the case cohort")

# Save attrition
write.csv(attr(longCovid_cases, "cohort_attrition"), paste0(dir_results,"attrition_longcovid_cases.csv"))
write.csv(attr(longCovid_controls, "cohort_attrition"), paste0(dir_results,"attrition_longcovid_controls.csv"))
write.csv(attr(lc_controls, "cohort_attrition"), paste0(dir_results,"attrition_longcovid_controls_clustering.csv"))

set.seed(11)
# Random sample roughly half the size of both cohorts combined
longCovid_random <- longCovid_cases |>
  dplyr::union_all(lc_controls) |>
  dplyr::slice_sample(n = 21962)

write.csv(lc_controls, paste0(dir_results,"longcovid_controls_clustering.csv"))
write.csv(longCovid_cases, paste0(dir_results,"longcovid_cases_clustering.csv"))
write.csv(longCovid_random, paste0(dir_results,"longcovid_random_clustering.csv"))

# Check again same here as in Marta's GitHub
