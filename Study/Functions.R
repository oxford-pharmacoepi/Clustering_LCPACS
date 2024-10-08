loadHesData <- function(dir_data){
  ukb         <- read.table(paste0(dir_data,"UKBiobank/ukb678283.tab"), header = TRUE, sep = "\t") |> as_tibble() |> select("eid" = "f.eid")
  hesin       <- read.table(paste0(dir_data,"UKBiobank/hesin.txt"), header = TRUE, sep = "\t") |> as_tibble()
  hesin_diag  <- read.table(paste0(dir_data,"UKBiobank/hesin_diag.txt"), header = TRUE, sep = "\t") |> as_tibble()
  
  hes <- ukb |>
    left_join(
      hesin |>
        select("eid", "ins_index", "epistart", "epiend", "admidate") |>
        inner_join(
          hesin_diag |>
            select("eid","ins_index","diag_icd10"),
          by = c("eid", "ins_index"))
    ) |>
    mutate(episode_date = if_else(is.na(epistart), admidate, epistart)) |>
    mutate(episode_date = as.Date(episode_date, format = "%d/%m/%Y")) |>
    select("eid","diag_icd10","episode_date")
  
  return(hes)
}

loadUKBData <- function(dir_data){
  bd <- read.table(paste0(dir_data,"UKBiobank/ukb678283.tab"), header=TRUE, sep="\t")
  bd$f.28001.0.0 <- as.Date(bd$f.28001.0.0)
  bd$f.28001.1.0 <- as.Date(bd$f.28001.1.0)
  bd$f.28001.2.0 <- as.Date(bd$f.28001.2.0)
  bd$f.28001.3.0 <- as.Date(bd$f.28001.3.0)
  bd$f.28001.4.0 <- as.Date(bd$f.28001.4.0)
  bd$f.28001.5.0 <- as.Date(bd$f.28001.5.0)
  bd$f.28005.0.0 <- as.Date(bd$f.28005.0.0)
  bd$f.28005.1.0 <- as.Date(bd$f.28005.1.0)
  bd$f.28005.2.0 <- as.Date(bd$f.28005.2.0)
  bd$f.28005.3.0 <- as.Date(bd$f.28005.3.0)
  bd$f.28005.4.0 <- as.Date(bd$f.28005.4.0)
  bd$f.28005.5.0 <- as.Date(bd$f.28005.5.0)
  bd$f.28006.0.0 <- as.Date(bd$f.28006.0.0)
  bd$f.28006.1.0 <- as.Date(bd$f.28006.1.0)
  bd$f.28006.2.0 <- as.Date(bd$f.28006.2.0)
  bd$f.28006.3.0 <- as.Date(bd$f.28006.3.0)
  bd$f.28006.4.0 <- as.Date(bd$f.28006.4.0)
  bd$f.28006.5.0 <- as.Date(bd$f.28006.5.0)
  bd$f.28008.0.0 <- as.Date(bd$f.28008.0.0)
  bd$f.28008.1.0 <- as.Date(bd$f.28008.1.0)
  bd$f.28008.2.0 <- as.Date(bd$f.28008.2.0)
  bd$f.28008.3.0 <- as.Date(bd$f.28008.3.0)
  bd$f.28008.4.0 <- as.Date(bd$f.28008.4.0)
  bd$f.28008.5.0 <- as.Date(bd$f.28008.5.0)
  bd$f.28009.0.0 <- as.Date(bd$f.28009.0.0)
  bd$f.28009.1.0 <- as.Date(bd$f.28009.1.0)
  bd$f.28009.2.0 <- as.Date(bd$f.28009.2.0)
  bd$f.28009.3.0 <- as.Date(bd$f.28009.3.0)
  bd$f.28009.4.0 <- as.Date(bd$f.28009.4.0)
  bd$f.28009.5.0 <- as.Date(bd$f.28009.5.0)
  bd$f.28030.0.0 <- as.Date(bd$f.28030.0.0)
  bd$f.28030.1.0 <- as.Date(bd$f.28030.1.0)
  bd$f.28030.2.0 <- as.Date(bd$f.28030.2.0)
  bd$f.28030.3.0 <- as.Date(bd$f.28030.3.0)
  bd$f.28030.4.0 <- as.Date(bd$f.28030.4.0)
  bd$f.28030.5.0 <- as.Date(bd$f.28030.5.0)
  bd$f.28030.6.0 <- as.Date(bd$f.28030.6.0)
  bd$f.28031.0.0 <- as.Date(bd$f.28031.0.0)
  bd$f.28031.1.0 <- as.Date(bd$f.28031.1.0)
  bd$f.28031.2.0 <- as.Date(bd$f.28031.2.0)
  bd$f.28031.3.0 <- as.Date(bd$f.28031.3.0)
  bd$f.28031.4.0 <- as.Date(bd$f.28031.4.0)
  bd$f.28031.5.0 <- as.Date(bd$f.28031.5.0)
  bd$f.28032.0.0 <- as.Date(bd$f.28032.0.0)
  bd$f.28032.1.0 <- as.Date(bd$f.28032.1.0)
  bd$f.28032.2.0 <- as.Date(bd$f.28032.2.0)
  bd$f.28032.3.0 <- as.Date(bd$f.28032.3.0)
  bd$f.28032.4.0 <- as.Date(bd$f.28032.4.0)
  bd$f.28032.5.0 <- as.Date(bd$f.28032.5.0)
  bd$f.28032.6.0 <- as.Date(bd$f.28032.6.0)
  bd$f.28032.7.0 <- as.Date(bd$f.28032.7.0)
  bd$f.28032.8.0 <- as.Date(bd$f.28032.8.0)
  bd$f.28033.2.0 <- as.Date(bd$f.28033.2.0)
  bd$f.28033.3.0 <- as.Date(bd$f.28033.3.0)
  bd$f.28033.4.0 <- as.Date(bd$f.28033.4.0)
  bd$f.28033.5.0 <- as.Date(bd$f.28033.5.0)
  bd$f.28033.6.0 <- as.Date(bd$f.28033.6.0)
  bd$f.28141.0.0 <- as.Date(bd$f.28141.0.0)
  bd$f.28142.0.0 <- as.Date(bd$f.28142.0.0)
  bd$f.28143.0.0 <- as.Date(bd$f.28143.0.0)
  bd$f.28146.0.0 <- as.Date(bd$f.28146.0.0)
  bd$f.28166.0.0 <- as.Date(bd$f.28166.0.0)
  bd$f.28167.0.0 <- as.Date(bd$f.28167.0.0)
  bd$f.28754.0.0 <- as.Date(bd$f.28754.0.0)
  bd$f.28756.0.0 <- as.Date(bd$f.28756.0.0)
  
  return(bd)
}

loadHealthQuestionnaire <- function(ukb){
  bd <- ukb |>
    select(
      "eid" = "f.eid",
      "questionnaire_started"            = "f.28754.0.0",
      "symptom_gastrointestinal_issues"          = "f.28606.0.0",
      "length_gastrointestinal_issues"           = "f.28607.0.0",
      "symptom_vision_problems"                  = "f.28609.0.0",
      "length_vision_problems"                   = "f.28610.0.0",
      "symptom_loss_or_change_in_sense_of_smell" = "f.28612.0.0",
      "length_loss_or_change_in_sense_of_smell"  = "f.28613.0.0",
      "symptom_loss_or_change_in_sense_of_taste" = "f.28615.0.0",
      "length_loss_or_change_in_sense_of_taste"  = "f.28616.0.0",
      "symptom_tinnitus"                         = "f.28624.0.0",
      "length_tinnitus"                          = "f.28625.0.0",
      "symptom_hearing_loss"                     = "f.28627.0.0",
      "length_hearing_loss"                      = "f.28628.0.0",
      "symptom_hearing_issues"                   = "f.28630.0.0",
      "length_hearing_issues"                    = "f.28631.0.0",
      "symptom_headaches"                        = "f.28633.0.0",
      "length_headaches"                         = "f.28634.0.0",
      "symptom_chest_pain"                       = "f.28642.0.0",
      "length_chest_pain"                        = "f.28643.0.0",
      "symptom_pain_on_breathing"                = "f.28645.0.0",
      "length_pain_on_breathing"                 = "f.28646.0.0",
      "symptom_abdominal_pain_tummy_ache"        = "f.28648.0.0",
      "length_abdominal_pain_tummy_ache"         = "f.28649.0.0",
      "symptom_muscle_pain_achy_muscles"         = "f.28654.0.0",
      "length_muscle_pain_achy_muscles"          = "f.28655.0.0",
      "symptom_joint_pain_or_swelling_of_joint"  = "f.28657.0.0",
      "length_joint_pain_or_swelling_of_joint"   = "f.28658.0.0",
      "symptom_persistent_cough"                 = "f.28663.0.0",
      "length_persistent_cough"                  = "f.28664.0.0",
      "symptom_tightness_in_the_chest"           = "f.28669.0.0",
      "length_tightness_in_the_chest"            = "f.28670.0.0",
      "symptom_chest_pressure"                   = "f.28672.0.0",
      "length_chest_pressure"                    = "f.28673.0.0",
      "symptom_postural_tachycardia"             = "f.28678.0.0",
      "length_postural_tachycardia"              = "f.28679.0.0",
      "symptom_dizziness_light_headedness"       = "f.28681.0.0",
      "length_dizziness_light_headedness"        = "f.28682.0.0",
      "symptom_shortness_of_breath_or_trouble_breathing" = "f.28684.0.0",
      "length_shortness_of_breath_or_trouble_breathing"  = "f.28685.0.0",
      "symptom_difficulty_sleeping"              = "f.28687.0.0",
      "length_difficulty_sleeping"               = "f.28688.0.0",
      "symptom_unrestful_sleep"                  = "f.28693.0.0",
      "length_unrestful_sleep"                   = "f.28694.0.0",
      "symptom_mild_fatigue"                     = "f.28696.0.0",
      "length_mild_fatigue"                      = "f.28697.0.0",
      "symptom_severe_fatigue"                   = "f.28699.0.0",
      "length_severe_fatigue"                    = "f.28700.0.0",
      "symptom_post_exertional_symptom_exacerbation" = "f.28702.0.0",
      "length_post_exertional_symptom_exacerbation"  = "f.28703.0.0",
      "symptom_new_allergy_or_intolerance"       = "f.28711.0.0",
      "length_new_allergy_or_intolerance"        = "f.28712.0.0",
      "symptom_fever"                            = "f.28714.0.0",
      "length_fever"                             = "f.28715.0.0",
      "symptom_problems_thinking"                = "f.28720.0.0",
      "length_problems_thinking"                 = "f.28721.0.0",
      "symptom_problems_communicating"           = "f.28723.0.0",
      "length_problems_communicating"            = "f.28724.0.0",
      "symptom_problems_relating_to_mood_anxiety_and_emotions" = "f.28726.0.0",
      "length_problems_relating_to_mood_anxiety_and_emotions"  = "f.28727.0.0",
      "symptom_numbness_or_tingling_somewhere_in_the_body"     = "f.28732.0.0",
      "length_numbness_or_tingling_somewhere_in_the_body"      = "f.28733.0.0"
    ) |>
    as_tibble()
  
  bd <- bd %>%
    mutate(questionnaire_started = if_else(rowSums(select(., starts_with("symptom_")) == -1) > 0, as.Date("1999-01-01"), questionnaire_started)) %>%
    mutate(questionnaire_started = if_else(rowSums(select(., starts_with("symptom_")) == -3) > 0, as.Date("1999-01-01"), questionnaire_started)) |>
    mutate(who_abdominal_pain    = symptom_abdominal_pain_tummy_ache,
           length_abdominal_pain = length_abdominal_pain_tummy_ache,
           who_menstrual_and_period_problems    = 0,
           length_menstrual_and_period_problems = 0,
           who_altered_smell_taste    = if_else(symptom_loss_or_change_in_sense_of_smell == 1 | symptom_loss_or_change_in_sense_of_taste == 1, 1, 0),
           length_altered_smell_taste = pmax(length_loss_or_change_in_sense_of_smell, length_loss_or_change_in_sense_of_taste, na.rm = TRUE),
           who_anxiety    = symptom_problems_relating_to_mood_anxiety_and_emotions,
           length_anxiety = length_problems_relating_to_mood_anxiety_and_emotions,
           who_blurred_vision = symptom_vision_problems,
           length_blurred_vision = length_vision_problems,
           who_chest_pain    = if_else(symptom_chest_pain == 1| symptom_pain_on_breathing == 1, 1, 0),
           length_chest_pain = pmax(length_chest_pain, length_pain_on_breathing, na.rm = TRUE),
           who_cognitive_dysfunction_brain_fog    = if_else(symptom_problems_communicating == 1 | symptom_problems_thinking == 1, 1, 0),
           length_cognitive_dysfunction_brain_fog = pmax(length_problems_communicating, length_problems_thinking, na.rm = TRUE),
           who_cough    = symptom_persistent_cough,
           length_cough = length_persistent_cough,
           who_depression    = 0,
           length_depression = 0,
           who_dizziness    = symptom_dizziness_light_headedness,
           length_dizziness = length_dizziness_light_headedness,
           who_fatigue    = if_else(symptom_mild_fatigue == 1 | symptom_severe_fatigue == 1, 1, 0),
           length_fatigue = pmax(length_mild_fatigue, length_severe_fatigue, na.rm = TRUE),
           who_intermittent_fever    = symptom_fever,
           length_intermittent_fever = length_fever,
           who_gastrointerestinal_issues    = symptom_gastrointestinal_issues,
           length_gastrointerestinal_issues = length_gastrointestinal_issues,
           who_headache    = symptom_headaches,
           length_headache = length_headaches,
           who_memory_issues    = 0,
           length_memory_issues = 0,
           who_joint_pain    = symptom_joint_pain_or_swelling_of_joint,
           length_joint_pain = length_joint_pain_or_swelling_of_joint,
           who_muscle_pain_spasms    = symptom_muscle_pain_achy_muscles,
           length_muscle_pain_spasms = length_muscle_pain_achy_muscles,
           who_neuralgias    = 0,
           length_neuralgias = 0,
           who_new_onset_allergies    = symptom_new_allergy_or_intolerance,
           length_new_onset_allergies = length_new_allergy_or_intolerance,
           who_pins_and_needles_sensations    = symptom_numbness_or_tingling_somewhere_in_the_body,
           length_pins_and_needles_sensations = length_numbness_or_tingling_somewhere_in_the_body,
           who_post_exertional_malaise    = symptom_post_exertional_symptom_exacerbation,
           length_post_exertional_malaise = length_post_exertional_symptom_exacerbation,
           who_shortness_of_breath    = symptom_shortness_of_breath_or_trouble_breathing,
           length_shortness_of_breath = length_shortness_of_breath_or_trouble_breathing,
           who_sleep_disorders    = if_else(symptom_difficulty_sleeping == 1 | symptom_unrestful_sleep == 1, 1, 0),
           length_sleep_disorders = pmax(length_difficulty_sleeping, length_unrestful_sleep, na.rm = TRUE),
           who_tachycardia_palpitations    = symptom_postural_tachycardia,
           length_tachycardia_palpitations = length_postural_tachycardia,
           who_tinnitus_and_other_hearing_issues    = if_else(symptom_tinnitus == 1 | symptom_hearing_loss == 1 | symptom_hearing_issues == 1, 1, 0),
           length_tinnitus_and_other_hearing_issues = pmax(length_tinnitus, length_hearing_loss, length_hearing_issues, na.rm = TRUE)
    ) |>
    select(-starts_with("symptom_")) |>
    rename_at(vars(starts_with("who_")), ~gsub("who","symptom",.))
  
  bd <- bd |>
    mutate_at(vars(starts_with("length_")), ~if_else(.  %in% c(-1,-3, 1), 0, .)) |>
    mutate_at(vars(starts_with("length_")), ~if_else(. == 2, 28, .)) |>
    mutate_at(vars(starts_with("length_")), ~if_else(. == 3, 365,.))
  return(bd)
  
}

loadCovid19Result <- function(ukb){
  covid19_result_england  <- read.table(paste0(dir_data,"UKBiobank/covid19_result_england.txt"), header = TRUE) |> as_tibble()
  covid19_result_scoltand <- read.table(paste0(dir_data,"UKBiobank/covid19_result_scotland.txt"), header = TRUE) |> as_tibble()
  covid19_result_wales    <- read.table(paste0(dir_data,"UKBiobank/covid19_result_wales.txt"), header = TRUE) |> as_tibble()
  
  covid19_result <- covid19_result_england |> mutate(laboratory = as.character(laboratory)) |>
    full_join(covid19_result_scoltand) |>
    full_join(covid19_result_wales |> mutate(laboratory = as.character(laboratory))) |>
    mutate(specdate = as.Date(specdate, format = "%d/%m/%Y"))
  
  return(covid19_result)
}

recordAttrition <- function(table, reason = NULL){
  if(is.null(reason)){
    attr(table, "cohort_attrition") <- tibble(
      "number_subjects" = table |> select("eid") |> distinct() |> nrow(),
      "number_records"      = table |> nrow(),
      "reason_id" = 1,
      "reason"       = "Initial qualifying events",
      "excluded_records" = 0,
      "excluded_subjects"      = 0
    )
  }else{
    n_records <- attr(table, "cohort_attrition") |>
      filter(row_number() == max(row_number())) |>
      pull(number_records)
    n_participants <- attr(table, "cohort_attrition") |>
      filter(row_number() == max(row_number())) |>
      pull(number_subjects)
    
    attr(table, "cohort_attrition") <- attr(table, "cohort_attrition") |>
      union_all(
        tibble(
          "number_subjects" = table |> select("eid") |> distinct() |> nrow(),
          "number_records"      = table |> nrow(),
          "reason_id" = max(attr(table, "cohort_attrition")[["reason_id"]])+1,
          "reason"       = reason,
          "excluded_subjects" = n_participants - (table |> select("eid") |> distinct() |> nrow()),
          "excluded_records" = n_records - (table |> nrow())
        )
      ) |>
      distinct()
  }
  return(table)
}

restoreAttrition <- function(table, table_old){
  attr(table, "cohort_attrition") <- attr(table_old, "cohort_attrition")
  return(table)
}

loadSequelaTable <- function(){
  sequela_table <- tibble(
    "organ_system" = as.character(),
    "sequela"      = as.character(),
    "icd10_code"   = as.character()
  ) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "acute_coronary_disease",  "icd10_code" = c("I24","I240","I241","I248","I249")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "angina",                  "icd10_code" = c("I20","I200","I201","I208","I209")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "atrial_fibrillation",     "icd10_code" = c("I480", "I481", "I482")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "atrial_flutter",          "icd10_code" = c("I483", "I484")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "bradycardia",             "icd10_code" = "R001") |>
    add_row("organ_system" = "cardiovascular", "sequela" = "cardiac_arrest",          "icd10_code" = c("I46","I460","I461","I469")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "cariogenic_shock",        "icd10_code" = "R570") |>
    add_row("organ_system" = "cardiovascular", "sequela" = "heart_failure",           "icd10_code" = c("I50","I500","I501","I509")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "ischemic_cardiomyopathy", "icd10_code" = "I255") |>
    add_row("organ_system" = "cardiovascular", "sequela" = "myocardial_infarction",   "icd10_code" = c("I21","I210","I211","I212","I213","I214","I219","I21X",
                                                                                                       "I22","I220","I221","I228","I229")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "myocarditis",             "icd10_code" = "I514") |>
    add_row("organ_system" = "cardiovascular", "sequela" = "non_ischemic_cardiomyopathy", "icd10_code" = c("I42","I420","I421","I422","I423","I424","I425","I426","I427","I428","I429",
                                                                                                           "I43","I430","I431","I432","I438",
                                                                                                           "B332")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "pericarditis",            "icd10_code" = c("I30","I300","I301","I308","I309",
                                                                                                       "I311","I312", "I313", "I318", "I319",
                                                                                                       "I32","I320","I321","I328")) |>
    add_row("organ_system" = "cardiovascular", "sequela" = "tachycardia",             "icd10_code" = "R000") |>
    add_row("organ_system" = "cardiovascular", "sequela" = "ventricular_arrhythmias", "icd10_code" = c("I490", "I470", "I471", "I472")) |>
    add_row("organ_system" = "coagulation", "sequela" = "anemia",               "icd10_code" = c("D60","D600","D601","D608","D609",
                                                                                                 "D61","D610","D611","D612","D613","D618","D619",
                                                                                                 "D63","D630","D631","D638",
                                                                                                 "D64","D640","D641","D642","D643","D644","D648","D649",
                                                                                                 "D62", "D62X")) |>
    add_row("organ_system" = "coagulation", "sequela" = "coagulation_defect",   "icd10_code" = c("D689")) |>
    add_row("organ_system" = "coagulation", "sequela" = "deep_vein_thrombosis", "icd10_code" = c("I801", "I802", "I803", "I81","I81X")) |>
    add_row("organ_system" = "coagulation", "sequela" = "pulmonary_embolism",   "icd10_code" = c("I26","I260","I269")) |>
    add_row("organ_system" = "coagulation", "sequela" = "venous_thrombotic_embolism",  "icd10_code" = c("I82","I820", "I822", "I823", "I828","I829"))
  return(sequela_table)
}

loadGeneticData <- function(ukb){
  ukb |>
    select("eid" = "f.eid", 
           "sex" = "f.31.0.0",
           "genetic_ethnic_grouping" = "f.22006.0.0",
           "genetic_sex" = "f.22001.0.0",
           "pc1" = "f.22009.0.1", "pc2" = "f.22009.0.2", "pc3" = "f.22009.0.3",
           "pc4" = "f.22009.0.4", "pc5" = "f.22009.0.5", "pc6" = "f.22009.0.6",
           "pc7" = "f.22009.0.7", "pc8" = "f.22009.0.8", "pc9" = "f.22009.0.9",
           "pc10" = "f.22009.0.10",
           "batch" = "f.22000.0.0", "f.22019.0.0",
           "year_birth" = "f.34.0.0",
           "sex_chromosome_aneuploidy" = "f.22019.0.0") |>
    as_tibble()
}
