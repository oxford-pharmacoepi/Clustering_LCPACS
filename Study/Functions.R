tibbleNameBiomarkers <- function(){
  tibble::tibble(
    "risk_factor"  = c("Cardiovascular",
                       "apolipoprotein_a", "apolipoprotein_b", "c_reactive_protein",
                       "cholesterol", "ldl_direct", "hdl_cholesterol", "lipoprotein_a",
                       "triglycerides",

                       "Bone and joint",
                       "alkaline_phosphatase","calcium", "rheumatoid_factor", "vitamin_d",

                       "Cancer",
                       "igf", "oestradiol", "shbg", "testosterone",

                       "Diabetes",
                       "glucose", "hbauc",

                       "Renal",
                       "creatinine", "cystatin_c", "phosphate", "total_protein","urate", "urea",

                       "Liver",
                       "alanine_aminotransferase", "albumin", "aspartate_aminotransferase",
                       "direct_bilirubin","gamma_glutamyltransferase", "total_bilirubin"
                       ),

    "risk_factor1" = c("Cardiovascular",
                       "Apolipoprotein A", "Apolipoprotein B", "C-reactive protein",
                       "Cholesterol", "LDL direct", "HDL-Cholesterol", "Lipoprotein (a)",
                       "Triglyceride",

                       "Bone and joint",
                       "Alkaline phosphatase", "Calcium", "Rheumatoid factor", "Vitamin D",

                       "Cancer",
                       "IGF-1", "Oestradiol", "SHBG", "Testosterone",

                       "Diabetes",
                       "Glucose", "HbA1c",

                       "Renal",
                       "Creatinine", "Cystatin C",  "Phosphate", "Total protein", "Urate", "Urea",

                       "Liver",
                       "Alanine aminotransferase", "Albumin", "Aspartate aminotransferase",
                       "Direct bilirubin", "Gamma glutamyltransferase", "Total bilirubin"
                       )
  )
}

tibbleNameComorbidities <- function(){
  tibble(
    "risk_factor"  = c("aids", "asthma", "cancer", "cancer_metastatic", "cerebrovascular_disease",
                       "congestive_heart_failure", "chronic_kidney_disease", "copd", "dementia", "diabetes",
                       "diabetes_organ_damage", "fracture", "hemiplegia", "liver_disease_mild",
                       "liver_disease_moderate_to_severe", "myocardial_infarction",
                       "peptic_ulcer","peripheral_vascular_disease", "rheumatoid_arthritis"),
    "risk_factor1" = c("AIDS", "Asthma", "Cancer", "Cancer - metastatic", "Cerebrovascular disease",
                       "Congestive heart failure", "Chronic kidney disease", "COPD", "Dementia",
                       "Diabetes", "Diabetes - organ damage", "Fracture", "Hemiplegia", "Liver disease - mild",
                       "Liver disease - moderate to severe",  "Myocardial infarction",
                       "Peptic ulcer", "Peripheral vascular disease","Rheumatoid arthritis"
    ))
}

loadHesData <- function(){
  ukb         <- read.table(paste0("UKBiobank/healthAndWellBeingQuestionnaire.tab"), header = TRUE, sep = "\t") |> as_tibble() |> dplyr::select("eid" = "f.eid")
  hesin       <- read.table(paste0("UKBiobank/hesin.txt"), header = TRUE, sep = "\t") |> as_tibble()
  hesin_diag  <- read.table(paste0("UKBiobank/hesin_diag.txt"), header = TRUE, sep = "\t") |> as_tibble()

  hes <- ukb |>
    left_join(
      hesin |>
        dplyr::select("eid", "ins_index", "epistart", "epiend", "admidate") |>
        inner_join(
          hesin_diag |>
            dplyr::select("eid","ins_index","diag_icd10"),
          by = c("eid", "ins_index"))
    ) |>
    mutate(episode_date = if_else(is.na(epistart), admidate, epistart)) |>
    mutate(episode_date = as.Date(episode_date, format = "%d/%m/%Y")) |>
    dplyr::select("eid","diag_icd10","episode_date")

  return(hes)
}

addCoding <- function(bd, variable){
  coding <- tibble(read.table(paste0("UKBiobank/",variable,"_coding.tsv"), header = TRUE, sep = "\t")) |>
    dplyr::select(!!variable := "coding", "meaning")

  bd <- bd |>
    left_join(
      coding,
      by = variable
    ) |>
    dplyr::select(-!!variable) |>
    rename(!!variable := "meaning") |>
    mutate(!!variable := if_else(.data[[variable]] %in% c("Prefer not to answer", "Do not know"), NA, .data[[variable]]))
  return(bd)
}

cont <- function(variable){
  paste0(prettyNum(round(mean(variable, na.rm = TRUE),2), big.mark = ","), " (", round(sd(variable, na.rm = TRUE),2), ")")
}

miss <- function(variable){
  if(sum(is.na(variable), na.rm = TRUE) == 0){
    x <- "0 (0)"
  }else{
    x <- paste0(prettyNum(sum(is.na(variable), na.rm = TRUE), big.mark = ","), " (", round(sum(is.na(variable), na.rm = TRUE)/length(variable)*100,2),")")
  }
  return(x)
}

bin  <- function(variable, level){
  paste0(prettyNum(sum(variable == level, na.rm = TRUE),big.mark = ","), " (", round(sum(variable == level, na.rm = TRUE)/length(variable)*100,2),")")
}

counts <- function(variable){
  paste0(formatC(sum(!is.na(variable), na.rm = TRUE), big.mark = ","), " (", round(sum(!is.na(variable), na.rm = TRUE)/length(variable)*100,2),")")
}

quantiles <- function(variable){
  return(
    paste0(
      formatC(format = "d", quantile(variable, probs = 0.05, na.rm = TRUE), big.mark = ",", digits = 2),", ",
      formatC(format = "d", quantile(variable, probs = 0.25, na.rm = TRUE), big.mark = ",", digits = 2),", ",
      formatC(format = "d", quantile(variable, probs = 0.5, na.rm = TRUE),  big.mark = ",", digits = 2),", ",
      formatC(format = "d", quantile(variable, probs = 0.75, na.rm = TRUE), big.mark = ",", digits = 2),", ",
      formatC(format = "d", quantile(variable, probs = 0.95, na.rm = TRUE), big.mark = ",", digits = 2)
    )

  )
}

quantilesBiomarkers <- function(variable){
  return(
    paste0(
      formatC(format = "f", quantile(variable, probs = 0.05, na.rm = TRUE), big.mark = ",", digits = 2),", ",
      formatC(format = "f", quantile(variable, probs = 0.25, na.rm = TRUE), big.mark = ",", digits = 2),", ",
      formatC(format = "f", quantile(variable, probs = 0.5, na.rm = TRUE),  big.mark = ",", digits = 2),", ",
      formatC(format = "f", quantile(variable, probs = 0.75, na.rm = TRUE), big.mark = ",", digits = 2),", ",
      formatC(format = "f", quantile(variable, probs = 0.95, na.rm = TRUE), big.mark = ",", digits = 2)
    )

  )
}

nam <- function(variable){
  return(" ")
}

summariseBiomarkers <- function(cohort){
  cohort |>
    summarise(across(
      .cols = -c("eid","X"),
      .fns = list(
        "Name" = nam,
        "Mean_(SD)" = cont,
        "q05,_q25,_q50,_q75,_q95" = quantilesBiomarkers,
        "Counts" = counts,
        "Missings" = miss),
      .names = "{.col} {.fn}"
    )) |>
    mutate_all(~as.character(.)) |>
    pivot_longer(everything(), names_to = "Risk factor", values_to = "UK Biobank dataset") |>
    mutate(Name = gsub(" .*","",`Risk factor`)) |>
    mutate(`Risk factor` = if_else(!grepl("Name",`Risk factor`), gsub(".* ","",`Risk factor`), `Risk factor`)) |>
    mutate(`Risk factor` = gsub("biomarker_","",`Risk factor`)) |>
    mutate(`Risk factor` = gsub("_", " ", `Risk factor`)) |>
    mutate(`Risk factor` = gsub("Sd","SD",`Risk factor`)) |>
    mutate(`Risk factor` = gsub("name","",`Risk factor`)) |>
    rowwise() |>
    mutate(z = list(cohort[[Name]][!is.na(cohort[[Name]])])) |>
    group_by(Name) |>
    mutate(z = if_else(row_number() == 1, z, list(NULL))) |>
    ungroup()
}

loadComorbidities <- function(comorbiditiesCrude){
  ukb <- read.table(paste0("UKBiobank/healthAndWellBeingQuestionnaire.tab"),
                    header = TRUE, sep = "\t") |> as_tibble() |> dplyr::select("eid" = "f.eid")

  cohort <- ukb |>
    left_join(
      as_tibble(read.csv(paste0(dir_results,"/LoadCovariates_commorbidities.csv"))) |>
        dplyr::select(-c("X","episode_date")) |>
        mutate(state = 1) |>
        pivot_wider(names_from = value, values_from = state) |>
        mutate(across(-c("eid"), ~if_else(is.na(.), 0, .))) |>
        distinct(),
      by = "eid") |>
    mutate(across(-c("eid"), ~if_else(is.na(.), 0,.))) |>
    mutate(fracture = if_else(
      ankle == 1 | ankle_open == 1 | elbow == 1 | elbow_open == 1 | femur_distal == 1 |
        femur_subtroch_shaft == 1 | femur_subtroch_shaft_history == 1 | knee == 1 |
        nhip_other == 1 | nhip_other_history == 1 | nhip_other_open == 1 | pelvis == 1 |
        pelvis_open == 1 | pelvis_spinalcord == 1 | radius_ulna_open == 1 | rib == 1 |
        shoulder == 1 | shoulder_open == 1 | spine == 1 | spine_history == 1 | spine_open == 1 |
        tibia == 1 | tibia_open == 1 | tibia_prox == 1 | wrist_forearm == 1 | wrist_forearm_open == 1 |
        hip == 1 | hip_open == 1 | foot == 1 | foot_open == 1 | lumbar_spine_pelvis == 1, 1, 0
    )) |>
    dplyr::select(-c("ankle", "ankle_open", "elbow", "elbow_open", "femur_distal", "femur_subtroch_shaft", "femur_subtroch_shaft_history",
              "knee", "nhip_other", "nhip_other_history", "nhip_other_open", "pelvis", "pelvis_open",
              "pelvis_spinalcord", "radius_ulna_open", "rib", "shoulder", "shoulder_open", "spine", "spine_history",
              "spine_open", "tibia", "tibia_open", "tibia_prox", "wrist_forearm", "wrist_forearm_open",
              "hip", "hip_open", "foot", "foot_open", "lumbar_spine_pelvis"))

  return(cohort)
}

summariseComorbidities <- function(commorbidities){
  commorbidities |>
    mutate(across(-c("eid"), ~if_else(. == 1, "Cases","Controls"))) |>
    summarise(
      "Acquired immunodeficiency syndrome (AIDS) (%)_aids" = " ",
      "\tCases_aids" = bin(aids, "Cases"),
      "\tControls_aids" = bin(aids, "Controls"),
      "Counts_aids"      = counts(aids),
      "Missings_aids" = miss(aids),

      "Asthma (%)_asthma" = " ",
      "\tCases_asthma" = bin(asthma, "Cases"),
      "\tControls_asthma" = bin(asthma, "Controls"),
      "Counts_asthma"      = counts(asthma),
      "Missings_asthma" = miss(asthma),

      "Cancer (%)_cancer" = " ",
      "\tCases_cancer" = bin(cancer, "Cases"),
      "\tControls_cancer" = bin(cancer, "Controls"),
      "Counts_cancer"      = counts(cancer),
      "Missings_cancer" = miss(cancer),

      "Cancer - metastatic (%)_cancer metastatic" = " ",
      "\tCases_cancer metastatic" = bin(cancer_metastatic, "Cases"),
      "\tControls_cancer metastatic" = bin(cancer_metastatic, "Controls"),
      "Counts_cancer metastatic"     = counts(cancer_metastatic),
      "Missings_cancer metastatic"   = miss(cancer_metastatic),

      "Cerebrovascular disease (%)_cerebrovascular disease" = " ",
      "\tCases_cerebrovascular disease" = bin(cerebrovascular_disease, "Cases"),
      "\tControls_cerebrovascular disease" = bin(cerebrovascular_disease, "Controls"),
      "Counts_cerebrovascular disease"     = counts(cerebrovascular_disease),
      "Missings_cerebrovascular disease"   = miss(cerebrovascular_disease),

      "Congestive heart failure (%)_congestive heart failure" = " ",
      "\tCases_congestive heart failure" = bin(congestive_heart_failure, "Cases"),
      "\tControls_congestive heart failure" = bin(congestive_heart_failure, "Controls"),
      "Counts_congestive heart failure"     = counts(congestive_heart_failure),
      "Missings_congestive heart failure"   = miss(congestive_heart_failure),

      "Chronic kidney disease (%)_chronic kidney disease" = " ",
      "\tCases_chronic kidney disease" = bin(chronic_kidney_disease, "Cases"),
      "\tControls_chronic kidney disease" = bin(chronic_kidney_disease, "Controls"),
      "Counts_chronic kidney disease"     = counts(chronic_kidney_disease),
      "Missings_chronic kidney disease"   = miss(chronic_kidney_disease),

      "Chronic obstructive pulmonary disease (COPD) (%)_copd" = " ",
      "\tCases_copd" = bin(copd, "Cases"),
      "\tControls_copd" = bin(copd, "Controls"),
      "Counts_copd"     = counts(copd),
      "Missings_copd"   = miss(copd),

      "Dementia (%)_dementia" = " ",
      "\tCases_dementia" = bin(dementia, "Cases"),
      "\tControls_dementia" = bin(dementia, "Controls"),
      "Counts_dementia"      = counts(dementia),
      "Missings_dementia" = miss(dementia),

      "Diabetes (%)_diabetes" = " ",
      "\tCases_diabetes" = bin(diabetes, "Cases"),
      "\tControls_diabetes" = bin(diabetes, "Controls"),
      "Counts_diabetes"      = counts(diabetes),
      "Missings_diabetes" = miss(diabetes),

      "Diabetes - organ damage (%)_diabetes organ damage" = " ",
      "\tCases_diabetes organ damage" = bin(diabetes_organ_damage, "Cases"),
      "\tControls_diabetes organ damage" = bin(diabetes_organ_damage, "Controls"),
      "Counts_diabetes organ damage"     = counts(diabetes_organ_damage),
      "Missings_diabetes organ damage"   = miss(diabetes_organ_damage),

      "Fracture (%)_fracture" = " ",
      "\tCases_fracture" = bin(fracture, "Cases"),
      "\tControls_fracture" = bin(fracture, "Controls"),
      "Counts_fracture"     = counts(fracture),
      "Missings_fracture"   = miss(fracture),

      "Hemiplegia (%)_hemiplegia" = " ",
      "\tCases_hemiplegia" = bin(hemiplegia, "Cases"),
      "\tControls_hemiplegia" = bin(hemiplegia, "Controls"),
      "Counts_hemiplegia"     = counts(hemiplegia),
      "Missings_hemiplegia"   = miss(hemiplegia),

      "Liver disease - mild (%)_liver disease mild" = " ",
      "\tCases_liver disease mild" = bin(liver_disease_mild, "Cases"),
      "\tControls_liver disease mild" = bin(liver_disease_mild, "Controls"),
      "Counts_liver disease mild"     = counts(liver_disease_mild),
      "Missings_liver disease mild"   = miss(liver_disease_mild),

      "Liver disease - moderate to severe (%)_liver disease moderate to severe" = " ",
      "\tCases_liver disease moderate to severe" = bin(liver_disease_moderate_to_severe, "Cases"),
      "\tControls_liver disease moderate to severe" = bin(liver_disease_moderate_to_severe, "Controls"),
      "Counts_liver disease moderate to severe"     = counts(liver_disease_moderate_to_severe),
      "Missings_liver disease moderate to severe"   = miss(liver_disease_moderate_to_severe),

      "Myocardial infarction (%)_myocardial infarction" = " ",
      "\tCases_myocardial infarction" = bin(myocardial_infarction, "Cases"),
      "\tControls_myocardial infarction" = bin(myocardial_infarction, "Controls"),
      "Counts_myocardial infarction"     = counts(myocardial_infarction),
      "Missings_myocardial infarction"   = miss(myocardial_infarction),

      "Peptic ulcer (%)_peptic ulcer" = " ",
      "\tCases_peptic ulcer" = bin(peptic_ulcer, "Cases"),
      "\tControls_peptic ulcer" = bin(peptic_ulcer, "Controls"),
      "Counts_peptic ulcer"     = counts(peptic_ulcer),
      "Missings_peptic ulcer"   = miss(peptic_ulcer),

      "Peripheral vascular disease (%)_peripheral vascular disease" = " ",
      "\tCases_peripheral vascular disease"     = bin(peripheral_vascular_disease, "Cases"),
      "\tControls_peripheral vascular disease"  = bin(peripheral_vascular_disease, "Controls"),
      "Counts_peripheral vascular disease"      = counts(peripheral_vascular_disease),
      "Missings_peripheral vascular disease"    = miss(peripheral_vascular_disease),

      "Rheumatoid arthritis (%)_rheumatoid arthritis" = " ",
      "\tCases_rheumatoid arthritis" = bin(rheumatoid_arthritis, "Cases"),
      "\tControls_rheumatoid arthritis" = bin(rheumatoid_arthritis, "Controls"),
      "Counts_rheumatoid arthritis"     = counts(rheumatoid_arthritis),
      "Missings_rheumatoid arthritis"   = miss(rheumatoid_arthritis)
    ) |>
    mutate_all(~as.character(.)) |>
    pivot_longer(everything(), names_to = "Risk factor", values_to = "UK Biobank dataset") |>
    mutate(Name = gsub(" ","_",gsub(".*_","",`Risk factor`))) |>
    mutate(`Risk factor` = gsub("_.*","",`Risk factor`)) |>
    rowwise() |>
    mutate(z = list(commorbidities[[Name]][!is.na(commorbidities[[Name]])])) |>
    ungroup() |>
    group_by(Name) |>
    mutate(z = if_else(row_number() == 1, z, list(NULL))) |>
    ungroup()

}

createFlextable <- function(table,name){
  y <- table |>
    flextable(col_keys = c("Risk factor", "UK Biobank dataset", "Distribution", "D")) |>
    bold(bold = TRUE, part = "header") |>
    bold(i = c(1:nrow(table ))[(table |> group_by(Name) |> mutate(row = row_number()) |> pull(row)) == 1], bold = TRUE, part = "body") |>
    hline(i = c(1:nrow(table))[table  |> group_by(Name) |> mutate(row = (row_number() == max(row_number()))) |> pull(row)], part = "body")

  l1 <- (c(1:nrow(table))[(table |> group_by(Name) |> mutate(row = row_number()) |> pull(row)) == 1])

  for(ii in 2:(length(l1)-1)){
    y <- y |>
      merge_at(i = l1[ii]:(l1[ii+1]-1), part = "body", j = 3) |>
      merge_at(i = l1[ii]:(l1[ii+1]-1), part = "body", j = 4) |>
      merge_at(i = l1[ii], part = "body", j = 1:2)
  }

  l1 <- c(1:nrow(table))[table$`Risk factor` %in% c("Baseline characteristics", "Biomarkers", "Comorbidities", "Vaccination status")]

  for(ii in 1:length(l1)){
    y <- y |>
      merge_at(i = l1[ii], part = "body", j = c(1:4)) |>
      bg(i = l1[ii], part = "body", bg = "#F2F2F2", j = c(1:4))
  }

  y <- y |>
    width(j = 1, width = 4.25, unit = "cm") |>
    width(j = 2, width = 3.25, unit = "cm") |>
    width(j = c(3,4), width = 4.5, unit = "cm") |>
    bg(i = 1, part = "header", bg = "#D9D9D9", j = c(1:4))

  save_as_docx(y, path = paste0(dir_results,"Tables/",name,".docx"))

  return(y)
}

getFigures <- function(y, folder){
  for(i in seq_len(nrow(y))){
    if(!is.null(y$z[[i]])){
      temp <- tibble(
        value = y$z[[i]],
        name  = gsub(" \\(.*","",y$`Risk factor`[[i]]),
        type  = y$type[[i]]
      )

      g <- gg_type(temp) +
        labs(x = "", y = "") +
        scale_y_continuous(expand = c(0,0), breaks = c()) +
        scale_x_continuous(breaks = c())

      ggsave(plot = g, filename = paste0(dir_results,"Figures/",folder,"/",temp$name[1],"_",temp$type[1],".png"), dpi = 400, height = 11, units = "cm", width = 20)
    }
  }

}

gg_type <- function(temp){
  if(temp$type[1] == "histogram"){
    order <- tibble(
      value = as.character(c(1,0))) |>
      full_join(
        tibble(value = c("Female", "Male")),
        by = "value"
      ) |>
      full_join(
        tibble(value = c("White", "British", "Irish", "Any other white background_ethnic background",
                         "Mixed", "White and Black Caribbean", "White and Black African",
                         "White and Asian", "Any other mixed background", "Asian or Asian British",
                         "Indian", "Pakistani", "Bangladeshi", "Any other Asian background",
                         "Black or Black British", "Caribbean", "African", "Any other Black background",
                         "Chinese", "Other ethnic group")),
        by = "value"
      ) |>
      full_join(
        tibble(value = c("very_deprived","deprived","average","affluent","very_affluent")),
        by = "value"
      ) |>
      full_join(tibble(value = c("Never", "Previous", "Current"),),
                by = "value") |>
      full_join(tibble(value = c("Cases","Controls")),
                by = "value") |>
      full_join(tibble(value = c("underweight", "average", "overweight", "obesity")),
                by = "value") |>
      full_join(tibble(value = c("White","Non-white")),
                by = "value")

    order <- order |>
      filter(value %in% unique(temp$value)) |>
      mutate(order = as.double(row_number()))

    data <- tibble(value = temp$value) |>
      group_by(value) |>
      tally() |>
      mutate(value = as.character(value)) |>
      inner_join(
        order,
        by = "value"
      )

  }else{
    data <- tibble("value" = temp$value)
  }

  geom <- switch(temp$type[1],
                 "line"       = geom_density(aes(x = value), color = "red", size = 1, show.legend = FALSE),
                 "histogram"  = geom_bar(aes(x = order, y = n, fill = order), stat = "identity", show.legend = FALSE),
                 "boxplot"    = geom_boxplot(aes(x = value), outlier.colour = "red", color = "red", fill = "white")
  )

  ggplot(data) +
    geom +
    theme_void()
}

loadHealthAndWellBeingQuestionnaire <- function(){
  bd <- read.table(paste0("UKBiobank/ukb678283.tab"), header = TRUE, sep = "\t")
  bd$f.28754.0.0 <- as.Date(bd$f.28754.0.0)
  bd$f.28756.0.0 <- as.Date(bd$f.28756.0.0)
 # bd$f.28755.0.0 <- as.Date(bd$f.28755.0.0)

  bd <- bd |>
    dplyr::select(
      "eid" = "f.eid",
      "questionnaire_started"            = "f.28754.0.0",
   #   "questionnaire_completed"            = "f.28755.0.0",
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
      "symptom_post_exertional_exacerbation" = "f.28702.0.0",
      "length_post_exertional_exacerbation"  = "f.28703.0.0",
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
      "length_numbness_or_tingling_somewhere_in_the_body"      = "f.28733.0.0",
      "symptom_decrease_in_appetite" = "f.28600.0.0",
      "length_decrease_in_appetite" = "f.28601.0.0",
      "symptom_nausea_and_or_vomiting" = "f.28603.0.0",
      "length_nausea_and_or_vomiting" = "f.28604.0.0",
      "symptom_sore_painful_throat" = "f.28618.0.0",
      "length_sore_painful_throat" = "f.28619.0.0",
      "symptom_nasal_congestion" = "f.28621.0.0",
      "length_nasal_congestion" = "f.28622.0.0",
      "symptom_neck_pain" = "f.28636.0.0",
      "length_neck_pain" = "f.28637.0.0",
      "symptom_back_pain" = "f.28639.0.0",
      "length_back_pain" = "f.28640.0.0",
      "symptom_leg_pain" = "f.28651.0.0",
      "length_leg_pain" = "f.28652.0.0",
      "symptom_bone_pain" = "f.28660.0.0",
      "length_bone_pain" = "f.28661.0.0",
      "symptom_chesty_cough" = "f.28666.0.0",
      "length_chesty_cough" = "f.28667.0.0",
      "symptom_palpitations" = "f.28675.0.0",
      "length_palpitations" = "f.28676.0.0",
      "symptom_night_sweats" = "f.28687.0.0",
      "length_night_sweats" = "f.28688.0.0",
      "symptom_skin_issues" = "f.28705.0.0",
      "length_skin_issues" = "f.28706.0.0",
      "symptom_sores_blisters_feet" = "f.28708.0.0",
      "length_sores_blisters_feet" = "f.28709.0.0",
      "symptom_chills" = "f.28717.0.0",
      "length_chills" = "f.28718.0.0",
      "symptom_weakness_muscles" = "f.28729.0.0",
      "length_weakness_muscles" = "f.28730.0.0"
    ) |>
    as_tibble()

  #bd <- bd %>%
  #  mutate(questionnaire_started = if_else(rowSums(select(., starts_with("symptom_")) == -1) > 0, as.Date("1999-01-01"), questionnaire_started)) %>%
  #  mutate(questionnaire_started = if_else(rowSums(select(., starts_with("symptom_")) == -3) > 0, as.Date("1999-01-01"), questionnaire_started))# |>
    # mutate(who_abdominal_pain    = symptom_abdominal_pain_tummy_ache,
    #        length_abdominal_pain = length_abdominal_pain_tummy_ache,
    #        who_menstrual_and_period_problems    = 0,
    #        length_menstrual_and_period_problems = 0,
    #        who_altered_smell_taste    = if_else(symptom_loss_or_change_in_sense_of_smell == 1 | symptom_loss_or_change_in_sense_of_taste == 1, 1, 0),
    #        length_altered_smell_taste = pmax(length_loss_or_change_in_sense_of_smell, length_loss_or_change_in_sense_of_taste, na.rm = TRUE),
    #        who_anxiety    = symptom_problems_relating_to_mood_anxiety_and_emotions,
    #        length_anxiety = length_problems_relating_to_mood_anxiety_and_emotions,
    #        who_blurred_vision = symptom_vision_problems,
    #        length_blurred_vision = length_vision_problems,
    #        who_chest_pain    = if_else(symptom_chest_pain == 1| symptom_pain_on_breathing == 1, 1, 0),
    #        length_chest_pain = pmax(length_chest_pain, length_pain_on_breathing, na.rm = TRUE),
    #        who_cognitive_dysfunction_brain_fog    = if_else(symptom_problems_communicating == 1 | symptom_problems_thinking == 1, 1, 0),
    #        length_cognitive_dysfunction_brain_fog = pmax(length_problems_communicating, length_problems_thinking, na.rm = TRUE),
    #        who_cough    = symptom_persistent_cough,
    #        length_cough = length_persistent_cough,
    #        who_depression    = 0,
    #        length_depression = 0,
    #        who_dizziness    = symptom_dizziness_light_headedness,
    #        length_dizziness = length_dizziness_light_headedness,
    #        who_fatigue    = if_else(symptom_mild_fatigue == 1 | symptom_severe_fatigue == 1, 1, 0),
    #        length_fatigue = pmax(length_mild_fatigue, length_severe_fatigue, na.rm = TRUE),
    #        who_intermittent_fever    = symptom_fever,
    #        length_intermittent_fever = length_fever,
    #        who_gastrointerestinal_issues    = symptom_gastrointestinal_issues,
    #        length_gastrointerestinal_issues = length_gastrointestinal_issues,
    #        who_headache    = symptom_headaches,
    #        length_headache = length_headaches,
    #        who_memory_issues    = 0,
    #        length_memory_issues = 0,
    #        who_joint_pain    = symptom_joint_pain_or_swelling_of_joint,
    #        length_joint_pain = length_joint_pain_or_swelling_of_joint,
    #        who_muscle_pain_spasms    = symptom_muscle_pain_achy_muscles,
    #        length_muscle_pain_spasms = length_muscle_pain_achy_muscles,
    #        who_neuralgias    = 0,
    #        length_neuralgias = 0,
    #        who_new_onset_allergies    = symptom_new_allergy_or_intolerance,
    #        length_new_onset_allergies = length_new_allergy_or_intolerance,
    #        who_pins_and_needles_sensations    = symptom_numbness_or_tingling_somewhere_in_the_body,
    #        length_pins_and_needles_sensations = length_numbness_or_tingling_somewhere_in_the_body,
    #        who_post_exertional_malaise    = symptom_post_exertional_symptom_exacerbation,
    #        length_post_exertional_malaise = length_post_exertional_symptom_exacerbation,
    #        who_shortness_of_breath    = symptom_shortness_of_breath_or_trouble_breathing,
    #        length_shortness_of_breath = length_shortness_of_breath_or_trouble_breathing,
    #        who_sleep_disorders    = if_else(symptom_difficulty_sleeping == 1 | symptom_unrestful_sleep == 1, 1, 0),
    #        length_sleep_disorders = pmax(length_difficulty_sleeping, length_unrestful_sleep, na.rm = TRUE),
    #        who_tachycardia_palpitations    = symptom_postural_tachycardia,
    #        length_tachycardia_palpitations = length_postural_tachycardia,
    #        who_tinnitus_and_other_hearing_issues    = if_else(symptom_tinnitus == 1 | symptom_hearing_loss == 1 | symptom_hearing_issues == 1, 1, 0),
    #        length_tinnitus_and_other_hearing_issues = pmax(length_tinnitus, length_hearing_loss, length_hearing_issues, na.rm = TRUE)
    # ) |>
  #  select(-starts_with("symptom_")) |>
 #   rename_at(vars(starts_with("who_")), ~gsub("who","symptom",.))

  bd <- bd |>
    mutate_at(vars(starts_with("length_")), ~if_else(.  %in% c(-1,-3, 1), 0, .)) |>
    mutate_at(vars(starts_with("length_")), ~if_else(. == 2, 28, .)) |>
    mutate_at(vars(starts_with("length_")), ~if_else(. == 3, 365,.))
  return(bd)
}

loadCovid19Result <- function(){
  covid19_result_england  <- read.table(paste0("UKBiobank/covid19_result_england.txt"), header = TRUE) |> as_tibble()
  covid19_result_scoltand <- read.table(paste0("UKBiobank/covid19_result_scotland.txt"), header = TRUE) |> as_tibble()
  covid19_result_wales    <- read.table(paste0("UKBiobank/covid19_result_wales.txt"), header = TRUE) |> as_tibble()

  covid19_result <- covid19_result_england |> mutate(laboratory = as.character(laboratory)) |>
    full_join(covid19_result_scoltand) |>
    full_join(covid19_result_wales |> mutate(laboratory = as.character(laboratory))) |>
    mutate(specdate = as.Date(specdate, format = "%d/%m/%Y"))

  return(covid19_result)
}

recordAttrition <- function(table, reason = NULL){
  if(is.null(reason)){
    attr(table, "cohort_attrition") <- tibble(
      "number_subjects" = table |> dplyr::select("eid") |> distinct() |> nrow(),
      "number_records"      = table |> nrow(),
      "reason_id" = 1,
      "reason"       = "Initial qualifying events",
      "excluded_records" = 0,
      "excluded_subjects"      = 0
    )
  }else{
    n_records <- attr(table, "cohort_attrition") |>
      dplyr::filter(row_number() == max(row_number())) |>
      dplyr::pull(number_records)
    n_participants <- attr(table, "cohort_attrition") |>
      dplyr::filter(row_number() == max(row_number())) |>
      dplyr::pull(number_subjects)

    attr(table, "cohort_attrition") <- attr(table, "cohort_attrition") |>
      dplyr::union_all(
        tibble(
          "number_subjects" = table |> dplyr::select("eid") |> distinct() |> nrow(),
          "number_records"      = table |> nrow(),
          "reason_id" = max(attr(table, "cohort_attrition")[["reason_id"]])+1,
          "reason"       = reason,
          "excluded_subjects" = n_participants - (table |> dplyr::select("eid") |> distinct() |> nrow()),
          "excluded_records" = n_records - (table |> nrow())
        )
      ) |>
      dplyr::distinct()
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

tableOneDemographics <- function(cohort){
  table <- cohort |>
    dplyr::select(-starts_with("biomarker_")) |>
    summarise("N_total"  = prettyNum(paste0(nrow(cohort)), big.mark = ","),

              "Age [Mean (SD)]" = cont(age_when_infected),

              "Sex (%)_sex" = "",
              "\tFemale_sex"       = bin(sex, "Female"),
              "\tMale_sex"         = bin(sex, "Male"),

              "BMI [kg/m2] [Mean (SD)]_body mass index" = cont(body_mass_index),

              "Index of multiple deprivation_index of multiple deprivation [Mean (SD)]" = cont(index_of_multiple_deprivation),

              "Ethnic background (%)_ethnic background" = "",
              "\tWhite_ethnic background" = bin(ethnic_background,"White"),
              "\tNon-white_ethnic background" = bin(ethnic_background,"Non-white"),

              "Smoking status (%)_smoking status" = "",
              "\tNever_smoking status" = bin(smoking_status, "Never"),
              "\tPrevious_smoking status" = bin(smoking_status, "Previous"),
              "\tCurrent_smoking status" = bin(smoking_status, "Current")) |>
    mutate_all(~as.character(.)) |>
    pivot_longer(everything(), names_to = "Risk factor", values_to = "UK Biobank dataset") |>
    mutate(Name = gsub(" ","_",gsub(".*_","",`Risk factor`))) |>
    mutate(`Risk factor` = gsub("_.*","",`Risk factor`)) |>
    rowwise() |>
    mutate("z" = if_else(
      Name == "other",
      list(cohort[["sex"]][!is.na(cohort[["sex"]])]),
      list(cohort[[Name]][!is.na(cohort[[Name]])]))) |>
    ungroup() |>
    group_by(Name) |>
    mutate(z = if_else(row_number() == 1, z, list(NULL))) |>
    ungroup()

  return(table)
}

tableOneBiomarkers <- function(cohort){
  t <- cohort |>
    summarise(across(
      .cols = -c("eid"),
      .fns = list(
        "Name" = cont),
      .names = "{.col} {.fn}"
    )) |>
    mutate_all(~as.character(.)) |>
    pivot_longer(everything(), names_to = "Risk factor", values_to = "UK Biobank dataset") |>
    mutate(Name = gsub(" .*","",`Risk factor`)) |>
    mutate(`Risk factor` = if_else(!grepl("Name",`Risk factor`), gsub(".* ","",`Risk factor`), `Risk factor`)) |>
    mutate(`Risk factor` = gsub("biomarker_","",`Risk factor`)) |>
    mutate(`Risk factor` = gsub("_", " ", `Risk factor`)) |>
    mutate(`Risk factor` = str_to_title(`Risk factor`)) |>
    mutate(`Risk factor` = gsub("Sd","SD",`Risk factor`)) |>
    mutate(`Risk factor` = gsub("Name","",`Risk factor`)) |>
    rowwise() |>
    mutate(z = list(cohort[[Name]][!is.na(cohort[[Name]])])) |>
    group_by(Name) |>
    mutate(z = if_else(row_number() == 1, z, list(NULL))) |>
    ungroup()
}

nameCategorised <- function(){
  tibble(
    "Risk factor" = c("ethnic_background", "body_mass_index", "body_mass_index",
                      "index_of_multiple_deprivation", "index_of_multiple_deprivation",
                      "index_of_multiple_deprivation", "index_of_multiple_deprivation",
                      "age_when_infected","age_when_infected","age_when_infected",
                      "age_when_infected","age_when_infected", "smoking_status",
                      "smoking_status", "sex"),
    "Variable"    = c("Non-white", "Overweight", "Obesity",
                      "Affluent", "Average", "Deprived", "Very deprived",
                      "55 to 59", "60 to 64","65 to 69","70 to 74", ">=75",  "Former", "Current",
                      "Male"
    )) |>
    group_by(`Risk factor`) |>
    mutate(n = row_number())

}

categoriseVariable <- function(x, variable){
  switch(variable,

         "body_mass_index" ={
           x <- x |>
             mutate(!!variable := case_when(
               .data[[variable]] < 25 ~ 0,
               .data[[variable]] >= 25 & .data[[variable]] < 30 ~ 1,
               .data[[variable]] >= 30 ~ 2
             ))
         },

         "index_of_multiple_deprivation" ={
           q <- quantile(x$index_of_multiple_deprivation, probs = c(0.2,0.4,0.6,0.8))

           x <- x |>
             mutate(!!variable := case_when(
               .data[[variable]] < q[[1]] ~ 0,
               .data[[variable]] >= q[[1]] & .data[[variable]] < q[[2]] ~ 1,
               .data[[variable]] >= q[[2]] & .data[[variable]] < q[[3]] ~ 2,
               .data[[variable]] >= q[[3]] & .data[[variable]] < q[[4]] ~ 3,
               .data[[variable]] >= q[[4]] ~ 4,
             ))
         },

         "age_when_infected" ={
           x <- x |>
             mutate("age" = age_when_infected) |>
             mutate(!!variable := case_when(
               .data[[variable]] < 55 ~ 0,
               .data[[variable]] >= 55 & .data[[variable]] < 60 ~ 1,
               .data[[variable]] >= 60 & .data[[variable]] < 65 ~ 2,
               .data[[variable]] >= 65 & .data[[variable]] < 70 ~ 3,
               .data[[variable]] >= 70 & .data[[variable]] < 75 ~ 4,
               .data[[variable]] >= 75 ~ 5
             ))
         }

  )

  x <- x |> mutate(!!variable := as.factor(.data[[variable]]))

  return(x)
}

tableOneCommorbidities <- function(commorbidities){
  commorbidities |>
    summarise(
      "Acquired immunodeficiency syndrome (AIDS) (%)_aids" =  bin(aids, "Cases"),
      "Asthma (%)_asthma" = bin(asthma, "Cases"),
      "Cancer (%)_cancer" =  bin(cancer, "Cases"),
      "Cancer - metastatic (%)_cancer metastatic" = bin(cancer_metastatic, "Cases"),
      "Cerebrovascular disease (%)_cerebrovascular disease" =  bin(cerebrovascular_disease, "Cases"),
      "Congestive heart failure (%)_congestive heart failure" =  bin(congestive_heart_failure, "Cases"),
      "Chronic kidney disease (%)_chronic kidney disease" =  bin(chronic_kidney_disease, "Cases"),
      "Chronic obstructive pulmonary disease (COPD) (%)_copd" =  bin(copd, "Cases"),
      "Dementia (%)_dementia" =  bin(dementia, "Cases"),
      "Diabetes (%)_diabetes" =  bin(diabetes, "Cases"),
      "Diabetes - organ damage (%)_diabetes organ damage" =  bin(diabetes_organ_damage, "Cases"),
      "Fracture (%)_fracture" =  bin(fracture, "Cases"),
      "Hemiplegia (%)_hemiplegia" = bin(hemiplegia, "Cases"),
      "Liver disease - mild (%)_liver disease mild" =  bin(liver_disease_mild, "Cases"),
      "Liver disease - moderate to severe (%)_liver disease moderate to severe" =  bin(liver_disease_moderate_to_severe, "Cases"),
      "Myocardial infarction (%)_myocardial infarction" =  bin(myocardial_infarction, "Cases"),
      "Peptic ulcer (%)_peptic ulcer" =  bin(peptic_ulcer, "Cases"),
      "Peripheral vascular disease (%)_peripheral vascular disease" =  bin(peripheral_vascular_disease, "Cases"),
      "Rheumatoid arthritis (%)_rheumatoid arthritis" = bin(rheumatoid_arthritis, "Cases")
    ) |>
    mutate_all(~as.character(.)) |>
    pivot_longer(everything(), names_to = "Risk factor", values_to = "UK Biobank dataset") |>
    mutate(Name = gsub(" ","_",gsub(".*_","",`Risk factor`))) |>
    mutate(`Risk factor` = gsub("_.*","",`Risk factor`))
}

addComorbidities <- function(cohort){

  comorbidities <- as_tibble(read.csv(paste0(dir_results,"/LoadCovariates_commorbidities.csv")))  |>
    dplyr::select(-c("X")) |>
    group_by(eid,value) |>
    mutate(episode_date = min(episode_date)) |>
    ungroup() |>
    distinct() |>
    pivot_wider(names_from = value, values_from = episode_date) |>
    mutate(across(-c("eid"), ~if_else(is.na(.), as.Date("3024-01-01"), as.Date(.)))) |>
    right_join(
      cohort |> dplyr::select("eid", "state", "specdate"),
      by = "eid"
    ) |>
    mutate(across(-c("eid", "state", "specdate"), ~if_else(. < specdate, 1, 0))) |>
    mutate(across(-c("eid", "state", "specdate"), ~if_else(is.na(.), 0, .))) |>
    mutate(fracture = if_else(
      ankle == 1 | ankle_open == 1 | elbow == 1 | elbow_open == 1 | femur_distal == 1 |
        femur_subtroch_shaft == 1 | femur_subtroch_shaft_history == 1 | knee == 1 |
        nhip_other == 1 | nhip_other_history == 1 | nhip_other_open == 1 | pelvis == 1 |
        pelvis_open == 1 | pelvis_spinalcord == 1 | radius_ulna_open == 1 | rib == 1 |
        shoulder == 1 | shoulder_open == 1 | spine == 1 | spine_history == 1 | spine_open == 1 |
        tibia == 1 | tibia_open == 1 | tibia_prox == 1 | wrist_forearm == 1 | wrist_forearm_open == 1 |
        hip == 1 | hip_open == 1 | foot == 1 | foot_open == 1 | lumbar_spine_pelvis == 1, 1, 0
    )) |>
    dplyr::select(-c("ankle", "ankle_open", "elbow", "elbow_open", "femur_distal", "femur_subtroch_shaft", "femur_subtroch_shaft_history",
              "knee", "nhip_other", "nhip_other_history", "nhip_other_open", "pelvis", "pelvis_open",
              "pelvis_spinalcord", "radius_ulna_open", "rib", "shoulder", "shoulder_open", "spine", "spine_history",
              "spine_open", "tibia", "tibia_open", "tibia_prox", "wrist_forearm", "wrist_forearm_open",
              "hip", "hip_open", "foot", "foot_open", "lumbar_spine_pelvis")) |>
    mutate(across(-c("eid","state","specdate"), ~if_else(.==1, "Cases","Controls"))) |>
    inner_join(
      cohort, by = c("eid", "state", "specdate")
    )

  return(comorbidities)
}

mergeAllCovariates <- function(cohort, baselineCharacteristics, biomarkers){
  cohort |>
    # Add baselinecharacteristics
    left_join(
      baselineCharacteristics,
      by = "eid"
    ) |>
    addCoding("sex") |>
    mutate("ethnic_background" = if_else(ethnic_background == 0, "White", "Non-white")) |>
    addCoding("smoking_status") |>
    # Calculate age when infected
    mutate(year_infected = year(specdate)) |>
    mutate(age_when_infected = year_infected - year_of_birth) |>
    dplyr::select("eid","specdate","state","ethnic_background","body_mass_index",
           "index_of_multiple_deprivation","sex","smoking_status","age_when_infected") |>
    # Add biomarkers
    left_join(
      biomarkers, by = "eid"
    ) |>
    # Add comorbidities
    addComorbidities()
}

tableOneStep1 <- function(x, baselineCharacteristics, biomarkers, name, name_cohort){

  cohort <- mergeAllCovariates(x, baselineCharacteristics, biomarkers)

  # Demographics table
  demographics_list <- list()
  biomarkers_list <- list()
  commorbidities_list <- list()

  for(ii in c(0:(length(name_cohort) -1))){
    cohort1 <- cohort |> filter(state == ii)

    cohort_demographics <- cohort1 |>
      tableOneDemographics() |>
      dplyr::select("Risk factor", "UK Biobank dataset") |>
      rename(!!name_cohort[ii+1] := "UK Biobank dataset")
    cohort_biomarkers <- cohort1 |>
      dplyr::select(all_of(colnames(biomarkers))) |>
      tableOneBiomarkers() |>
      dplyr::select("Risk factor", "UK Biobank dataset") |>
      rename(!!name_cohort[ii+1] := "UK Biobank dataset")
    cohort_commorbidities <- cohort1 |>
      tableOneCommorbidities() |>
      dplyr::select("Risk factor", "UK Biobank dataset") |>
      rename(!!name_cohort[ii+1] := "UK Biobank dataset")

    demographics_list[[ii+1]] <- cohort_demographics |>
      mutate(order = row_number())
    biomarkers_list[[ii+1]] <- cohort_biomarkers |>
      mutate(order = row_number())
    commorbidities_list[[ii+1]] <- cohort_commorbidities |>
      mutate(order = row_number())
  }

  if(length(demographics_list) == 2) {
    # Merge all the tables --- Long covid
    x_cohort <- tibble("Risk factor" = "Sociodemographic factors", "Controls" = " ", "order" = 0, "Cases" = " ") |>
      rbind(
        demographics_list[[1]] |>
          inner_join(demographics_list[[2]], by = c("order","Risk factor"))
      ) |>
      add_row(
        "Risk factor" = "Comorbidities [Cases (%)]", "Controls" = " ", "order" = 0, "Cases" = " ") |>
      rbind(
        commorbidities_list[[1]] |>
          inner_join(commorbidities_list[[2]], by = c("order", "Risk factor")) |>
          mutate(`Risk factor` = gsub(" \\(\\%\\)","",`Risk factor`))
      ) |>
      add_row("Risk factor" = "Biomarkers [Mean (SD)]","Controls" = " ","order" = 0,"Cases" = " ") |>
      rbind(
        biomarkers_list[[1]] |>
          inner_join(biomarkers_list[[2]], by = c("order", "Risk factor")) |>
          mutate(`Risk factor` = gsub(" \\(\\%\\)","",`Risk factor`))
      ) |>
      filter(!`Risk factor` %in% c("Counts [N (%)]","Missings", "Counts", "\t\tCounts", "\t\tQ05, Q25, Q50, Q75, Q95", "\t\tMissings")) |>
      dplyr::select(-c("order")) |>
      rename(!!paste0(name,"_Controls") := "Controls", !!paste0(name,"_Cases") := "Cases")
  } else {
    # Merge all the tables --- Long covid
    x_cohort <- tibble("Risk factor" = "Sociodemographic factors", "Cluster 1" = " ", "order" = 0, "Cluster 2" = " ", "Cluster 3" = " ") |>
      rbind(
        demographics_list[[1]] |>
          inner_join(demographics_list[[2]], by = c("order","Risk factor")) |>
          inner_join(demographics_list[[3]], by = c("order","Risk factor"))
      ) |>
      add_row(
        "Risk factor" = "Comorbidities [Cases (%)]", "Cluster 1" = " ", "order" = 0, "Cluster 2" = " ", "Cluster 3" = " ") |>
      rbind(
        commorbidities_list[[1]] |>
          inner_join(commorbidities_list[[2]], by = c("order", "Risk factor")) |>
          mutate(`Risk factor` = gsub(" \\(\\%\\)","",`Risk factor`)) |>
          inner_join(commorbidities_list[[3]], by = c("order", "Risk factor")) |>
          mutate(`Risk factor` = gsub(" \\(\\%\\)","",`Risk factor`))
      ) |>
      add_row("Risk factor" = "Biomarkers [Mean (SD)]","Cluster 1" = " ","order" = 0,"Cluster 2" = " ","Cluster 3" = " ") |>
      rbind(
        biomarkers_list[[1]] |>
          inner_join(biomarkers_list[[2]], by = c("order", "Risk factor")) |>
          mutate(`Risk factor` = gsub(" \\(\\%\\)","",`Risk factor`)) |>
          inner_join(biomarkers_list[[3]], by = c("order", "Risk factor")) |>
          mutate(`Risk factor` = gsub(" \\(\\%\\)","",`Risk factor`))
      ) |>
      filter(!`Risk factor` %in% c("Counts [N (%)]","Missings", "Counts", "\t\tCounts", "\t\tQ05, Q25, Q50, Q75, Q95", "\t\tMissings")) |>
      dplyr::select(-c("order")) |>
      rename(!!paste0(name,"_Cluster 1") := "Cluster 1", !!paste0(name,"_Cluster 2") := "Cluster 2", !!paste0(name,"_Cluster 3") := "Cluster 3")
  }
 

  return(x_cohort)
}

reverseCoding <- function(bd, variable){
  coding <- tibble(read.table(paste0("UKBiobank/",variable,"_coding.tsv"), header = TRUE, sep = "\t")) |>
    dplyr::select(!!variable := "meaning", "coding") |>
    mutate(!!variable := if_else(.data[[variable]] %in% c("Prefer not to answer", "Do not know", "Prefer not to say"), NA, .data[[variable]])) |>
    mutate(coding = if_else(is.na(.data[[variable]]), NA, coding)) |>
    distinct()

  bd <- bd |>
    inner_join(
      coding,
      by = variable
    ) |>
    dplyr::select(-!!variable) |>
    rename(!!variable := "coding")

  return(bd)
}

confint_adjust_mah <- function(object, parm, k, level = 0.95, method = "none") {
  # match method
  method <- match.arg(method, c("none", "bonferroni", "wh"))
  # estimated coefficients
  cf <- stats::coef(object)
  # estimated standard errors
  ses <- sqrt(diag(stats::vcov(object)))
  # select variables
  pnames <- names(ses)
  if (is.matrix(cf)) {
    cf <- stats::setNames(as.vector(cf), pnames)
  }
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  # alpha/2, 1-alpha/2
  a <- (1 - level)/2
  # number of intervals
  p <- object$rank
  if (method == "none") {
    fac <- stats::qt(c(a, 1 - a), object$df.residual)
    adj_level = max(1 - k * (1 - level), 0)
  } else if (method == "bonferroni") {
    fac <- stats::qt(c(a/k, 1 - a/k), object$df.residual)
    adj_level = level
  } else if (method == "wh") {
    fac <- c(-1, 1) * sqrt(p * stats::qf(level, df1 = p, df2 = object$df.residual))
    adj_level = level
  }

  # format returned object
  pct <- format_perc_api2lm(c(a, 1 - a), 3)
  ci <- array(NA_real_,
              dim = c(length(parm), 2L),
              dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci <- data.frame(term = row.names(ci),
                   estimate = cf[parm],
                   lwr = ci[,1],
                   upr = ci[,2])
  attributes(ci)$method <- method
  attributes(ci)$adj_level <- adj_level
  class(ci) <- c("confint_adjust", class(ci))
  ci
}

categoriseBiomarker <- function(x1){
  q <- quantile(x1$risk_factor, probs = c(0.2,0.4,0.6,0.8))

  x1 <- x1 |>
    mutate("risk_factor" = case_when(
      risk_factor < q[1] ~ 0,
      risk_factor < q[2] & risk_factor >= q[1] ~ 1,
      risk_factor < q[3] & risk_factor >= q[2] ~ 2,
      risk_factor < q[4] & risk_factor >= q[3] ~ 3,
      risk_factor >= q[4] ~ 4,
    )) |>
    mutate(risk_factor = as.factor(risk_factor))

  return(x1)
}

getForestPlot <- function(dt, dir_results, cohort_name, nam, l1, col1, col2, extra = ""){
  tm <- forest_theme(base_size = 10,
                     base_family = "Calibri",
                     refline_gp = gpar(lty = "solid"),
                     ci_pch = c(15, 15),
                     ci_col = c(col1,col2),
                     ci_Theight = 0,
                     legend_name = "",
                     legend_value = c("Crude", "Adjusted"),
                     legend_position = "top",
                     legend_gp = gpar(fontfamily = "Calibri", cex = 1),
                     colhead =  list(fg_params = list(fontsize = 12)))


  p <- forest(dt[,c(1,3,15,13,14)],
              est = list(dt$OR_crude,
                         dt$OR_adj),
              lower = list(dt$CI_LOW_crude,
                           dt$CI_LOW_adj),
              upper = list(dt$CI_HIGH_crude,
                           dt$CI_HIGH_adj),
              ci_column = 3,
              ref_line = 1,
              nudge_y = 0.4,
              x_trans = "log",
              xlim = c(0.5,l1),
              theme = tm,
              ticks_at = c(0.5, 0.75, 1, 1.5, seq(2,l1,1)),
              xlab = "OR"
  ) |>
    edit_plot(row = c(1:nrow(dt))[is.na(dt$OR_crude) & dt$Type != ""], gp = gpar(fontface = "bold.italic", fontsize = 11)) |>
    add_border(row = c(1:nrow(dt))[is.na(dt$OR_crude) & dt$Type != ""], where = "top", gp = gpar(cex = 0.25, col = "#C6C9C9")) |>
    add_border(part = "header", row = 1, where = "bottom")

  # p_sc <- get_scale(plot = p, width_wanted = 4, height_wanted = 0.15*nrow(dt), unit = "in")
  # ggplot2::ggsave(filename = if_else(analysis == "secondary",
  #                                    paste0(dir_results,"SecondaryAnalysis/",cohort_name,"_forestplot_main_",nam,extra,".png"),
  #                                    paste0(dir_results,"Figures/Main/",cohort_name,"_forestplot_main_",nam,extra,".png")),
  #                 plot = p,
  #                 dpi = 500,
  #                 units = "in",
  #                 width = 4, height = 5,
  #                 scale = p_sc
  #                 )

  return(p)
}

getUngroupedForestPlot <- function(dt, dir_results, cohort_name, nam, l1, col, ...){
  tm <- forest_theme(base_size = 10,
                     base_family = "Calibri",
                     refline_gp = gpar(lty = "solid"),
                     ci_pch = c(15),
                     ci_col = col,
                     ci_Theight = 0,
                     legend_gp = gpar(fontfamily = "Calibri", cex = 1),
                     colhead =  list(fg_params = list(fontsize = 12)))

  p <- forest(dt[,c(1,3,11,9,10)],
              est = dt$OR_adj,
              lower = dt$CI_LOW_adj,
              upper = dt$CI_HIGH_adj,
              ci_column = 3,
              ref_line = 1,
              nudge_y = 0,
              x_trans = "log",
              xlim = c(0.5,l1),
              theme = tm,
              ticks_at = c(0.5, 0.75, 1, 1.5, seq(2,l1,1)),
              xlab = "OR",
              ...
  ) |>
    edit_plot(row = c(1:nrow(dt))[is.na(dt$OR_adj) & dt$Type != ""], gp = gpar(fontface = "bold.italic", fontsize = 11)) |>
    add_border(row = c(1:nrow(dt))[is.na(dt$OR_adj) & dt$Type != ""], where = "top", gp = gpar(cex = 0.25, col = "#C6C9C9")) |>
    add_border(part = "header", row = 1, where = "bottom")

  return(p)
}

getManhattanPlot <- function(gwas, 
                             y_lim,
                             dot_size = 0.1,
                             colors = c("#92C5DE","#4393C3","#2166AC"),
                             chr_len = 22,
                             reduce_dataset = 0.01,
                             h_line_color = "red"
){
  # Reduce dataset
  gwas_top <- gwas %>% filter(LOG10P > 1)
  gwas_low <- gwas %>% filter(LOG10P <= 1) %>% sample_frac(reduce_dataset)
  
  gwas1 <- gwas_top %>% full_join(gwas_low) 
  
  # Manhattan plot -------
  don <- gwas1 %>%
    # Compute chromosome size
    group_by(CHROM) %>%
    summarise(chr_len = max(GENPOS)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwas1, ., by = c("CHROM" ="CHROM")) %>%
    # Add the cumulative position of each SNP
    arrange(CHROM, GENPOS) %>%
    mutate(GENPOScum = GENPOS+tot)
  
  # x axis
  axisdf = don %>% group_by(CHROM) %>% summarize(center=(max(GENPOScum) + min(GENPOScum) ) / 2 )
  
  plot1 <- ggplot(don, aes(x=GENPOScum, y=LOG10P)) +
    # Show all points
    geom_point(aes(color=as.factor(CHROM)), size = dot_size) +
    scale_color_manual(values = rep(colors, chr_len)) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center, expand = c(0.015,0.015)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,y_lim,2)) +     # remove space between plot area and x axis
    coord_cartesian(ylim = c(0,y_lim)) +
    # # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.text.x = element_text(margin = margin(t = 1), size = 9),
      axis.text.y  = element_text(size = 9),
      axis.title = element_text(size = 9),
    ) +
    #Plot a red horizontal line at 5e-8
    geom_hline(yintercept = -log10(5e-8), color = h_line_color, linewidth = 0.3) +
    labs(x = 'Chromosome', y = expression(-log[10](P)))
  
  return(plot1)
}

getQQPlot <- function(gwas,x_lim,y_lim,
                      color = "#2166AC",
                      dot_size = 0.1,
                      reduce_dataset = 0.01
){
  
  # QQ plot ----------------------------------------------------------------------
  n <- nrow(gwas)
  ci <- .95
  
  dat <- data.frame(
    observed = sort(gwas$LOG10P),
    expected = sort(-log10(ppoints(n))),
    clower   = sort(-log10(qbeta(p = (1 - ci) / 2, shape1 = seq(n), shape2 = rev(seq(n))))),
    cupper   = sort(-log10(qbeta(p = (1 + ci) / 2, shape1 = seq(n), shape2 = rev(seq(n)))))
  )
  
  # Reduce dataset size
  data_top <- dat %>% filter(observed > 1)
  data_low <- dat %>% filter(observed <= 1) %>% sample_frac(reduce_dataset)
  dat <- data_top %>% full_join(data_low)
  
  # Customize qqplot
  plot2 <- ggplot(dat, aes(x = expected, y = observed)) +
    scale_x_continuous(limits = c(0,x_lim), expand = c(0,0), breaks = seq(0,7,2)) + 
    scale_y_continuous(limits = c(0,y_lim), expand = c(0,0), breaks = seq(0,y_lim,4)) +
    geom_segment(data = . %>% filter(expected == max(expected)),
                 aes(x = 0, xend = x_lim, y = 0, yend = x_lim),
                 size = 0.25, color = "grey30", lineend = "round",alpha = 0.7) +
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    geom_point(color = color,  size = dot_size) +
    labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
         y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_rect(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 9)
    ) 
  
  return(plot2)
}
