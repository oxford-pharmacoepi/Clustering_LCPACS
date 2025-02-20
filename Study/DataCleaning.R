# Load data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"/LoadCovariates_baselineCharacteristics.csv")))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"/LoadCovariates_biomarkers.csv")))
comorbidities <- as_tibble(read.csv(paste0(dir_results,"/LoadCovariates_commorbidities.csv")))

# Merge tables ----
baselineCharacteristics <- baselineCharacteristics |> dplyr::select(-c("X")) |> mutate(smoking_status = if_else(smoking_status == -3, NA, smoking_status))
biomarkers              <- biomarkers |> dplyr::select(-c("X", "rheumatoid_factor", "oestradiol"))
comorbidities           <- comorbidities |> loadComorbidities()

data <- baselineCharacteristics |>
  inner_join(biomarkers, by = "eid") |>
  inner_join(comorbidities, by = "eid")

# Create table to record missing and outliers ----
cleaning_table <- tibble("variables" = colnames(data),
                         "outliers"   = 0, "outliers_pc" = 0)

data_transformed <- data |> dplyr::select(-c("eid"))

# Multiple imputation ----
x <- mice(data_transformed, seed = 1, m = 5, print = TRUE)
data_transformed1 <- complete(x) |> as_tibble()

cleaning_table$missings <- c(colSums(is.na(data)))
cleaning_table$missings_pc <- c(colSums(is.na(data))/nrow(data)*100)

# Remove outliers ----
data_transformed2 <- data_transformed1 |>
  mutate(across(
    .cols = -c("sex", "year_of_birth","ethnic_background","smoking_status", any_of(colnames(comorbidities))),
    .fns  = ~if_else(. > (mean(., na.rm = TRUE) + 3*sd(., na.rm = TRUE)),
                     mean(., na.rm = TRUE) + 3*sd(., na.rm = TRUE),
                     .)
  ))

for (i in seq_len(length(colnames(data_transformed1)))){
  cleaning_table$outliers[i+1]    <- length(setdiff(data_transformed1[,i], data_transformed2[,i]) |> pull())
  cleaning_table$outliers_pc[i+1] <- (length(setdiff(data_transformed1[,i], data_transformed2[,i]) |> pull()))/nrow(data_transformed1)*100
}

# Z transformation ----
data_transformed3 <- data_transformed2 |>
  mutate(across(
    -c("sex", "year_of_birth","ethnic_background","smoking_status", "body_mass_index", "index_of_multiple_deprivation", any_of(colnames(comorbidities))),
    ~ (.-mean(.))/sd(.))) |>
  mutate(eid = data$eid) |>
  relocate(eid)

# Further transformations ----
data_transformed4 <- data_transformed3 |>
  # Ethnic background
  mutate("ethnic_background" = if_else(ethnic_background %in% c(1, 1001, 1002, 1003), 0, 1))

# Tidy datasets ----
data_transformed4 |>
  dplyr::select("eid", all_of(colnames(baselineCharacteristics))) |>
  write.csv(paste0(dir_results, "/CleanData_baselineCharacteristics.csv"))

data_transformed4 |>
  dplyr::select("eid", all_of(colnames(biomarkers))) |>
  write.csv(paste0(dir_results, "/CleanData_biomarkers.csv"))

data_transformed4 |>
  dplyr::select("eid", all_of(colnames(comorbidities))) |>
  write.csv(paste0(dir_results, "/CleanData_comorbidities.csv"))

# Table of missing values and outliers -----
tibble(variables = c(colnames(baselineCharacteristics), colnames(biomarkers), colnames(comorbidities))) |>
  filter(!variables %in% c("X", "eid")) |>
  inner_join(
    cleaning_table, by = "variables"
  ) |>
  mutate(variables = gsub("biomarker_","",variables)) |>
  mutate(variables = str_to_sentence(gsub("_"," ",variables))) |>
  filter(variables != "Eid") |>
  mutate(`Number of outliers (%)` = paste0(formatC(outliers, big.mark = ",")," (",round(outliers_pc,2),")")) |>
  mutate(`Number of missings (%)` = paste0(formatC(missings, big.mark = ",",format = "d")," (",round(missings_pc,2),")"))  |>
  dplyr::select(-c("outliers", "outliers_pc", "missings", "missings_pc")) |>
  rename("Risk factors" = "variables") |>
  flextable() |>
  width(j = 1, width = 6.25, unit = "cm") |>
  width(j = c(2,3), width = 3, unit = "cm") |>
  bold(i = 1, part = "header") |>
  align(j = c(2,3), align = "center", part = "all") |>
  save_as_docx(path = paste0(dir_results,"/Tables/cleaning_table.docx"))


rm(list = c("baselineCharacteristics", "biomarkers", "comorbidities", "data", "cleaning_table",
            "data_transformed", "x", "data_transformed1", "data_transformed2", "i", "data_transformed3",
            "data_transformed4"))
