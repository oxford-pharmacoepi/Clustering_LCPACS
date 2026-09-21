library(readxl)
library(dplyr)

x <- read_excel('Draft/DPA_Replication_VL_20260120.xlsx', sheet = 'DPA_SNPs_FinnGen')

x2 <- x %>%
  filter(!is.na(FinnGen_Analysis)) %>%
  filter(rsID != 'rs12976386') %>%
  mutate(
    p_value = 10^(-as.numeric(LOG10P)),
    tier = case_when(
      p_value <= 1e-3 ~ 'p\\leq 0.001',
      p_value <= 1e-2 ~ 'p\\leq 0.01',
      p_value < 0.05 ~ 'p<0.05',
      TRUE ~ 'not nominal'
    )
  ) %>%
  arrange(rsID, FinnGen_Analysis) %>%
  select(rsID, DPA_Analysis, NearestGene, FinnGen_Analysis, BETA, SE, LOG10P, p_value, tier)

esc <- function(v) gsub('_', '\\\\_', v, fixed = TRUE)

out <- character(nrow(x2))
for (i in seq_len(nrow(x2))) {
  r <- x2[i, ]
  out[i] <- paste0(
    r$rsID, ' & ',
    esc(as.character(r$DPA_Analysis)), ' & ',
    esc(as.character(r$NearestGene)), ' & ',
    esc(as.character(r$FinnGen_Analysis)), ' & ',
    sprintf('%.4f', as.numeric(r$BETA)), ' & ',
    sprintf('%.4f', as.numeric(r$SE)), ' & ',
    sprintf('%.4f', as.numeric(r$LOG10P)), ' & ',
    sprintf('%.2e', as.numeric(r$p_value)), ' & ',
    as.character(r$tier), ' \\\\'
  )
}

writeLines(out, 'Draft/fingen_full_rows.tex')
cat('Rows written:', length(out), '\n')
