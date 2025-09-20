# Get a big gene table from FUMA for all GWASes run

# Clust 1 vs 2&3
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vs23", "clust1vs23_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C1vs23","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C1vs23","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C1vs23","snps.txt"))

results1vs23 <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 1", controls = "Clust 2 + Clust 3",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 2 vs 1&3
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vs13", "clust2vs13_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C2vs13","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C2vs13","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C2vs13","snps.txt"))

results2vs13 <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 2", controls = "Clust 1 + Clust 3",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 3 vs 1&2
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust3vs12", "clust3vs12_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C3vs12","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C3vs12","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C3vs12","snps.txt"))

results3vs12 <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 3", controls = "Clust 1 + Clust 2",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 1 vs no
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vsno", "clust1vsno_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C1vsno","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C1vsno","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C1vsno","snps.txt"))

results1vsno <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 1", controls = "no COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 2 vs no
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vsno", "clust2vsno_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C2vsno","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C2vsno","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C2vsno","snps.txt"))

results2vsno <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 2", controls = "no COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 3 vs no
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust3vsno", "clust3vsno_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C3vsno","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C3vsno","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C3vsno","snps.txt"))

results3vsno <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 3", controls = "no COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 1 vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vsall", "clust1vsall_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C1vsall","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C1vsall","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C1vsall","snps.txt"))

results1vsall <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 1", controls = "Clust 2 + Clust 3 + no COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 2 vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vsall", "clust2vsall_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C2vsall","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C2vsall","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C2vsall","snps.txt"))

results2vsall <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 2", controls = "Clust 1 + Clust 3 + no COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 3 vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust3vsall", "clust3vsall_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"C3vsall","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"C3vsall","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"C3vsall","snps.txt"))

results3vsall <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 3", controls = "Clust 1 + Clust 2 + no COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Put it all together
results <- results1vs23 |>
  dplyr::union_all(results2vs13) |>
  dplyr::union_all(results3vs12) |>
  dplyr::union_all(results1vsno) |>
  dplyr::union_all(results2vsno) |>
  dplyr::union_all(results3vsno) |>
  dplyr::union_all(results1vsall) |>
  dplyr::union_all(results2vsall) |>
  dplyr::union_all(results3vsall)

# This can be used to inspect manually and annotate traits with clinical knowledge
results |>
  write_csv(file = paste0(dir_gwas_results,"/GWAS/", "snps_genes_annotated.csv"))

# Add info separately to non-lead SNPs
input1 <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vs23", "clust1vs23_FUMA.txt"))
input2 <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vs13", "clust2vs13_FUMA.txt"))
input3 <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vsno", "clust2vsno_FUMA.txt"))
input4 <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vsall", "clust2vsall_FUMA.txt"))

 refseq_map <- c(
   "1"="NC_000001.11","2"="NC_000002.12","3"="NC_000003.12",
   "4"="NC_000004.12","5"="NC_000005.10","6"="NC_000006.12",
   "7"="NC_000007.14","8"="NC_000008.11","9"="NC_000009.12",
   "10"="NC_000010.11","11"="NC_000011.10","12"="NC_000012.12",
   "13"="NC_000013.11","14"="NC_000014.9","15"="NC_000015.10",
   "16"="NC_000016.10","17"="NC_000017.11","18"="NC_000018.10",
   "19"="NC_000019.10","20"="NC_000020.11","21"="NC_000021.9",
   "22"="NC_000022.11","23"="NC_000023.11"
 )

non_lead_1_23 <- results |>
  dplyr::filter(is.na(chr), controls == "Clust 2 + Clust 3") |>
  dplyr::select(-c("chr", "pos", "P", "REF", "ALT", "OR", "SE")) |>
  tidyr::separate_rows(rsID, sep = ";") |>
  dplyr::left_join(input1, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA),
                chr = names(refseq_map)[match(CHR_REFSEQ, refseq_map)],
                pos = NA) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE", "gene", "type", "cases", "controls", "lead")

non_lead_2_13 <- results |>
  dplyr::filter(is.na(chr), controls == "Clust 1 + Clust 3") |>
  dplyr::select(-c("chr", "pos", "P", "REF", "ALT", "OR", "SE")) |>
  tidyr::separate_rows(rsID, sep = ";") |>
  dplyr::left_join(input2, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA),
                chr = names(refseq_map)[match(CHR_REFSEQ, refseq_map)],
                pos = NA) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE", "gene", "type", "cases", "controls", "lead")

non_lead_2_no <- results |>
  dplyr::filter(is.na(chr), controls == "no COVID-19") |>
  dplyr::select(-c("chr", "pos", "P", "REF", "ALT", "OR", "SE")) |>
  tidyr::separate_rows(rsID, sep = ";") |>
  dplyr::left_join(input3, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA),
                chr = names(refseq_map)[match(CHR_REFSEQ, refseq_map)],
                pos = NA) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE", "gene", "type", "cases", "controls", "lead")

non_lead_2_all <- results |>
  dplyr::filter(is.na(chr), controls == "Clust 1 + Clust 3 + no COVID-19") |>
  dplyr::select(-c("chr", "pos", "P", "REF", "ALT", "OR", "SE")) |>
  tidyr::separate_rows(rsID, sep = ";") |>
  dplyr::left_join(input4, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA),
                chr = names(refseq_map)[match(CHR_REFSEQ, refseq_map)],
                pos = NA) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE", "gene", "type", "cases", "controls", "lead")

# One can now inspect all these SNPs and add to the aforementioned big table if deemed necessary

# Second part, now get reported traits programatically
# Set the function 
genesymbol2gwas <- function(gene, column){
  url <- paste0(
    "https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:", gene,
    "&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true"
  )
  traits <- tryCatch({
    fread(url) |>
      pull(column) |>
      unique() |>
      paste(collapse = "; ")
  }, error = function(e) NA_character_) 
  
  return(traits)
}

# Use the function for all genes
results <- results %>%
  mutate(reported_traits = sapply(gene, genesymbol2gwas, column = "DISEASE/TRAIT")) 

# Good format
traits_long <- results %>%
  separate_rows(reported_traits, sep = ";\\s*") %>%
  filter(!is.na(reported_traits), reported_traits != "")

# Some plots of all traits in GWAS catalog
# Overall counts
trait_counts <- traits_long %>%
  count(reported_traits, sort = TRUE)

# Counts stratified by case/control
trait_counts_group <- traits_long %>%
  count(cases, controls, reported_traits, sort = TRUE)

# Overall barplot
p1 <- trait_counts %>%
  slice_max(n, n = 20) %>%
  ggplot(aes(x = n, y = fct_reorder(reported_traits, n))) +
  geom_col(fill = "steelblue") +
  labs(x = "Count", y = "Trait", title = "Top reported traits across all genes") +
  theme_minimal()
ggsave(paste0(dir_gwas_results,"Traits","top_traits_overall.png"), p1, width = 8, height = 6)

# Faceted by case
p2 <- trait_counts_group %>%
  group_by(cases) %>%
  slice_max(n, n = 15) %>%
  ggplot(aes(x = n, y = fct_reorder(reported_traits, n), fill = controls)) +
  geom_col() +
  facet_wrap(~cases, scales = "free_y") +
  labs(x = "Count", y = "Trait", title = "Top reported traits by case group") +
  theme_minimal()
ggsave(paste0(dir_gwas_results,"Traits","top_traits_by_case.png"), p2, width = 20, height = 6)

# Faceted by case and control
p3 <- trait_counts_group %>%
  mutate(group = paste0(cases, " vs ", controls)) |>
  group_by(group) %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = n, y = fct_reorder(reported_traits, n), fill = group)) +
  geom_col() +
  facet_wrap(~group, scales = "free_y") +
  labs(x = "Count", y = "Trait", title = "Top reported traits by case and control group") +
  theme(legend.position="none")
ggsave(paste0(dir_gwas_results,"Traits","top_traits_by_case_control.png"), p3, width = 20, height = 12)

# Wordcloud verall
png(paste0(dir_gwas_results,"Traits","wordcloud_overall.png"), width = 800, height = 600)
set.seed(123)
with(trait_counts, wordcloud(words = reported_traits, freq = n,
                             min.freq = 2, max.words = 200,
                             random.order = FALSE,
                             colors = brewer.pal(8, "Dark2")))
dev.off()

# By case/control grid
traits_long_count <- traits_long %>%
  mutate(study = paste0(cases," vs ", controls)) %>%
  group_by(cases, controls, study, reported_traits) %>%
  tally(name = "freq") %>%
  ungroup() %>%
  mutate(
    control_type = if_else(controls == "no COVID-19", "no COVID-19", "Other subtypes"),
    case_type = gsub("Clust", "Cluster", cases),
    control_type = if_else(grepl("no COVID-19", controls) & control_type == "Other subtypes",
                           "Other subtypes + no COVID-19", control_type)
  )

case_clusters <- c("Cluster 1", "Cluster 2", "Cluster 3")
control_types <- c("Other subtypes", "Other subtypes + no COVID-19", "no COVID-19")

n_rows <- length(case_clusters)
n_cols <- length(control_types)

png(paste0(dir_gwas_results,"Traits","wordcloud_grid.png"), width = 1800, height = 1200)
layout_matrix <- matrix(1:((n_rows+1)*(n_cols+1)), nrow = n_rows + 1, ncol = n_cols + 1, byrow = TRUE)
layout(layout_matrix, widths = c(2, rep(6, n_cols)), heights = c(1.5, rep(6, n_rows)))
par(mar = c(0,0,0,0))

# Top-left corner empty
plot.new()

# Column titles
for(ct in control_types){
  plot.new()
  text(0.5, 0.5, ct, cex = 3, font = 2)
}

# Row titles + wordclouds
for(case in case_clusters){
  plot.new()
  text(0.5, 0.5, case, srt = 90, cex = 3, font = 2)
  
  for(control in control_types){
    par(mar = c(0,0,1,0))
    dat <- subset(traits_long_count, case_type == case & control_type == control)
    
    if(nrow(dat) > 0){
      with(dat, wordcloud(words = reported_traits, freq = freq,
                          min.freq = 2, max.words = 40,
                          random.order = FALSE,
                          scale = c(4,0.5),
                          colors = brewer.pal(8, "Dark2")))
    } else plot.new()
  }
}
dev.off()

# Select top traits programatically
# Get unique traits first
unique_traits <- traits_long_count %>%
  distinct(reported_traits) 

# Needs this file in local (gwas catalog information for traits)
gwas_file <- "gwas_catalog_v1.0-associations_e115_r2025-09-15.tsv"

# Read the file
gwas <- fread(gwas_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Define the traits of interest
traits_of_interest <- unique_traits |> pull()

# Count number of associations per trait
trait_counts <- gwas %>%
  as_tibble() %>%
  filter(`DISEASE/TRAIT` %in% traits_of_interest) %>%
  group_by(`DISEASE/TRAIT`) %>%
  tally() |>
  rename("reported_traits" = "DISEASE/TRAIT")

# Join back to main table
traits_long_count <- traits_long_count %>%
  left_join(trait_counts, by = "reported_traits")

traits_long_count |>
  dplyr::select(c("case_type", "control_type", "reported_traits", "n_study" = "freq", "n_catalog" = "n")) |>
  arrange(-n_study) |>
  write_csv(paste0(dir_gwas_results,"Traits","traits_GWAS_long.csv"))

# total counts per study
study_totals <- traits_long_count %>%
  group_by(study) %>%
  summarise(N_study = sum(freq))

# Compute frequency by total GWAS catalog counts
df <- traits_long_count %>%
  left_join(study_totals, by = "study") %>%
  mutate(
    freq_weight = freq/n
  )

# Select top traits per study
top_traits <- df %>%
  group_by(study) %>%
  arrange(freq_weight) %>%
  slice_max(freq_weight, n = 10) %>%
  ungroup()

# Plot
p <- ggplot(top_traits, aes(x = freq_weight, y = reorder(reported_traits, freq_weight),
                            size = freq, color = freq_weight)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~study, scales = "free_y") +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Counts weighted by GWAS catalog", y = "Trait",
       title = "Top relevant traits per study (adjusted for GWAS catalog)", 
       size = "Count in study") +
  theme_bw()
ggsave(paste0(dir_gwas_results,"Traits","top_traits_enrichment.png"), width = 30, height = 8, dpi = 300)

top_percent <- 0.1 # chan change for other % cutoff
catalog_cutoff <- quantile(df$n, probs = 1 - top_percent, na.rm = TRUE)

# Select top traits per study after lowest percentile cutoff
top_traits_g1 <- df %>%
  dplyr::filter(n >= catalog_cutoff) %>%
  group_by(study) %>%
  arrange(freq_weight) %>%
  slice_max(freq_weight, n = 10) %>%
  ungroup()

# Plot
p <- ggplot(top_traits_g1, aes(x = freq_weight, y = reorder(reported_traits, freq_weight),
                            size = freq, color = freq_weight)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~study, scales = "free_y") +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Counts weighted by GWAS catalog", y = "Trait",
       title = "Top relevant traits per study (adjusted for GWAS catalog)", 
       size = "Count in study") +
  theme_bw()
ggsave(paste0(dir_gwas_results,"Traits","top_traits_enrichment_cutoff.png"), width = 30, height = 8, dpi = 300)

# Binomial cumulative procedure
total_catalog_counts <- sum(trait_counts$n)

df <- traits_long_count %>%
  group_by(study) %>%
  mutate(
    n_study = sum(freq),
    p_expected = n / total_catalog_counts,
    # cumulative binomial probability
    p_binom = 1 - pbinom(freq - 1, n_study, p_expected),
    impressiveness = -log10(p_binom + 1e-10)  # log-transform for visualisation
  ) %>%
  ungroup()

# select top traits per study
top_traits <- df %>%
  group_by(study) %>%
  slice_max(impressiveness, n = 10) %>%
  ungroup()

# plot
ggplot(top_traits, aes(x = impressiveness, y = reorder(reported_traits, impressiveness),
                       size = freq, color = impressiveness)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~study, scales = "free_y") +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "-log10(cumulative binomial p)", y = "Trait",
       title = "Top impressive traits per study (adjusted for catalog popularity)",
       size = "Count in study") +
  theme_bw()

ggsave(paste0(dir_gwas_results,"Traits","top_traits_impressiveness.png"), width = 30, height = 8, dpi = 300)

# With small counts filter
min_count <- 10  # or another number reasonable
df_filtered <- df %>% 
  filter(n >= min_count)

top_traits <- df_filtered %>%
  group_by(study) %>%
  slice_max(impressiveness, n = 10) %>%
  ungroup()

ggplot(top_traits, aes(x = impressiveness, y = reorder(reported_traits, impressiveness),
                       size = freq, color = impressiveness)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~study, scales = "free_y") +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Impressiveness (-log10 binomial p)", y = "Trait",
       title = "Top relevant traits per study (adjusted for GWAS catalog)",
       size = "Count in study") +
  theme_bw()
ggsave(paste0(dir_gwas_results,"Traits","top_traits_impressiveness_cutoff.png"), width = 30, height = 8, dpi = 300)
