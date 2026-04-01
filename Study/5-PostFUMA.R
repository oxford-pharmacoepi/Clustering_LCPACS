# Get a big gene table from FUMA for all GWASes run

# Clust 1 vs 2&3
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vs23/", "clust1vs23_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C1vs23/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C1vs23/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C1vs23/","snps.txt"))

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
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vs13/", "clust2vs13_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C2vs13/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C2vs13/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C2vs13/","snps.txt"))

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
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust3vs12/", "clust3vs12_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C3vs12/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C3vs12/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C3vs12/","snps.txt"))

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
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vsno/", "clust1vsno_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C1vsno/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C1vsno/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C1vsno/","snps.txt"))

results1vsno <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 1", controls = "COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 2 vs no
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vsno/", "clust2vsno_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C2vsno/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C2vsno/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C2vsno/","snps.txt"))

results2vsno <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 2", controls = "COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 3 vs no
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust3vsno/", "clust3vsno_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C3vsno/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C3vsno/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C3vsno/","snps.txt"))

results3vsno <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 3", controls = "COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 1 vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust1vsall/", "clust1vsall_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C1vsall/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C1vsall/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C1vsall/","snps.txt"))

results1vsall <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 1", controls = "Clust 2 + Clust 3 + COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 2 vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust2vsall/", "clust2vsall_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C2vsall/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C2vsall/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C2vsall/","snps.txt"))

results2vsall <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 2", controls = "Clust 1 + Clust 3 + COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

# Clust 3 vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "Clust3vsall/", "clust3vsall_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/C3vsall/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/C3vsall/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/C3vsall/","snps.txt"))

results3vsall <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "Clust 3", controls = "Clust 1 + Clust 2 + COVID-19",
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
  dplyr::filter(is.na(chr), controls == "COVID-19") |>
  dplyr::select(-c("chr", "pos", "P", "REF", "ALT", "OR", "SE")) |>
  tidyr::separate_rows(rsID, sep = ";") |>
  dplyr::left_join(input3, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA),
                chr = names(refseq_map)[match(CHR_REFSEQ, refseq_map)],
                pos = NA) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE", "gene", "type", "cases", "controls", "lead")

non_lead_2_all <- results |>
  dplyr::filter(is.na(chr), controls == "Clust 1 + Clust 3 + COVID-19") |>
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
genesymbol2gwas_pairs <- function(gene,
                                  cols = c("DISEASE/TRAIT", "MAPPED_TRAIT")) {
  url <- paste0(
    "https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:", gene,
    "&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true"
  )
  dt <- tryCatch(fread(url), error = function(e) NULL)
  if (is.null(dt)) {
    return(tibble(reported = character(), mapped = character()))
  }
  
  out <- as_tibble(dt) %>%
    transmute(
      reported = .data[[cols[1]]],
      mapped   = .data[[cols[2]]]
    ) %>%
    mutate(
      reported = na_if(str_squish(as.character(reported)), ""),
      mapped   = na_if(str_squish(as.character(mapped)), "")
    ) %>%
    filter(!is.na(reported), !is.na(mapped)) %>%
    # remove exact duplicate PAIRS (not independently)
    distinct(reported, mapped, .keep_all = TRUE)
  
  out
}

# Build pairwise lists per gene, then collapse each side AFTER pairing
pairs_by_gene <- results %>%
  distinct(gene) %>%
  mutate(pairs = lapply(gene, genesymbol2gwas_pairs))

results <- results %>%
  left_join(
    pairs_by_gene %>%
      mutate(
        reported_traits = map_chr(
          pairs, ~ paste(na.omit(.x$reported), collapse = "&&& ")
        ),
        mapped_traits = map_chr(
          pairs, ~ paste(na.omit(.x$mapped), collapse = "&&& ")
        )
      ) %>%
      dplyr::select(gene, reported_traits, mapped_traits),
    by = "gene"
  )

# Inspection on mapping traits
traits_long_aligned <- pairs_by_gene %>%
  dplyr::select(gene, pairs) %>%
  unnest(pairs) %>%              # columns: gene, reported, mapped
  mutate(k = row_number(), .by = gene)

# Good format
rep_long <- results %>%
  transmute(
    gene,
    reported = str_split(coalesce(reported_traits, ""), "\\s*&&&\\s*")
  ) %>%
  unnest_longer(reported, indices_include = TRUE, indices_to = "k") %>%
  mutate(reported = na_if(str_trim(reported), ""))

map_long <- results %>%
  transmute(
    gene,
    mapped = str_split(coalesce(mapped_traits, ""), "\\s*&&&\\s*")
  ) %>%
  unnest_longer(mapped, indices_include = TRUE, indices_to = "k") %>%
  mutate(mapped = na_if(str_trim(mapped), ""))

# 2) Align by (gene, position)
traits_aligned <- full_join(rep_long, map_long, by = c("gene","k")) %>%
  # drop rows where both sides are empty
  filter(!(is.na(reported) & is.na(mapped)))

traits_long <- results |>
  dplyr::select(-c("reported_traits", "mapped_traits")) |>
  dplyr::left_join(traits_aligned |>
                     dplyr::select(c("gene", "reported", "mapped"))) |>
  dplyr::filter(!is.na(reported)) |>
  dplyr::filter(!is.na(mapped)) |>
  dplyr::distinct()

# Some plots of all traits in GWAS catalog
# Overall counts
trait_counts <- traits_long %>%
  dplyr::count(reported, sort = TRUE)

trait_counts_mapped <- traits_long %>%
  dplyr::count(mapped, sort = TRUE)

# Counts stratified by case/control
trait_counts_group <- traits_long %>%
  dplyr::count(cases, controls, reported, sort = TRUE)

# Counts stratified by case/control
trait_counts_group_mapped <- traits_long %>%
  dplyr::count(cases, controls, mapped, sort = TRUE)

# Overall barplot
p1 <- trait_counts %>%
  slice_max(n, n = 20) %>%
  ggplot(aes(x = n, y = fct_reorder(reported, n))) +
  geom_col(fill = "steelblue") +
  labs(x = "Count", y = "Trait", title = "Top reported traits across all genes") +
  theme_minimal()
ggsave(paste0(dir_gwas_results,"/Traits/","top_traits_overall.png"), p1, width = 8, height = 6)

p1 <- trait_counts_mapped %>%
  slice_max(n, n = 20) %>%
  ggplot(aes(x = n, y = fct_reorder(mapped, n))) +
  geom_col(fill = "steelblue") +
  labs(x = "Count", y = "Trait", title = "Top mapped traits across all genes") +
  theme_minimal()
ggsave(paste0(dir_gwas_results,"/Traits/","top_traits_mapped_overall.png"), p1, width = 8, height = 6)

# Faceted by case
p2 <- trait_counts_group %>%
  group_by(cases) %>%
  slice_max(n, n = 15) %>%
  ggplot(aes(x = n, y = fct_reorder(reported, n), fill = controls)) +
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
  ggplot(aes(x = n, y = fct_reorder(reported, n), fill = group)) +
  geom_col() +
  facet_wrap(~group, scales = "free_y") +
  labs(x = "Count", y = "Trait", title = "Top reported traits by case and control group") +
  theme(legend.position="none")
ggsave(paste0(dir_gwas_results,"Traits","top_traits_by_case_control.png"), p3, width = 20, height = 12)

# Wordcloud overall
png(paste0(dir_gwas_results,"/Traits/","wordcloud_overall.png"), width = 800, height = 600)
set.seed(123)
with(trait_counts, wordcloud(words = reported, freq = n,
                             min.freq = 2, max.words = 200,
                             random.order = FALSE,
                             colors = brewer.pal(8, "Dark2")))
dev.off()

# By case/control grid
traits_long_count <- traits_long %>%
  mutate(study = paste0(cases," vs ", controls)) %>%
  group_by(cases, controls, study, reported) %>%
  tally(name = "freq") %>%
  ungroup() %>%
  mutate(
    control_type = if_else(controls == "COVID-19", "COVID-19", "Other subtypes"),
    case_type = gsub("Clust", "Cluster", cases),
    control_type = if_else(grepl("COVID-19", controls) & control_type == "Other subtypes",
                           "Other subtypes + COVID-19", control_type)
  )

case_clusters <- c("Cluster 1", "Cluster 2", "Cluster 3")
control_types <- c("Other subtypes", "Other subtypes + COVID-19", "COVID-19")

n_rows <- length(case_clusters)
n_cols <- length(control_types)

png(here(dir_gwas_results,"Traits","wordcloud_grid.png"), width = 1600, height = 1000)
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
    par(mar = c(0,0,0,0))
    dat <- subset(traits_long_count, case_type == case & control_type == control)
    
    if(nrow(dat) > 0){
      with(dat, wordcloud(words = reported, freq = freq,
                          min.freq = 1, max.words = 40,
                          random.order = FALSE,
                          scale = c(4,0.5),
                          colors = brewer.pal(8, "Dark2")))
    } else plot.new()
  }
}
dev.off()


# Mapped
traits_long_count_mapped <- traits_long %>%
  mutate(study = paste0(cases," vs ", controls)) %>%
  group_by(cases, controls, study, mapped) %>%
  tally(name = "freq") %>%
  ungroup() %>%
  mutate(
    control_type = if_else(controls == "COVID-19", "COVID-19", "Other subtypes"),
    case_type = gsub("Clust", "Cluster", cases),
    control_type = if_else(grepl("COVID-19", controls) & control_type == "Other subtypes",
                           "Other subtypes + COVID-19", control_type)
  )
png(here(dir_gwas_results,"Traits","wordcloud_grid_mapped.png"), width = 1800, height = 1200)
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
    dat <- subset(traits_long_count_mapped, case_type == case & control_type == control)
    
    if(nrow(dat) > 0){
      with(dat, wordcloud(words = mapped, freq = freq,
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
  distinct(reported) 

# Needs this file in local (gwas catalog information for traits)
gwas_file <- "gwas_catalog_v1.0-associations_e115_r2025-09-15.tsv"

# Read the file
gwas <- fread(gwas_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Define the traits of interest
traits_of_interest <- unique_traits |> pull()

# Binomial cumulative procedure
# First restrict the universe to disease-related traits
catalog <- gwas %>% as_tibble()
colnames(catalog) <- gsub("\\s+", "_", colnames(catalog))

# Check how many of the traits in my tibble of interest appear in the catalog file
traits_catalog <- catalog |>
  dplyr::pull(`DISEASE/TRAIT`) |>
  unique()
sum(traits_of_interest %in% traits_catalog)/length(traits_of_interest)
# More than 99%

# ===============================
# Disease / Measurement vocabularies (clean, deduplicated, updated)
# ===============================

# --- Core disease keywords ---
disease_core <- c(
  # General / categories
  "disease","disorder","syndrome","cancer","neoplasm","tumou?r","carcinoma",
  "adenocarcinoma","sarcoma","glioma","glioblastoma","meningioma","melanoma",
  "leukemi","lymphom","myeloma","neuroblastoma","pituitary adenoma",
  "malignant","benign neoplasm", "insomnia", "sleep", "alcohol",
  
  # Infection / inflammatory / autoimmune
  "infection","infectious","sepsis","pneumonia","meningitis","hepatitis","tubercul",
  "influenza","cytomegalovirus","epstein","herpes","hiv","aids","malaria",
  "leprosy","measles","mumps","rubella","chickenpox","shingles",
  "scarlet fever","rheumatic fever","dysentery","typhoid","enteric fever",
  "buruli ulcer","leishmaniasis","mononucleosis","bacteremia","bacteraemia",
  "bullous pemphigoid","pemphigus","sarcoidosis","vasculitis","thyroiditis",
  "autoimmune","inflammatory","arthritis","osteoarth","spondylo","psoriasis",
  "eczema","dermatitis","ibd","crohn","ulcerative colitis","rhinosinusitis",
  "urticaria","atopy","asthma","allergy","eosinophilic","myasthenia","myositis",
  "scleroderma","systemic sclerosis","lupus","wegener","ankylosing spondylitis",
  "autoimmune hemolytic anemia","antiphospholipid","selective iga deficiency",
  
  # Haematologic / metabolic
  "anemia","thalassemia","sickle cell","polycythemia","hemochromatosis",
  "iron deficiency","iron overload","thrombocytopenia","thrombocytosis",
  "thrombophilia","cytopenia","heparin-induced","coagulation","hyperuricemia",
  "diabetes","gout","metabolic syndrome","obes","dyslip","hyperlip","fatty liver",
  "steatohepat","nash","hypoglycemia","hyperglycemia","monoclonal gammopathy",
  
  # Cardiovascular / cerebrovascular
  "cardio","coronary","myocard","heart failure","myocardial infarction","angina",
  "arrhythmia","atrial fibrillation","ventricular ectopy","tachycardia","bradycardia",
  "stroke","ischemic stroke","intracerebral hemorrhage","subarachnoid hemorrhage",
  "thrombosis","thromboembol","venous thromboembolism","aneurysm","amyloid angiopathy",
  "arteriosclerosis","atherosclerosis","endocarditis","pericarditis","vasculopathy",
  "aortic valve","aortic root","aortic stiffness","hypertension","hypotension",
  "orthostatic hypotension","varicose veins","pulmonary embolism",
  "intracranial artery stenosis","vascular disease","fibrosis","arterial","venous",
  
  # Respiratory
  "respiratory","bronchiectasis","bronchopulmonary dysplasia","pulmonary fibrosis",
  "emphysema","copd","asthma","pneumothorax","bronchitis","sleep apnea",
  "obstructive sleep apnea","airflow obstruction","chronic mucus","interstitial lung",
  "lymphangioleiomyomatosis","pneumococcal",
  
  # Renal / urological
  "renal","kidney","ckd","nephrop","nephritis","pyelonephritis","glomerulonephritis",
  "glomerulosclerosis","nephropathy","acute kidney injury","urinary incontinence",
  "neuropathic bladder","urolithiasis","nephrolithiasis","kidney stone",
  "vesicoureteric reflux","hematuria","proteinuria","cystitis",
  
  # Gastrointestinal / hepatic
  "hepatic","liver injury","cirrhosis","biliary","gastritis","gastric ulcer",
  "duodenal ulcer","cholecystitis","cholelithiasis","gallbladder polyp",
  "gastric polyp","esophageal varix","enterocolitis","colitis","crohn","ulcerative colitis",
  "barrett'?s (oesophagus|esophagus)","hepatitis","pancreatitis","gastric atrophy",
  "esophageal ulcer","polyp of stomach and duodenum",
  
  # Neurological / psychiatric
  "neurolog","alzheimer","parkinson","epilep","multiple sclerosis","\\bms\\b",
  "amyotrophic lateral sclerosis","als","migraine","dementia","cluster headache",
  "headache","neuropathy","dystonia","psychosis","bipolar","schizophren",
  "depress","anxiety","suicid","ocd","obsessive[- ]compulsive","adhd","autism",
  "pathological gambling","tardive dyskinesia","narcolepsy","essential tremor",
  "motion sickness","vestibular neuritis","vertigo","tinnitus","delirium",
  "developmental stuttering","hallucination","impulsivity","hyperactivity",
  "inattentive symptoms","bulimia","anorexia nervosa","substance dependence",
  "opioid dependence","cocaine dependence","cannabis dependence","heroin dependence",
  "stimulant dependence",
  
  # Musculoskeletal / connective tissue
  "arthritis","arthropathy","arthrosis","osteoporosis","osteonecrosis","fracture",
  "myositis","fibromyalgia","scoliosis","spinal stenosis","tendinopathy","ossification",
  "rotator cuff","clubfoot","trigger finger","frozen shoulder","spinal canal stenosis",
  "ankylosing","spondylitis","achilles tendinopathy",
  
  # Endocrine / reproductive
  "thyroid","goiter","hyperthyroid","hypothyroid","polycystic","pcos","menstrual",
  "fibroids","premature ovarian failure","erectile dysfunction","uterine","infertility",
  "azoospermia","androgen","menopause","puberty","endometriosis","breast","pituitary",
  "ovarian cyst","mastopathy","endocrine","thyrotoxic hypokalemic periodic paralysis",
  "hyperemesis gravidarum","ovarian reserve","male infertility",
  
  # Dermatologic
  "acne","eczema","dermatitis","psoriasis","vitiligo","alopecia","pemphigus",
  "urticaria","rosacea","keloid","melanoma","scar","pigmentation","skin cancer",
  "cutaneous lupus","bullous pemphigoid","atopic","plantar warts","dry skin",
  "lichen planus","café[- ]au[- ]lait macule","age spots","excessive sweating",
  
  # ENT / ocular
  "otitis","sinusitis","hearing loss","deafness","tonsillitis","myringotomy",
  "conjunctivitis","keratitis","uveitis","glaucoma","retinopathy","macular degeneration",
  "iridocyclitis","iritis","otosclerosis","nasal polyps","colorectal polyp","colonic polyp",
  
  # Rare / congenital / developmental
  "craniosynostosis","craniofacial","dysplasia","malformation","syndrome","situs",
  "trisomy","microcephaly","polydactyly","fibrodysplasia","cleft (lip|palate)",
  "orofacial cleft","hypospadias","classic bladder exstrophy",
  "tetralogy of fallot","infantile hypertrophic pyloric stenosis",
  "neural tube defect","agenesis",
  
  # Other specific entities
  "amyloidosis","amyloid angiopath","amyloid beta deposition","amyloid status",
  "neuralgia","rheumatic fever","gastroenteritis","cold sores","covid","covid-19",
  "sars-cov-2","monoclonal gammopathy of undetermined significance","mgus",
  "mastocytosis","myeloid clonal hematopoiesis","clonal hematopoiesis",
  "systemic mastocytosis","lymphangioleiomyomatosis", "concussion", "septic",
  "hernia", "stress", "pain", "seizures", "smell", "taste",
  "icd10", "phecode"
)

# --- Drug treatment / toxicity / survival endpoints (keep if disease context) ---
disease_treatment <- c(
  "response to","treatment response","drug[- ]induced","toxicity","adverse response",
  "adverse drug reaction","side effect","vaccine-related adverse event","subjective response",
  "objective response","placebo response","drug hypersensitivity","treatment outcome",
  "progression[- ]?free survival","overall survival","mortality","prognosis","time to event",
  "maintenance dose","dosage","metformin","lithium","clopidogrel","warfarin","acenocoumarol",
  "asparaginase","thiopurine","interferon","methylphenidate","paliperidone","duloxetine",
  "radiotherapy","chemotherapy","paclitaxel","trastuzumab","carboplatin","cisplatin",
  "tamoxifen","fenofibrate","statin","ssri","bronchodilator","glucocorticoid","inhaled glucocorticoid",
  "antineoplastic","immune checkpoint inhibitor","nevirapine","capecitabine","epirubicin"
)

# --- Disease progression patterns (keep) ---
disease_progression <- c(
  "disease progression","progression to","progression in","disease relapse","remission"
)

# --- Measurement / imaging / behavioral / biomarker keywords (drop) ---
measurement_like <- c(
  # Quantitative / derived variables
  "level","levels","concentration","ratio","index","score","z[ -]?score","percent(age|ile)?",
  "count","percentage","fraction","mean","variance","variability","change over time", "area",
  "width", "number of", "depth",
  
  # Biomarkers / lab assays
  "hdl|ldl|cholesterol|triglyceride|glucose|insulin|hba1c|egfr|creatinine|bilirubin|ferritin|cystatin|transaminase|urea",
  "c[- ]?reactive protein|fibrinogen|thrombin[- ]?antithrombin|nt[- ]?probnp|uric acid|albumin",
  "protein level","proteomic","glycosylation","cytokine","il-","ifn","tnf","chemokine",
  "biochemical measure","metabolite","lipidomic","metabolomic","metabolite", "measurement",
  
  # Imaging / morphometry / IDPs
  "mri|fmri|ct|scan|imaging","brain volume","cortical","white matter","grey matter",
  "ventricle","hippocamp","nucleus","gyrus","morpholog","thickness","surface area","volume",
  "optic disc","corneal curvature","bone mineral density","bone (content|size)",
  "idp","fast rois","aseg","diffusion","tract","t2star","brain shape","amyloid status",
  "tau burden","cerebral blood flow","connectivity","cytoarchitecture",
  
  # ECG / cardiac intervals & dimensions
  "\\bqt\\b|qt interval|pr interval|qrs complex|rr interval|jt interval|lvef",
  "left ventricular (mass|dimension)","left atrial","aortic root size","p wave terminal force",
  
  # Anthropometrics / body composition
  "height|weight|bmi|circumference|waist|hip|body fat|grip strength|lean mass|fat mass|trunk fat",
  "adiposity|anthropometric|sarcopenia|skeletal muscle mass|visceral fat",
  "appendicular lean mass|lean body mass|bone density|adipose tissue",
  
  # Behavioral / social / psychometric
  "intelligence","cognitive","memory","learning","language","reading","executive","reaction time",
  "processing speed","reasoning","personality","temperament","neuroticism","mood instability",
  "psychopathology","emotion","empathy","subjective well-being","hedonic","eudaimonic",
  "risk[- ]?taking","alcohol consumption","drinking","smok","nicotine","caffeine","addiction","habit",
  "education","income","occupation","self-rated health","job-related exhaustion","social","lonely",
  
  # Appearance / pigmentation / morphology
  "eye color","iris","freckles","hair (color|shape|texture|greying)","bald","monobrow","chin dimple",
  "skin (pigmentation|color|darkness|tanning|aging|reflectance|autofluorescence)",
  "nose","ear","philtrum","outercanthal","nasal","eyebrow","eyelid","face shape","cranial width","cranial length",
  "consumption", "acid", "lipids", "content", "blood",
  
  # Sleep / circadian / lifestyle
  "chronotype","nap","napping","daytime sleepiness","snoring","morning person",
  "resting heart rate","heart rate variability","activity","step count","exercise","training","sports",
  
  # Exposures / environment / omics
  "microbiome|microbiota|bacterial|otu","beta diversity","urine","serum","plasma","cadmium|arsenic|mercur|nickel|pesticide",
  "blood pressure","biomarker","methylation","telomere","recombination","chromosome y", "stool",
  
  # Derived traits / composite
  "healthspan","morbidity[- ]free survival","longevity","aging","lifespan","parental longevity",
  
  # Procedures / non-disease operations
  "tonsillectomy","appendectomy","myringotomy","colonoscopy","biopsy","resection",
  
  # Other task or measurement constructs
  "task","performance","study participation","survey","questionnaire","symptom score",
  "subjective tiredness","self-reported tiredness","voice acoustics","gait","balance",
  "reaction time","psychometric","test","trait","index","function", "seropositivity",
  "intensity-contrast", "cell", "(?<![A-Za-z])[LR](?=\\s*($|[^A-Za-z]))", "deficiency", "fasting", "postprandial", "abundance",
  "attainment", "strain", "colour", "HLA", "monocyte", "granulocyte", "lymphocyte", "total", 
  "size", "math", "age", "medication", "dietary", "pattern", "ability", "density", "feeling", 
  "liking", "diameter", "length", "microbial", "intake", "net100", "dti ", "lead", 
  "inv-norm transformed", "shape", "synchrony", "amplitude", "ukb data field",
  "principal component", "frequency", "odor", "prevalence", "bone mass", "hla dr",
  "heteroplasmy", "spherical", "cost", "pollutants", "english", "meal", "modulation",
  "killer", "fvc", "political", "eyes", "perceived"
)


# -----------------------------
# 2) Classify with priority
# -----------------------------
# Priority rules:
#   1) DROP if clearly measurement-like (even if disease words present),
#   2) else KEEP if treatment/toxicity/survival in disease,
#   3) else KEEP if core disease,
#   4) else UNCERTAIN.

re_or <- function(x) paste(x, collapse = "|")
matches <- function(x, pattern_vec) grepl(re_or(pattern_vec), tolower(x), perl = TRUE)

normalize_txt <- function(x) {
  x <- gsub("[\\xA0\\x{2007}\\x{202F}]", " ", x, perl = TRUE)
  x <- gsub("[\\x{2010}-\\x{2015}\\x{2212}]", "-", x, perl = TRUE)            # hyphen-like → "-"
  x <- gsub("[ \t]+", " ", x, perl = TRUE)                    # collapse multiple spaces
  x <- gsub("\\s+([,;:.])", "\\1", x, perl = TRUE)            # remove space before punct
  x <- trimws(x)
  x
}

catalog_full <- catalog %>%
  transmute(
    trait      = `DISEASE/TRAIT`,
    trait_norm = normalize_txt(`DISEASE/TRAIT`)
  ) %>%
  filter(!is.na(trait_norm), trait_norm != "")

# 1b) Create one row per unique trait to assign labels (avoid 1-per-row flags)
trait_labels <- catalog_full %>%
  distinct(trait_norm) %>%
  mutate(
    flag_measure = matches(trait_norm, measurement_like),
    flag_treat   = matches(trait_norm, disease_treatment),
    flag_prog    = matches(trait_norm, disease_progression),
    flag_core    = matches(trait_norm, disease_core) | flag_prog,
    decision = dplyr::case_when(
      flag_measure ~ "DROP_MEASUREMENT",
      flag_treat   ~ "KEEP_DISEASE_TREATMENT",
      flag_core    ~ "KEEP_DISEASE",
      TRUE         ~ "UNCERTAIN"
    )
  ) %>%
  dplyr::select(trait_norm, decision)

# 1c) Join labels back to the FULL catalog rows, then keep disease-like rows
catalog_labeled <- catalog %>%
  left_join(trait_labels, by = c(`DISEASE/TRAIT` = "trait_norm"))

kept_rows <- catalog_labeled %>%
  mutate("trait_norm" = `DISEASE/TRAIT`) |>
#  filter(decision %in% c("KEEP_DISEASE","KEEP_DISEASE_TREATMENT")) %>%
  filter(decision %in% c("KEEP_DISEASE")) %>%
  mutate(reported = trait_norm)  # keep original display text

# 1d) Background counts per trait from the FULL catalog (no distinct())
trait_counts_bg <- kept_rows %>%
  dplyr::count(reported, name = "n")

total_catalog_counts <- sum(trait_counts_bg$n)

# Keep only traits that are in the disease-like background and attach catalog counts
df <- traits_long_count %>%
  inner_join(trait_counts_bg, by = "reported") %>%   # adds catalog n (popularity)
  mutate(reported = as.character(reported))

# total catalog count already computed above:
# total_catalog_counts <- sum(trait_counts_bg$n)

# (Optional) drop studies with zero/NA freq just in case
df <- df %>% filter(!is.na(freq), freq > 0)

# df <- df |>
#   dplyr::select(-"n.y") |>
#   dplyr::rename("n" = "n.x")

# --- Binomial enrichment per study ---
df_binom <- df %>%
  group_by(study) %>%
  mutate(
    n_study    = sum(freq, na.rm = TRUE),
    p_expected = n / total_catalog_counts,
    # cumulative binomial: P(X >= freq) with X ~ Binom(n_study, p_expected)
    p_binom    = 1 - pbinom(q = pmax(freq, 0) - 1, size = n_study, prob = p_expected),
    impressiveness = -log10(p_binom + 1e-300)  # robust floor for log
  ) %>%
  ungroup()

# Top 10 per study by impressiveness (ties included)
top_traits <- df_binom %>%
  group_by(study) %>%
  slice_max(order_by = impressiveness, n = 9, with_ties = TRUE) %>%
  ungroup()

# --- Plot (per study facets) ---
p <- ggplot(
  top_traits,
  aes(x = impressiveness,
      y = reorder(reported, impressiveness),
      size = freq)
) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ study, scales = "free_y") +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits per study (adjusted for GWAS catalog popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial.png"), width = 30, height = 8, dpi = 300)

# Reported traits less than 5 times in catalog out
top_traits_5 <- df_binom %>%
  filter(n > 5) %>%
  group_by(study) %>%
  slice_max(order_by = impressiveness, n = 10, with_ties = TRUE) %>%
  ungroup()

p_5 <- ggplot(
  top_traits_5,
  aes(x = impressiveness,
      y = reorder(reported, impressiveness),
      size = freq)
) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ study, scales = "free_y") +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits per study (adjusted for GWAS catalog popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial_cutoff.png"), width = 30, height = 8, dpi = 300)

# Reported traits less than 5 times in catalog out and more than twice in study
top_traits_5_2 <- df_binom %>%
  filter(n > 5,
         freq > 1) %>%
  group_by(study) %>%
  slice_max(order_by = impressiveness, n = 10, with_ties = TRUE) %>%
  ungroup()

p_5 <- ggplot(
  top_traits_5_2,
  aes(x = impressiveness,
      y = reorder(reported, impressiveness),
      size = freq)
) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ study, scales = "free_y") +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits per study (adjusted for GWAS catalog popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial_cutoff_freq.png"), width = 30, height = 8, dpi = 300)


# Now with mapped

# 1d) Background counts per trait from the FULL catalog (no distinct())
trait_counts_bg_mapped <- trait_counts_bg |>
  left_join(traits_aligned |>
              dplyr::select(c("reported", "mapped"))) |>
  filter(!grepl("neuroticism", mapped)) |>
  filter(!grepl("consumption", mapped)) |>
  filter(!grepl("distribution", mapped)) |>
  filter(!grepl("measurement", mapped)) |>
  filter(!grepl("amount", mapped)) |>
  dplyr::distinct() |>
  dplyr::select(-"reported") |>
  group_by(mapped) |>
  summarise(n = sum(n))

# Keep only traits that are in the disease-like background and attach catalog counts
df_m <- traits_long_count_mapped %>%
  inner_join(trait_counts_bg_mapped, by = "mapped") %>%   # adds catalog n (popularity)
  mutate(mapped = as.character(mapped))

# total catalog count already computed above:
total_catalog_counts_m <- sum(trait_counts_bg_mapped$n)

# (Optional) drop studies with zero/NA freq just in case
df_m <- df_m %>% filter(!is.na(freq), freq > 0)

# df <- df |>
#   dplyr::select(-"n.y") |>
#   dplyr::rename("n" = "n.x")

# --- Binomial enrichment per study ---
df_binom_m <- df_m %>%
  group_by(study) %>%
  mutate(
    n_study    = sum(freq, na.rm = TRUE),
    p_expected = n / total_catalog_counts,
    # cumulative binomial: P(X >= freq) with X ~ Binom(n_study, p_expected)
    p_binom    = 1 - pbinom(q = pmax(freq, 0) - 1, size = n_study, prob = p_expected),
    impressiveness = -log10(p_binom + 1e-300)  # robust floor for log
  ) %>%
  ungroup()

# Top 10 per study by impressiveness (ties included)
top_traits_m <- df_binom_m %>%
  group_by(study) %>%
  slice_max(order_by = impressiveness, n = 10, with_ties = TRUE) %>%
  ungroup()

# --- Plot (per study facets) ---
p <- ggplot(
  top_traits_m,
  aes(x = impressiveness,
      y = reorder(mapped, impressiveness),
      size = freq)
) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ study, scales = "free_y") +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits per study (adjusted for GWAS catalog popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial_mapped.png"), width = 30, height = 8, dpi = 300)

# Reported traits less than 5 times in catalog out
top_traits_5_m <- df_binom_m %>%
  filter(n > 5) %>%
  group_by(study) %>%
  slice_max(order_by = impressiveness, n = 10, with_ties = TRUE) %>%
  ungroup()

p_5 <- ggplot(
  top_traits_5_m,
  aes(x = impressiveness,
      y = reorder(mapped, impressiveness),
      size = freq)
) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ study, scales = "free_y") +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits per study (adjusted for GWAS catalog popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial_cutoff_mapped.png"), width = 30, height = 8, dpi = 300)

# Reported traits less than 5 times in catalog out and more than twice in study
top_traits_5_2_m <- df_binom_m %>%
  filter(n > 5,
         freq > 1) %>%
  group_by(study) %>%
  slice_max(order_by = impressiveness, n = 10, with_ties = TRUE) %>%
  ungroup()

p_5 <- ggplot(
  top_traits_5_2_m,
  aes(x = impressiveness,
      y = reorder(mapped, impressiveness),
      size = freq)
) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ study, scales = "free_y") +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits per study (adjusted for GWAS catalog popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial_cutoff_freq_mapped.png"), width = 30, height = 8, dpi = 300)

# top_traits |>
#   write_csv(here::here("Traits", "traits_all.csv"))
# top_traits_5 |>
#   write_csv(here::here("Traits", "traits_5.csv"))
# top_traits_5_2 |>
#   write_csv(here::here("Traits", "traits_5_2.csv"))
# top_traits_m |>
#   write_csv(here::here("Traits", "traits_all_m.csv"))
# top_traits_5_m |>
#   write_csv(here::here("Traits", "traits_all_5_m.csv"))
# top_traits_5_2_m |>
#   write_csv(here::here("Traits", "traits_all_5_2_m.csv"))



# PCC all
# HERE (K) then check counts for mapped also work, and check results for this one. Pick honestly whtever looks best LOL.
# Decide on method, correct all code, redo plots and tables.
# HERE (K) finally do everything also for full PCC. Put all plots in LaTeX and paste dummy paragraphs of what has been done, to revise and write better soon.
# Check if genes have any drugs, or if PCC has had any DR done in itself.
# Tallk to Anthony about (a) validation (or think), and (b) DR and methods (check vs DL or similar). Then, two days to have suggestions on both LR validation and DR, do in Oct LR and protocol, Nov all DR analyses.

# Clust PCC vs all
input <- read_table(paste0(dir_gwas_results,"/GWAS/", "PCCvsCovid/", "pccvscovid_FUMA.txt"))
genes <- read_table(paste0(dir_gwas_results,"/PCCvsCovid/","genes.txt"))
leads <- read_table(paste0(dir_gwas_results,"/PCCvsCovid/","leadSNPs.txt"))
snps <- read_table(paste0(dir_gwas_results,"/PCCvsCovid/","snps.txt"))

resultspccvscovid <- leads |>
  dplyr::left_join(input, by = c("rsID" = "RSID")) |>
  dplyr::mutate(OR = exp(BETA)) |>
  dplyr::select("chr",  "pos", "rsID", "P", "REF", "ALT", "OR", "SE") |>
  dplyr::full_join(genes |>
                     dplyr::select("rsID" = "IndSigSNPs", "gene" = "symbol", "type"),
                   by = "rsID") |>
  dplyr::mutate(cases = "PCC", controls = "COVID-19",
                lead = dplyr::if_else(is.na(chr), "No", "Yes")) |>
  dplyr::left_join(snps |>
                     dplyr::select(rsID, func))

genes_used <- resultspccvscovid %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)

pairs_by_gene_pcc <- genes_used %>%
  mutate(pairs = lapply(gene, genesymbol2gwas_pairs))

resultspccvscovid <- resultspccvscovid %>%
  dplyr::filter(!is.na(gene)) %>%
  left_join(
    pairs_by_gene_pcc %>%
      mutate(
        reported_traits = map_chr(
          pairs, ~ paste(na.omit(.x$reported), collapse = "&&& ")
        ),
        mapped_traits = map_chr(
          pairs, ~ paste(na.omit(.x$mapped), collapse = "&&& ")
        )
      ) %>%
      dplyr::select(gene, reported_traits, mapped_traits),
    by = "gene"
  )

# Long table (reported traits), normalised text
traits_aligned_pcc <- pairs_by_gene_pcc %>%
  tidyr::unnest(pairs) %>%
  filter(!is.na(reported), reported != "") %>%
  transmute(
    gene,
    reported
  ) %>%
  distinct(gene, reported)

traits_long_pcc <- resultspccvscovid |>
  dplyr::select(-c("reported_traits", "mapped_traits")) |>
  dplyr::left_join(traits_aligned |>
                     dplyr::select(c("gene", "reported", "mapped"))) |>
  dplyr::filter(!is.na(reported)) |>
  dplyr::filter(!is.na(mapped)) |>
  dplyr::distinct()

# Keep only traits that are in the disease-like background and attach catalog counts
df_pcc <- traits_long_pcc %>%
  group_by(reported) %>%
  tally(name = "freq") %>%
  inner_join(trait_counts_bg, by = "reported") %>%   # adds catalog n (popularity)
  mutate(reported = as.character(reported))

# --- Binomial enrichment per study ---
df_binom_pcc <- df_pcc %>%
  mutate(
    n_study    = sum(freq, na.rm = TRUE),
    p_expected = n / total_catalog_counts,
    # cumulative binomial: P(X >= freq) with X ~ Binom(n_study, p_expected)
    p_binom    = 1 - pbinom(q = pmax(freq, 0) - 1, size = n_study, prob = p_expected),
    impressiveness = -log10(p_binom + 1e-300)  # robust floor for log
  ) %>%
  ungroup()

# Top 10 per study by impressiveness (ties included)
top_traits_5_2_pcc <- df_binom_pcc %>%
  filter(n > 5,
         freq > 1) %>%
  slice_max(order_by = impressiveness, n = 10, with_ties = TRUE) %>%
  ungroup()

# --- Plot (per study facets) ---
p <- ggplot(
  top_traits_5_2_pcc,
  aes(x = impressiveness,
      y = reorder(reported, impressiveness),
      size = freq)  
) +
  geom_point(alpha = 0.9) +
  labs(
    x = "-log10(cumulative binomial p)",
    y = "Trait",
    title = "Top enriched disease traits in All PCC vs COVID-19",
    subtitle = "Color = GWAS catalog count (trait popularity)",
    size = "Count in study"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    
    plot.title = element_text(face = "bold")
  )

ggsave(here(dir_gwas_results,"Traits","top_traits_binomial_pcc.png"), width = 12, height = 8, dpi = 300)
