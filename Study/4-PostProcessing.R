# Quality control function
quality_control <- function(gwas) {
  gwas |> 
    dplyr::filter(!is.na(LOG10P)) |> 
    dplyr::filter(nchar(NON_EFFECT_ALLELE) == 1, nchar(EFFECT_ALLELE) == 1) |>
    dplyr::filter(INFO >= 0.8) |>
    dplyr::mutate(dif = 1 - A1FREQ) |>
    dplyr::filter(!(A1FREQ < 0.01 | dif < 0.01)) |>
    dplyr::select(-dif)
}

# QQ plot function
plot_qq <- function(gwas_data, output_name, x_lim = 6.5, y_lim = 8,
                    color = "#2166AC",
                    dot_size = 0.1,
                    reduce_dataset = 0.01) {
  n <- nrow(gwas_data)
  ci <- .95
  
  dat <- data.frame(
    observed = sort(gwas_data$LOG10P),
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
  
  ggsave(plot = plot2, width = 30, height = 15, dpi = 300, units = "cm",
         filename = paste0(dir_gwas_results,"/GWAS/", output_name))
}

# FUMA format function (get rsId)
# only need to do the next part (generate rds of rsids) once
# get all of the gwas results files for the chrom/positions available
gwas <- read_table(paste0(dir_gwas_results,"/GWAS/clust1_vs23.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas2 <- read_table(paste0(dir_gwas_results,"/GWAS/clust2_vs13.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas3 <- read_table(paste0(dir_gwas_results,"/GWAS/clust3_vs12.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas4 <- read_table(paste0(dir_gwas_results,"/GWAS/clust1_vsno.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas5 <- read_table(paste0(dir_gwas_results,"/GWAS/clust2_vsno.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas6 <- read_table(paste0(dir_gwas_results,"/GWAS/clust3_vsno.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas7 <- read_table(paste0(dir_gwas_results,"/GWAS/clust1_vsall.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas8 <- read_table(paste0(dir_gwas_results,"/GWAS/clust2_vsall.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()
gwas9 <- read_table(paste0(dir_gwas_results,"/GWAS/clust3_vsall.txt")) |>
  dplyr::filter(P < 0.05) |>
  dplyr::select(c("CHROM", "POSITION", "NON_EFFECT_ALLELE", "NON_EFFECT_ALLELE")) |>
  dplyr::distinct()

gwas_all <- gwas |>
  dplyr::union_all(gwas2) |>
  dplyr::union_all(gwas3) |>
  dplyr::union_all(gwas4) |>
  dplyr::union_all(gwas5) |>
  dplyr::union_all(gwas6) |>
  dplyr::union_all(gwas7) |>
  dplyr::union_all(gwas8) |>
  dplyr::union_all(gwas9) |>
  dplyr::distinct()

# path to the VCF
vcf_file <- "GCF_000001405.40.gz"

# define your variants as a GRanges object
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

gwas_all$CHROM <- as.character(gwas_all$CHROM)
gwas_all$CHR_REFSEQ <- refseq_map[as.character(gwas_all$CHROM)]

# tibble with all chrom, positions and rsids of interest
results <- gwas_all %>%
  group_split(CHR_REFSEQ) %>%
  lapply(function(subset_df) {
    message("Querying ", unique(subset_df$CHR_REFSEQ), " with ", nrow(subset_df), " positions...")

    # Build one GRanges per chromosome
    gr <- GRanges(subset_df$CHR_REFSEQ, IRanges(subset_df$POSITION, subset_df$POSITION))

    # One call per chromosome
    param <- ScanVcfParam(which = gr, info = "RS")
    vcf <- readVcf(vcf_file, genome = "hg38", param = param)
    rr <- rowRanges(vcf)
    ref <- as.character(mcols(rr)$REF)
    alt_list <- lapply(mcols(rr)$ALT, as.character)

    # Extract RSID + positions
    tibble(
      CHR_REFSEQ = as.character(seqnames(rr)),
      POSITION   = start(rr),
      RSID       = names(rr),
      REF        = ref,
      ALT        = alt_list
    )  %>%
      tidyr::unnest_longer(ALT)
  }) %>%
  bind_rows()

saveRDS(results, paste0(dir_gwas_results,"/GWAS/","all_rsids.rds"), compress = "gzip")

# Read rsIds tibble
allRsIds <- readRDS(paste0(dir_gwas_results,"/GWAS/","all_rsids.rds"))

# Function to get fuma format (with rsId)
get_fuma_format <- function(gwas_data, output_file) {
  gwas_data$CHROM <- as.character(gwas_data$CHROM)
  gwas_data$CHR_REFSEQ <- refseq_map[as.character(gwas_data$CHROM)]
  
  gwas_data_with_rsid <- gwas_data |>
    dplyr::filter(P < 0.05) |>
    left_join(allRsIds, by = c("CHR_REFSEQ", "POSITION")) %>%
    filter(
      (toupper(NON_EFFECT_ALLELE) == REF & toupper(EFFECT_ALLELE) == ALT) |
        (toupper(NON_EFFECT_ALLELE) == ALT & toupper(EFFECT_ALLELE) == REF)
    ) %>%
    distinct()
  
  gwas_data_with_rsid  |>
    dplyr::rename(
      "CHR" = "CHROM",
      "BP" = "POSITION",
      "A1" = "EFFECT_ALLELE",
      "A2" = "NON_EFFECT_ALLELE",
    ) |>
    dplyr::select(-c("CHR", "BP", "SNP")) |>
    dplyr::filter(!is.na(RSID)) |>
    write.table(paste0(dir_gwas_results,"/GWAS/",output_file), row.names = FALSE, quote = FALSE)
}

# Dataset list, apply to all
datasets <- tibble::tribble(
  ~file,                                 ~manhattan_png,     ~qq_png,        ~fuma_file,       ~filter_dragen,  ~pval_file,   ~outcome_folder, 
  paste0(dir_gwas_results,"/GWAS/clust1_vs23.txt"), "clust1vs23.png", "clust1vs23_QQ.png", "clust1vs23_FUMA.txt", TRUE,    "clust1vs23_pval.txt",   "Clust1vs23",
  paste0(dir_gwas_results,"/GWAS/clust2_vs13.txt"), "clust2vs13.png", "clust2vs13_QQ.png", "clust2vs13_FUMA.txt", FALSE,   "clust2vs13_pval.txt",   "Clust2vs13",   
  paste0(dir_gwas_results,"/GWAS/clust3_vs12.txt"), "clust3vs12.png", "clust3vs12_QQ.png", "clust3vs12_FUMA.txt", FALSE,   "clust3vs12_pval.txt",   "Clust3vs12",
  paste0(dir_gwas_results,"/GWAS/clust1_vsno.txt"), "clust1vsno.png", "clust1vsno_QQ.png", "clust1vsno_FUMA.txt", FALSE,   "clust1vsno_pval.txt",   "Clust1vsno",
  paste0(dir_gwas_results,"/GWAS/clust2_vsno.txt"), "clust2vsno.png", "clust2vsno_QQ.png", "clust2vsno_FUMA.txt", TRUE,   "clust2vsno_pval.txt",   "Clust2vsno",
  paste0(dir_gwas_results,"/GWAS/clust3_vsno.txt"), "clust3vsno.png", "clust3vsno_QQ.png", "clust3vsno_FUMA.txt", FALSE,   "clust3vsno_pval.txt",   "Clust3vsno",
  paste0(dir_gwas_results,"/GWAS/clust1_vsall.txt"), "clust1vsall.png", "clust1vsall_QQ.png", "clust1vsall_FUMA.txt", FALSE,   "clust1vsall_pval.txt",   "Clust1vsall",
  paste0(dir_gwas_results,"/GWAS/clust2_vsall.txt"), "clust2vsall.png", "clust2vsall_QQ.png", "clust2vsall_FUMA.txt", TRUE,   "clust2vsall_pval.txt",   "Clust2vsall",
  paste0(dir_gwas_results,"/GWAS/clust3_vsall.txt"), "clust3vsall.png", "clust3vsall_QQ.png", "clust3vsall_FUMA.txt", FALSE,   "clust3vsall_pval.txt",   "Clust3vsall"
)

for (i in seq_len(nrow(datasets))) {
  params <- datasets[i, ]
  
  gwas <- read_table(params$file)
  
  # Manhattan
  if (params$filter_dragen) {
    gwas <- gwas |> filter(grepl("DRAGEN", SNP))
  }
  
  gwas <- quality_control(gwas)
  
  folder <- params$outcome_folder
  dir.create(paste0(dir_gwas_results,"/GWAS/",folder))
  
  # Reduce dataset for plotting speed
  gwas_reduced <- bind_rows(
    gwas |> filter(LOG10P > 3),
    gwas |> filter(LOG10P <= 3, LOG10P > 2) |> sample_frac(0.5),
    gwas |> filter(LOG10P <= 2, LOG10P > 1) |> sample_frac(0.1),
    gwas |> filter(LOG10P <= 1) |> sample_frac(0.01)
  )
  
  # Calculate cumulative positions
  don <- gwas_reduced |>
    group_by(CHROM) |>
    summarise(chr_len = max(POSITION), .groups = "drop") %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(gwas_reduced, ., by = "CHROM") %>%
    arrange(CHROM, POSITION) |>
    mutate(BPcum = POSITION + tot)
  
  axisdf <- don %>% group_by(CHROM) %>% 
    summarise(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")
  
  # Plot
  colors <- c("#92C5DE", "#4393C3", "#2166AC")
  y_lim <- 9
  
  p <- ggplot(don, aes(x = BPcum, y = LOG10P)) +
    geom_point(aes(color = as.factor(CHROM)), size = 1) +
    scale_color_manual(values = rep(colors, 22)) +
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center, expand = c(0.015, 0.015)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, y_lim, 2)) +
    coord_cartesian(ylim = c(0, y_lim)) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.text.x = element_text(margin = margin(t = 1), size = 9),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.3) +
    geom_hline(yintercept = -log10(5e-6), color = "grey", linewidth = 0.3) +
    labs(x = 'Chromosome', y = expression(-log[10](P)))
  
  ggsave(plot = p, width = 30, height = 15, dpi = 300, units = "cm",
         filename = here("GWAS", params$outcome_folder, params$manhattan_png))
  # ggsave(plot = p, width = 30, height = 15, units = "cm", filename = here("GWAS", output_name), device = "pdf") # PDF version
  
  pvals <- list(
    suggestive = don |> filter(LOG10P > -log10(5e-06)),
    genomewide = don |> filter(LOG10P > -log10(5e-08))
  )

  write.csv(pvals$suggestive, paste0(dir_gwas_results,"/GWAS/", params$pval_file))
  
  # QQ
  plot_qq(gwas, paste0(params$outcome_folder,"/", params$qq_png))
  
  # FUMA
  get_fuma_format(gwas, paste0(params$outcome_folder,"/", params$fuma_file))

  print(paste0("Done ", params$file))
}
