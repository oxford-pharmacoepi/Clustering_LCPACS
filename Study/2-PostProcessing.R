gwas <- read.table(here("Results/LongCovid_assoc.regenie.merged.txt"), skip = 1)
gwas <- gwas |> as_tibble()
colnames(gwas) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                    "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

gwas_old <- read.delim(here("Results/LongCovid_assoc.regenie.merged_20241029.txt"))
gwas_old <- gwas_old |> as_tibble()

gwas_wide <- read.delim(here("Results/LongCovid_wide_assoc.regenie.merged.txt"), skip = 1)
gwas_wide <- gwas_wide |> as_tibble()
colnames(gwas_wide) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                    "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

gwas_pacs <- read.delim(here("Results/Pacs_assoc.regenie.merged.txt"), skip = 1)
gwas_pacs <- gwas_pacs |> as_tibble()
colnames(gwas_pacs) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                    "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

gwas_pacs_a <- read.delim(here("Results/Pacs_arterial_assoc.regenie.merged.txt"), skip = 1)
gwas_pacs_a <- gwas_pacs_a |> as_tibble()
colnames(gwas_pacs_a) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                         "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

gwas_pacs_v <- read.delim(here("Results/Pacs_venous_assoc.regenie.merged.txt"), skip = 1)
gwas_pacs_v <- gwas_pacs_v |> as_tibble()
colnames(gwas_pacs_v) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                         "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")

plot_helper <- function(gwas, name) {
  # Separate chromosomes x and xy
  gwas23 <- gwas |> filter(CHROM == 23) 
  x <- which(gwas23$GENPOS[1:length(gwas23$GENPOS)-1] > gwas23$GENPOS[2:length(gwas23$GENPOS)])
  
  gwas <- gwas |>
    filter(CHROM != 23) |>
    rbind(
      gwas23 |> mutate(CHROM = if_else(row_number() > x, 24, .data$CHROM))
    ) 
  
  gwas |>
    write_delim(here(paste0("Results/",name,".txt")), delim = "\t")
  
  count_snps <- tibble(
    step = "initial",
    counts = gwas |> tally() |> pull()
  )
  
  gwas <- gwas |>
    dplyr::mutate(P = 10^(-LOG10P)) %>%
    filter(!is.na(LOG10P)) 
  
  count_snps <- union_all(
    count_snps,
    tibble(
      step = "LOG10P not NA",
      counts = gwas |> tally() |> pull()
    ))
  
  if(("INFO" %in% colnames(gwas))) {
    gwas <- gwas |>
      filter(INFO > 0.8)   
    
    count_snps <- union_all(
      count_snps,
      tibble(
        step = "INFO greater than 0.8",
        counts = gwas |> tally() |> pull()
      ))
    
    gwas <- gwas |>
      mutate(dif = 1-A1FREQ) |>
      filter(!(A1FREQ < 0.01 | dif < 0.01)) |>
      select(-"dif")   
    
    count_snps <- union_all(
      count_snps,
      tibble(
        step = "MAF greater than 0.01",
        counts = gwas |> tally() |> pull()
      ))
  } 
  
  gwas |>
    filter(LOG10P > 7) |>
    print()

    count_snps |>
      write_csv(here(paste0("Results/",name,"_post_processing_filtering.txt")))
  
  plot1 <- getManhattanPlot(gwas, 8) #ylim?
  plot2 <- getQQPlot(gwas, x_lim = 6.5, y_lim = 8)
  
  plot1 + plot2
  ggsave(here(paste0("Results/",name,".png")), width = 16)

  }

plot_helper(gwas, "LongCovid") # SNP rs200784681 (CHROM 12, GENPOS 129417693)
plot_helper(gwas_wide, "LongCovid_wide") # SNP rs6453106 (CHROM 5, GENPOS 74298946)
plot_helper(gwas_old, "LongCovid_old") # SNP rs200784681 (CHROM 12, GENPOS 129417693)
plot_helper(gwas_pacs, "Pacs") # no SNP 
plot_helper(gwas_pacs_a, "Pacs_arterial") # no SNP 
plot_helper(gwas_pacs_v, "Pacs_venous") # no SNP 

# no difference 0.7 SNPs ??

