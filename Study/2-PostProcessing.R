gwas <- read.delim(here("Results/LongCovid_assoc.regenie.merged.txt"))
gwas <- gwas |> as_tibble()

# Separate chromosomes x and xy
gwas23 <- gwas |> filter(CHROM == 23) 
x <- which(gwas23$GENPOS[1:length(gwas23$GENPOS)-1] > gwas23$GENPOS[2:length(gwas23$GENPOS)])

gwas |>
  filter(CHROM != 23) |>
  rbind(
    gwas23 |> mutate(CHROM = if_else(row_number() > x, 24, CHROM))
  ) |>
  write_delim(here("Results/LongCovid.txt"), delim = "\t")
