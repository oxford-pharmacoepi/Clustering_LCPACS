#!/usr/bin/env Rscript

message("Reading HapMap3 SNP list")
hm3 <- read.table("Draft/ldsc/eur_w_ld_chr/w_hm3.snplist", header = TRUE, stringsAsFactors = FALSE)
hm3_ids <- unique(hm3$SNP)

message("Reading full rsID map from GWAS/all_rsids.rds")
rsmap <- readRDS("GWAS/all_rsids.rds")

message("Filtering to biallelic HapMap3 SNPs")
rsmap <- rsmap[rsmap$RSID %in% hm3_ids, c("CHR_REFSEQ", "POSITION", "RSID", "REF", "ALT")]
rsmap <- rsmap[nchar(rsmap$REF) == 1 & nchar(rsmap$ALT) == 1, ]
rsmap <- rsmap[rsmap$REF %in% c("A", "C", "G", "T") & rsmap$ALT %in% c("A", "C", "G", "T"), ]
rsmap <- unique(rsmap)

message("Writing Draft/ldsc/hm3_rsid_position_allele_map_full.tsv")
write.table(
  rsmap,
  file = "Draft/ldsc/hm3_rsid_position_allele_map_full.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Rows written: ", nrow(rsmap))
