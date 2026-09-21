library(readxl)
library(data.table)

x <- read_excel('Draft/DPA_Replication_VL_20260120.xlsx', sheet = 'DPA_SNPs_FinnGen')
x <- x[!is.na(x$rsID) & !is.na(x$LOG10P), ]

x_dt <- as.data.table(x)
setorder(x_dt, rsID, -LOG10P)
best <- x_dt[, .SD[1], by = rsID]

map_file <- function(a) {
  if (a == 'AllPCCvsGenPop') return('GWAS/pcc_vsallfil.txt.gz')
  if (a == 'Subtype1vsPopCtrl') return('GWAS/clust1_vsallfil.txt.gz')
  if (a == 'Subtype2vsPopCtrl') return('GWAS/clust2_vsallfil.txt.gz')
  if (a == 'Subtype3vsCOVID') return('GWAS/clust3_vsno.txt.gz')
  NA_character_
}

get_beta <- function(file, snp) {
  if (is.na(file) || !file.exists(file)) return(NA_real_)
  d <- fread(file, select = c('SNP', 'BETA'))
  h <- d[SNP == snp, BETA]
  if (length(h) == 0) return(NA_real_)
  as.numeric(h[1])
}

best[, file := vapply(DPA_Analysis, map_file, character(1))]
best[, my_beta := mapply(get_beta, file, rsID)]
best[, direction_same := ifelse(is.na(my_beta) | is.na(BETA), NA, sign(my_beta) == sign(BETA))]
best[, fingen_p := 10^(-LOG10P)]
best[, replication_tier := fifelse(LOG10P >= 3, 'p<=0.001',
                            fifelse(LOG10P >= 2, 'p<=0.01',
                            fifelse(LOG10P >= 1.30103, 'p<0.05', 'not_nominal')))]

out <- best[, .(
  SNP = rsID,
  analysis = DPA_Analysis,
  gene = NearestGene,
  fingen_endpoint = FinnGen_Analysis,
  fingen_beta = BETA,
  fingen_se = SE,
  fingen_log10p = LOG10P,
  fingen_p = fingen_p,
  my_beta = my_beta,
  direction_same = direction_same,
  replication_tier = replication_tier
)]
setorder(out, -fingen_log10p)

fwrite(out, 'Draft/fingen_replication_summary.csv')

cat('n_index_snps=', nrow(out), '\n', sep = '')
cat('nominal_replications=', sum(out$fingen_log10p >= 1.30103, na.rm = TRUE), '\n', sep = '')
cat('p_le_0.01=', sum(out$fingen_log10p >= 2, na.rm = TRUE), '\n', sep = '')
cat('p_le_0.001=', sum(out$fingen_log10p >= 3, na.rm = TRUE), '\n', sep = '')
cat('direction_concordant_n=', sum(out$direction_same %in% TRUE, na.rm = TRUE), '\n', sep = '')
cat('direction_discordant_n=', sum(out$direction_same %in% FALSE, na.rm = TRUE), '\n', sep = '')
cat('direction_missing_n=', sum(is.na(out$direction_same)), '\n', sep = '')
print(out)
