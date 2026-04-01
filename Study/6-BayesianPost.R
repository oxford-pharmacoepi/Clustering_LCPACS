############################################################
## Long COVID vs acute COVID genetic disentanglement
## and LC vs COVID-no-LC overlap (Bayesian-ish)
##
## - Your LC GWAS: .txt (tab-separated)
## - HGI B2/C2:    .tsv.gz (as downloaded)
############################################################

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

############################################################
## 0. USER PATHS — EDIT THESE
############################################################

# Your LC vs population summary stats (TXT, tab-separated)
lc_pop_path        <- here::here("Bayes/pcc_vsallfil.txt")

# Your LC vs COVID-no-LC summary stats (TXT, tab-separated)
lc_covid_ctrl_path <- here::here("Bayes/pcc_vscovid.txt")

# HGI Round 7 B2 (hospitalised vs population), TSV.GZ
hgi_b2_path <- here::here("Bayes/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz")

# (Optional) HGI Round 7 C2 (COVID vs population), TSV.GZ
hgi_c2_path <- here::here("Bayes/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz")

############################################################
## 1. Helper functions
############################################################

## Harmonise two GWAS on chr+pos and allele orientation
harmonise_chrpos <- function(df1, df2,
                             chr_col1 = "chr", pos_col1 = "pos",
                             chr_col2 = "chr", pos_col2 = "pos",
                             ea1 = "effect_allele", oa1 = "other_allele",
                             ea2 = "effect_allele", oa2 = "other_allele",
                             beta1 = "beta", se1 = "se",
                             beta2 = "beta", se2 = "se") {
  dt1 <- as.data.table(df1)
  dt2 <- as.data.table(df2)
  
  setnames(dt1, c(chr_col1, pos_col1, ea1, oa1, beta1, se1),
           c("chr_1", "pos_1", "ea_1", "oa_1", "beta_1", "se_1"))
  setnames(dt2, c(chr_col2, pos_col2, ea2, oa2, beta2, se2),
           c("chr_2", "pos_2", "ea_2", "oa_2", "beta_2", "se_2"))
  
  setkey(dt1, chr_1, pos_1)
  setkey(dt2, chr_2, pos_2)
  
  m <- dt1[dt2, nomatch = 0]
  
  m[, `:=`(chr = chr_1, pos = pos_1)]
  
  m[, `:=`(
    ea_1 = toupper(ea_1),
    oa_1 = toupper(oa_1),
    ea_2 = toupper(ea_2),
    oa_2 = toupper(oa_2)
  )]
  
  same <- m$ea_1 == m$ea_2 & m$oa_1 == m$oa_2
  flip <- m$ea_1 == m$oa_2 & m$oa_1 == m$ea_2
  
  keep <- same | flip
  m <- m[keep]
  flip <- flip[keep]
  
  # Flip beta_2 where alleles reversed
  m$beta_2[flip] <- -m$beta_2[flip]
  
  # Return clean tibble
  as_tibble(m) %>%
    dplyr::select(chr, pos,
           ea_1, oa_1, beta_1, se_1,
           ea_2, oa_2, beta_2, se_2,
           everything())
}

## Weighted regression through origin: y ~ alpha * x (IVW-style)
estimate_alpha <- function(beta_x, beta_y, se_y) {
  w <- 1 / (se_y^2)
  sum(w * beta_x * beta_y, na.rm = TRUE) / sum(w * beta_x^2, na.rm = TRUE)
}

## Delta-method SE for ratio R = y/x
ratio_se <- function(beta_x, beta_y, se_x, se_y) {
  sqrt( (se_y^2 / beta_x^2) + (beta_y^2 * se_x^2 / beta_x^4) )
}

## Wakefield Approximate Bayes Factor
wakefield_abf <- function(beta, se, W = 0.04) {
  V <- se^2
  z <- beta / se
  r <- W / (W + V)
  sqrt(1 - r) * exp(0.5 * r * z^2)
}


############################################################
## 2. Load your LC GWAS (.txt) and HGI (.tsv.gz)
############################################################

## 2.1 LC vs population
lc_pop_raw <- fread(lc_pop_path)

## If column names differ, adjust the rename() here
lc_pop <- lc_pop_raw %>%
  dplyr::rename(
    chr           = CHROM,
    pos           = POSITION,
    effect_allele = EFFECT_ALLELE,
    other_allele  = NON_EFFECT_ALLELE,
    beta          = BETA,
    se            = SE,
    p             = P
  )

## 2.2 LC vs COVID-no-LC
lc_covid_ctrl_raw <- fread(lc_covid_ctrl_path)

lc_covid_ctrl <- lc_covid_ctrl_raw %>%
  dplyr::rename(
    chr           = CHROM,
    pos           = POSITION,
    effect_allele = EFFECT_ALLELE,
    other_allele  = NON_EFFECT_ALLELE,
    beta          = BETA,
    se            = SE,
    p             = P
  )

## 2.3 HGI B2 (hospitalised vs population)
hgi_b2 <- fread(hgi_b2_path)

hgi_b2_small <- hgi_b2 %>%
  dplyr::transmute(
    chr           = `#CHR`,
    pos           = POS,
    rsid = rsid,
    effect_allele = ALT,  # HGI: beta is for ALT
    other_allele  = REF,
    beta          = all_inv_var_meta_beta,
    se            = all_inv_var_meta_sebeta,
    p             = all_inv_var_meta_p
  )

## 2.4 HGI C2 (COVID vs population) — optional
hgi_c2 <- fread(hgi_c2_path)
hgi_c2_small <- hgi_c2 %>%
  dplyr::transmute(
    chr           = `#CHR`,
    pos           = POS,
    rsid = rsid,
    effect_allele = ALT,
    other_allele  = REF,
    beta          = all_inv_var_meta_beta,
    se            = all_inv_var_meta_sebeta,
    p             = all_inv_var_meta_p
  )


############################################################
## 3. Lammi-style acute → LC slope and LC-specific loci
##    (HGI B2/C2 × LC vs population)
############################################################

## 3.1 Use only HGI GWS SNPs (fast & faithful to paper)
b2_gws <- hgi_b2_small %>%
  filter(p < 5e-8)

c2_gws <- hgi_c2_small %>%
  filter(p < 5e-8)

## 3.2 Harmonise LC vs population with HGI B2 GWS variants (by chr+pos)
lc_pop_b2 <- harmonise_chrpos(
  df1 = lc_pop,
  df2 = b2_gws,
  chr_col1 = "chr", pos_col1 = "pos",
  chr_col2 = "chr", pos_col2 = "pos",
  ea1 = "effect_allele", oa1 = "other_allele",
  ea2 = "effect_allele", oa2 = "other_allele",
  beta1 = "beta", se1 = "se",
  beta2 = "beta", se2 = "se"
)

# beta_1, se_1 = LC vs population
# beta_2, se_2 = HGI B2

## 3.3 Estimate alpha_B2 (global LC ~ alpha * severity)
alpha_b2 <- estimate_alpha(
  beta_x = lc_pop_b2$beta_2,
  beta_y = lc_pop_b2$beta_1,
  se_y   = lc_pop_b2$se_1
)

cat("Estimated alpha_B2 (LC ~ alpha * B2):", alpha_b2, "\n")

## 3.4 Do the same for C2, if available (infection susceptibility)
if (!is.null(c2_gws)) {
  lc_pop_c2 <- harmonise_chrpos(
    df1 = lc_pop,
    df2 = c2_gws,
    chr_col1 = "chr", pos_col1 = "pos",
    chr_col2 = "chr", pos_col2 = "pos",
    ea1 = "effect_allele", oa1 = "other_allele",
    ea2 = "effect_allele", oa2 = "other_allele",
    beta1 = "beta", se1 = "se",
    beta2 = "beta", se2 = "se"
  )
  
  alpha_c2 <- estimate_alpha(
    beta_x = lc_pop_c2$beta_2,
    beta_y = lc_pop_c2$beta_1,
    se_y   = lc_pop_c2$se_1
  )
  
  cat("Estimated alpha_C2 (LC ~ alpha * C2):", alpha_c2, "\n")
}

## 3.5 For each HGI B2 locus, test “LC-specificity” vs B2
b2_lammi <- lc_pop_b2 %>%
  mutate(
    beta_LC  = beta_1,
    se_LC    = se_1,
    beta_B2  = beta_2,
    se_B2    = se_2,
    ratio    = beta_LC / beta_B2,
    ratio_se = ratio_se(beta_B2, beta_LC, se_B2, se_LC),
    z_ratio  = (ratio - alpha_b2) / ratio_se,
    p_ratio  = 2 * pnorm(-abs(z_ratio)),
    beta_LC_pred = alpha_b2 * beta_B2,
    residual     = beta_LC - beta_LC_pred
  )

# HGI B2 loci whose LC effect is “too big” vs severity-mediated expectation
b2_lc_heavy <- b2_lammi %>%
  arrange(p_ratio)

cat("Top B2 loci with LC-heavy effects (smallest p_ratio):\n")
print(head(b2_lc_heavy, 10))

## 3.6 Quick diagnostic plot (B2 vs LC, with slope alpha_B2)
ggplot(b2_lammi, aes(x = beta_B2, y = beta_LC)) +
  geom_point() +
  geom_abline(intercept = 0, slope = alpha_b2, linetype = "dashed") +
  labs(
    x = "HGI B2 beta (hospitalised vs population)",
    y = "LC vs population beta",
    title = paste0("LC vs B2 relationship (alpha_B2 = ",
                   round(alpha_b2, 3), ")")
  )

## 3.4bis Variant-level “LC-heavy vs C2” check (Lammi-style for C2)

c2_lammi <- lc_pop_c2 %>%
  mutate(
    beta_LC  = beta_1,
    se_LC    = se_1,
    beta_C2  = beta_2,
    se_C2    = se_2,
    ratio    = beta_LC / beta_C2,
    ratio_se = ratio_se(beta_C2, beta_LC, se_C2, se_LC),
    z_ratio  = (ratio - alpha_c2) / ratio_se,
    p_ratio  = 2 * pnorm(-abs(z_ratio)),
    beta_LC_pred = alpha_c2 * beta_C2,
    residual     = beta_LC - beta_LC_pred
  )

# C2 loci with strongest deviation from the global C2->LC slope
c2_lc_heavy <- c2_lammi %>%
  arrange(p_ratio)

cat("Top C2 loci with strongest deviation from global alpha_C2 (smallest p_ratio):\n")
print(head(c2_lc_heavy, 10))

## Plot LC vs C2 with slope alpha_C2
ggplot(c2_lammi, aes(x = beta_C2, y = beta_LC)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = alpha_c2, linetype = "dashed") +
  labs(
    x = "HGI C2 beta (COVID vs population)",
    y = "LC vs population beta",
    title = paste0(
      "LC vs C2 relationship (alpha_C2 = ",
      round(alpha_c2, 3), ")"
    )
  )

## Add gene names for special loci (if you’ve annotated them)
annot_c2 <- c2_lammi %>%
  mutate(
    gene = case_when(
      rsid %in% c("rs73062389", "rs73062394") ~ "FOXP4",
      rsid %in% c("rs601338", "rs632111", "rs633372", "rs492602") ~ "FUT2",
      TRUE ~ NA_character_
    )
  )

gg_c2 <- annot_c2 %>%
  ggplot(aes(x = beta_C2, y = beta_LC)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(aes(alpha = 1 / (1 + p)), size = 2) +
  geom_abline(intercept = 0, slope = alpha_c2, linetype = "dashed") +
  geom_text_repel(
    data = subset(annot_c2, !is.na(gene)),
    aes(label = gene),
    size = 3.2
  ) +
  labs(
    x = "β (HGI C2: infection susceptibility)",
    y = "β (LC vs population)",
    title = paste0(
      "LC vs C2: variant-level enrichment (α = ",
      round(alpha_c2, 3), ")"
    )
  ) +
  theme_bw()

gg_c2





############################################################
## 4. LC vs population × LC vs COVID-no-LC
##    (shared LC effect vs control choice)
############################################################

## 4.1 (Performance-friendly) pre-filter:
##     Only keep SNPs with some signal in at least one GWAS
lc_pop_small <- lc_pop %>% filter(p < 1e-3)
lc_covid_small <- lc_covid_ctrl %>% filter(p < 1e-2)

## 4.2 Harmonise LC vs pop and LC vs COVID-no-LC by chr+pos
lcpc <- harmonise_chrpos(
  df1 = lc_pop_small,
  df2 = lc_covid_small,
  chr_col1 = "chr", pos_col1 = "pos",
  chr_col2 = "chr", pos_col2 = "pos",
  ea1 = "effect_allele", oa1 = "other_allele",
  ea2 = "effect_allele", oa2 = "other_allele",
  beta1 = "beta", se1 = "se",
  beta2 = "beta", se2 = "se"
)

lcpc <- lcpc %>%
  mutate(
    beta_LC_pop   = beta_1,
    se_LC_pop     = se_1,
    beta_LC_covid = beta_2,
    se_LC_covid   = se_2,
    p_LC_pop      = 2 * pnorm(-abs(beta_LC_pop   / se_LC_pop)),
    p_LC_covid    = 2 * pnorm(-abs(beta_LC_covid / se_LC_covid))
  )

## 4.3 Simple overlapping suggestive SNPs
overlap_suggestive <- lcpc %>%
  filter(
    p_LC_pop   < 1e-5,
    p_LC_covid < 1e-3,
    sign(beta_LC_pop) == sign(beta_LC_covid)
  ) %>%
  arrange(p_LC_pop)

cat("Number of overlapping suggestive SNPs (LCpop<1e-5 & LCcovid<1e-3, same direction):",
    nrow(overlap_suggestive), "\n")
print(head(overlap_suggestive, 10))

## 4.4 Bayesian “non-null in both” ranking (Wakefield ABF model)
W <- 0.04   # prior var on logOR ~ N(0, 0.04) => sd ~ 0.2
pi_00 <- 0.97
pi_10 <- 0.01
pi_01 <- 0.01
pi_11 <- 0.01

lcpc_abf <- lcpc %>%
  mutate(
    ABF_pop   = wakefield_abf(beta_LC_pop,   se_LC_pop,   W),
    ABF_covid = wakefield_abf(beta_LC_covid, se_LC_covid, W),
    BF_00 = 1,
    BF_10 = ABF_pop,
    BF_01 = ABF_covid,
    BF_11 = ABF_pop * ABF_covid
  ) %>%
  rowwise() %>%
  mutate(
    denom   = pi_00*BF_00 + pi_10*BF_10 + pi_01*BF_01 + pi_11*BF_11,
    post_00 = (pi_00*BF_00) / denom,
    post_10 = (pi_10*BF_10) / denom,
    post_01 = (pi_01*BF_01) / denom,
    post_11 = (pi_11*BF_11) / denom
  ) %>%
  ungroup()

shared_lc_snps <- lcpc_abf %>%
  arrange(desc(post_11))

cat("Top SNPs with high posterior of being non-null in BOTH LC GWAS (post_11):\n")
print(head(shared_lc_snps, 10))

# If you want a “hard” cut, e.g. post_11 > 0.5:
shared_strong <- shared_lc_snps %>% filter(post_11 > 0.5)
cat("Number of SNPs with post_11 > 0.5:", nrow(shared_strong), "\n")

# Nicer tables and plots for this part

## Categorise by max posterior class
lcpc_abf_summary <- lcpc_abf %>%
  mutate(
    category = case_when(
      post_11 == pmax(post_00, post_10, post_01, post_11) ~ "Non-null in both (11)",
      post_10 == pmax(post_00, post_10, post_01, post_11) ~ "LCpop only (10)",
      post_01 == pmax(post_00, post_10, post_01, post_11) ~ "LCcovid only (01)",
      TRUE                                                 ~ "Null (00)"
    )
  )

## Table of top “shared” SNPs
shared_table <- lcpc_abf_summary %>%
  arrange(desc(post_11)) %>%
  transmute(
    chr,
    pos,
    rsid = ifelse("rsid" %in% names(.), rsid, NA_character_),
    effect_allele = ea_1,
    other_allele  = oa_1,
    beta_LC_pop,
    se_LC_pop,
    p_LC_pop,
    beta_LC_covid,
    se_LC_covid,
    p_LC_covid,
    ABF_pop,
    ABF_covid,
    post_00, post_10, post_01, post_11,
    category
  )

head(shared_table, 20)

#readr::write_tsv(shared_table, "LC_vs_COVID_shared_ABF.tsv")

shared_strict <- shared_table %>% filter(post_11 > 0.999)

library(ggrepel)

lcpc_abf_summary <- lcpc_abf %>%
  mutate(
    category = case_when(
      post_11 == pmax(post_00, post_10, post_01, post_11) ~ "Non-null in both (11)",
      post_10 == pmax(post_00, post_10, post_01, post_11) ~ "LCpop only (10)",
      post_01 == pmax(post_00, post_10, post_01, post_11) ~ "LCcovid only (01)",
      TRUE                                                 ~ "Null (00)"
    )
  )

gg_lc_lc <- ggplot(lcpc_abf_summary,
                   aes(x = beta_LC_pop, y = beta_LC_covid)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(aes(colour = category,
                 size = -log10(p_LC_pop * p_LC_covid + 1e-300)),
             alpha = 0.7) +
  scale_size_continuous(name = "-log10(p1 * p2)", range = c(1, 5)) +
  scale_colour_manual(
    values = c(
      "Non-null in both (11)" = "firebrick",
      "LCpop only (10)"       = "steelblue",
      "LCcovid only (01)"     = "darkorange",
      "Null (00)"             = "grey70"
    )
  ) +
  geom_abline(intercept = 0, slope = coef(fit_lc_lc)[["beta_1"]],
              linetype = "dashed") +
  labs(
    x = "β (LC vs population)",
    y = "β (LC vs COVID-no-LC)",
    colour = "Posterior category",
    title = "LC GWAS with different controls: strong effect-size concordance"
  ) +
  theme_bw()

gg_lc_lc


## 1. Correlations & sign agreement
cor_beta <- cor(lcpc$beta_1, lcpc$beta_2, use = "complete.obs")
sign_agree <- mean(sign(lcpc$beta_1) == sign(lcpc$beta_2), na.rm = TRUE)

cat("Correlation of betas (LC vs pop, LC vs COVID):", cor_beta, "\n")
cat("Proportion with same sign:", sign_agree, "\n")

fit_lc_lc <- lm(beta_2 ~ 0 + beta_1, data = lcpc)  # through origin
summary(fit_lc_lc)

shared_strict <- lcpc_abf_summary %>%
  filter(post_11 > 0.999) %>%
  arrange(desc(post_11))

head(shared_strict, 20)

############################################################
## End of script
############################################################
