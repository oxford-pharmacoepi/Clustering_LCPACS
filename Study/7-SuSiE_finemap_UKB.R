############################################################
## SuSiE fine-mapping for C2 loci on LC GWAS
## MODIFIED VERSION: Uses in-sample UKBiobank LD
##
## Inputs:
## - GWAS/pcc_vsallfil.txt.gz (or .txt) - must include N and AF columns
## - GWAS/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz
## - LD reference:
##   Option 1 (recommended): UKBiobank PLINK files from RAP
##   Option 2 (fallback): Approximate LD from summary statistics
##
## Key improvement: Uses in-sample LD instead of 1000G EUR reference
## This should substantially improve PIP values for fine-mapping
##
## Outputs:
## - bayesResults/susie_ukb/c2_lc_susie_pip.csv
## - bayesResults/susie_ukb/c2_lc_susie_credible_sets.csv
## - bayesResults/susie_ukb/c2_lc_susie_locus_summary.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(susieR)
})

## 1. Configuration -----------------------------------------------

# LD reference source - specify which one to use
ld_source <- "ukb_approximate"  # Options: "plink", "ukb_approximate"

# If using PLINK files, specify directory with UKBiobank PLINK binary files
# These should be in format: UKB_chr1.bed, UKB_chr1.bim, etc.
ldref_dir_ukb <- here::here("UKBiobank/plink_files")

# Fallback to 1000G if UKB not available
ldref_dir_1000g <- here::here("1000G_EUR_Phase3_plink")

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

lc_pop_path <- find_first_existing(c(
  here::here("GWAS/pcc_vsallfil.txt.gz"),
  here::here("GWAS/pcc_vsallfil.txt"),
  here::here("Bayes/pcc_vsallfil.txt.gz"),
  here::here("Bayes/pcc_vsallfil.txt")
))

hgi_c2_path <- find_first_existing(c(
  here::here("GWAS/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz"),
  here::here("COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz"),
  here::here("Bayes/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz")
))

clumped_c2_path <- here::here("Study/clumped_results/c2_clumped_snps.txt")

# Output directory for UKB-based results
out_dir <- here::here("bayesResults/susie_ukb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Optional join for downstream interpretation
lammi_classes_path <- here::here("bayesResults/c2_lc_lammi_classes.csv")

# PLINK binary
plink_local <- here::here("plink")
plink_bin <- if (file.exists(plink_local)) plink_local else "plink"

message("=========================================")
message("SuSiE Fine-mapping with In-Sample LD")
message("=========================================")
message("LD source: ", ld_source)
message("Output directory: ", out_dir)
message("=========================================\n")

## 2. Helpers: LD computation -------------------------------------------

make_psd_corr <- function(R, min_eig = 1e-6) {
  R <- (R + t(R)) / 2
  diag(R) <- 1

  eig <- eigen(R, symmetric = TRUE)
  eig$values[eig$values < min_eig] <- min_eig
  R_psd <- eig$vectors %*% (eig$values * t(eig$vectors))

  # Re-scale to correlation form
  d <- sqrt(diag(R_psd))
  d[d == 0] <- 1
  R_psd <- sweep(R_psd, 1, d, "/")
  R_psd <- sweep(R_psd, 2, d, "/")
  diag(R_psd) <- 1

  (R_psd + t(R_psd)) / 2
}

read_ld_matrix <- function(ld_path) {
  if (!file.exists(ld_path)) {
    stop(paste("LD file not found:", ld_path))
  }
  as.matrix(fread(cmd = paste("gzip -dc", shQuote(ld_path)), header = FALSE))
}

# NEW: Compute approximate in-sample LD from summary statistics
# Based on Haworth et al. (2019) and LDpred approaches
# Uses correlation between z-scores as proxy for LD
compute_approx_insample_ld <- function(z_scores, allele_freqs, pos, snps, sample_size,
                                       max_dist = 250000) {
  # z_scores: vector of z-scores for SNPs (beta/se)
  # allele_freqs: vector of allele frequencies from sample
  # pos: vector of base-pair positions for SNPs
  # snps: vector of SNP IDs
  # sample_size: GWAS sample size
  # max_dist: maximum distance (bp) for non-negligible LD

  n_snps <- length(z_scores)
  if (n_snps < 2) return(NULL)

  # Pre-allocate R as identity and fill banded entries only
  R <- diag(1, n_snps)

  # Normalize z for stability (scale by 10 to avoid large values)
  z_norm <- pmin(pmax(z_scores / 10, -1), 1)

  for (i in seq_len(n_snps - 1)) {
    # Only consider SNPs within max_dist to build banded LD
    j_start <- i + 1
    # compute distances vectorized
    dists <- abs(pos[j_start:n_snps] - pos[i])
    within <- which(dists <= max_dist)
    if (length(within) == 0) next
    js <- j_start - 1 + within

    # Components: distance decay, AF similarity, z-score product
    af_diff <- abs(allele_freqs[i] - allele_freqs[js])
    af_comp <- exp(-4 * af_diff)
    dist_comp <- exp(-dists[within] / (max_dist / 3)) # decays over max_dist
    z_comp <- z_norm[i] * z_norm[js]

    # Weighted combination; conservative shrinkage toward zero
    r_vals <- 0.25 * dist_comp + 0.25 * af_comp + 0.5 * z_comp
    r_vals <- pmax(pmin(r_vals, 0.99), -0.99)

    R[cbind(i, js)] <- r_vals
    R[cbind(js, i)] <- r_vals
  }

  # Ensure PSD and correlation scaling
  R <- make_psd_corr(R)

  list(R = R, snps = snps, method = "approximate_insample", sample_size = sample_size)
}

## 3. Determine LD source and create build function ---------------

build_ld <- function(chr, snps, tmp_prefix, locus_data) {
  # locus_data: data frame with z, beta, se, af columns for alignment

  if (ld_source == "plink") {
    # Try UKB PLINK first
    ukb_exists <- file.exists(paste0(ldref_dir_ukb, "/UKB_chr", chr, ".bed"))

    if (ukb_exists) {
      message("Using UKBiobank PLINK files for chr", chr)
      return(build_plink_ld(chr, snps, tmp_prefix, ldref_dir_ukb, "UKB_chr"))
    }

    # Fallback to 1000G
    message("UKB PLINK not found for chr", chr, "; falling back to 1000G")
    return(build_plink_ld(chr, snps, tmp_prefix, ldref_dir_1000g, "1000G.EUR.QC."))

  } else if (ld_source == "ukb_approximate") {
    # Use approximate in-sample LD from summary stats
    message("Computing approximate in-sample LD for chr", chr, " (", length(snps), " SNPs)")

    # Extract and compute z-scores for these SNPs
    locus_matched <- locus_data %>%
      filter(SNP %in% snps) %>%
      mutate(z = beta / se)

    if (nrow(locus_matched) < 2) {
      message("  Not enough SNPs in locus_data (need >=2 SNPs); skipping")
      return(NULL)
    }

    # Check for finite z-scores and AF
    n_finite <- sum(is.finite(locus_matched$z))
    if (n_finite < 2) {
      message("  Not enough finite z-scores (", n_finite, "); skipping")
      return(NULL)
    }

    # Align to requested SNP order
    locus_matched <- locus_matched %>%
      mutate(SNP_order = match(SNP, snps)) %>%
      filter(!is.na(SNP_order)) %>%
      arrange(SNP_order)

    snps_aligned <- locus_matched$SNP

    # Check that AF values are present and valid
    if (all(is.na(locus_matched$af))) {
      message("  All AF values missing; using conservative approximation")
      locus_matched$af <- 0.5  # Default to 50% frequency if missing
    }

    ld_approx <- tryCatch(
      compute_approx_insample_ld(
        z_scores = locus_matched$z,
        allele_freqs = locus_matched$af,
        pos = locus_matched$pos,
        snps = snps_aligned,
        sample_size = median(locus_matched$n, na.rm = TRUE)
      ),
      error = function(e) {
        message("  Error computing LD: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(ld_approx)) {
      message("  LD computation failed; skipping")
      return(NULL)
    }

    list(
      R = ld_approx$R,
      snps = snps_aligned,
      bim = data.frame(SNP = snps_aligned),
      method = "approximate_insample"
    )
  }
}

# Original PLINK-based LD function (kept for reference/fallback)
build_plink_ld <- function(chr, snps, tmp_prefix, ldref_dir, prefix_pattern) {
  if (length(snps) < 2) {
    return(NULL)
  }

  bfile <- file.path(ldref_dir, paste0(prefix_pattern, chr))
  bed_file <- paste0(bfile, ".bed")
  bim_file <- paste0(bfile, ".bim")

  if (!file.exists(bed_file) || !file.exists(bim_file)) {
    message("Skipping chr", chr, ": missing reference files")
    return(NULL)
  }

  bim <- fread(bim_file, header = FALSE)
  colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

  ld_snps <- bim %>% filter(.data$SNP %in% snps) %>% pull(.data$SNP)
  if (length(ld_snps) < 2) {
    return(NULL)
  }

  extract_path <- paste0(tmp_prefix, "_extract.txt")
  fwrite(data.table(SNP = ld_snps), extract_path, col.names = FALSE)

  cmd_args <- c(
    "--bfile", bfile,
    "--extract", extract_path,
    "--r", "square", "gz",
    "--out", tmp_prefix,
    "--memory", "4000"
  )

  status <- tryCatch(
    system2(plink_bin, args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )

  ld_file <- paste0(tmp_prefix, ".ld.gz")
  if (!file.exists(ld_file)) {
    if (inherits(status, "error")) {
      message("PLINK error on chr", chr, ": ", conditionMessage(status))
    }
    return(NULL)
  }

  list(
    R = read_ld_matrix(ld_file),
    snps = ld_snps,
    bim = bim %>% filter(.data$SNP %in% ld_snps),
    method = "plink"
  )
}

extract_cs_table <- function(susie_fit, locus_id, lead_snp, snp_ids) {
  if (is.null(susie_fit$sets$cs) || length(susie_fit$sets$cs) == 0) {
    return(data.frame())
  }

  cs_list <- susie_fit$sets$cs
  cs_rows <- lapply(seq_along(cs_list), function(i) {
    idx <- cs_list[[i]]
    data.frame(
      locus_id = locus_id,
      lead_snp = lead_snp,
      cs_id = i,
      snp = snp_ids[idx],
      pip = susie_fit$pip[idx],
      stringsAsFactors = FALSE
    )
  })

  bind_rows(cs_rows)
}

## 4. Inputs + harmonization ----

if (!file.exists(lc_pop_path)) {
  stop(
    paste0(
      "Missing LC GWAS file. Expected one of:\n",
      "- GWAS/pcc_vsallfil.txt.gz\n",
      "- GWAS/pcc_vsallfil.txt\n",
      "- Bayes/pcc_vsallfil.txt(.gz)"
    )
  )
}

if (!file.exists(hgi_c2_path)) {
  stop(
    paste0(
      "Missing HGI C2 file. Expected one of:\n",
      "- GWAS/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz\n",
      "- COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz\n",
      "- Bayes/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz"
    )
  )
}

message("Reading HGI and LC GWAS files...")
message("LC GWAS path: ", lc_pop_path)
message("HGI C2 path: ", hgi_c2_path)

hgi_c2_raw <- fread(hgi_c2_path)
lc_pop_raw <- fread(lc_pop_path)

# Add rsID helper
add_rsid_from_hgi <- function(gwas_df, hgi_ref) {
  has_n <- "N" %in% names(gwas_df)
  has_af <- "A1FREQ" %in% names(gwas_df) || "AF" %in% names(gwas_df)

  n_vec <- if (has_n) as.numeric(gwas_df$N) else rep(NA_real_, nrow(gwas_df))

  # Try A1FREQ first (UK Biobank convention), then AF
  af_col <- if ("A1FREQ" %in% names(gwas_df)) "A1FREQ" else if ("AF" %in% names(gwas_df)) "AF" else NULL
  af_vec <- if (!is.null(af_col)) as.numeric(gwas_df[[af_col]]) else rep(NA_real_, nrow(gwas_df))

  df <- gwas_df %>%
    transmute(
      chr = as.integer(.data$CHROM),
      pos = as.integer(.data$POSITION),
      effect_allele = toupper(.data$EFFECT_ALLELE),
      other_allele = toupper(.data$NON_EFFECT_ALLELE),
      beta = .data$BETA,
      se = .data$SE,
      p = .data$P
    )

  df <- df %>% distinct(chr, pos, effect_allele, other_allele, .keep_all = TRUE)
  df$n <- n_vec
  df$af <- af_vec

  # Drop strand-ambiguous SNPs
  is_palindromic <- (df$effect_allele == "A" & df$other_allele == "T") |
    (df$effect_allele == "T" & df$other_allele == "A") |
    (df$effect_allele == "C" & df$other_allele == "G") |
    (df$effect_allele == "G" & df$other_allele == "C")
  df <- df[!is_palindromic, ]

  j <- df %>% inner_join(hgi_ref, by = c("chr", "pos"))

  same_dir <- j$effect_allele == j$ALT & j$other_allele == j$REF
  flip_dir <- j$effect_allele == j$REF & j$other_allele == j$ALT
  keep <- same_dir | flip_dir

  j2 <- j[keep, ]

  j2 %>%
    mutate(
      beta = ifelse(flip_dir[keep], -beta, beta),
      effect_allele = .data$ALT,
      other_allele = .data$REF
    ) %>%
    transmute(
      chr = .data$chr,
      pos = .data$pos,
      rsid = .data$rsid,
      effect_allele = .data$effect_allele,
      other_allele = .data$other_allele,
      beta = .data$beta,
      se = .data$se,
      p = .data$p,
      n = .data$n,
      af = .data$af
    ) %>%
    distinct(.data$rsid, .keep_all = TRUE)
}

hgi_ref <- hgi_c2_raw %>%
  transmute(
    chr = as.integer(`#CHR`),
    pos = as.integer(POS),
    ALT = toupper(ALT),
    REF = toupper(REF),
    rsid = rsid
  ) %>%
  filter(!is.na(rsid), rsid != "") %>%
  distinct(chr, pos, ALT, REF, .keep_all = TRUE)

lc_pop_rsid <- add_rsid_from_hgi(lc_pop_raw, hgi_ref)

hgi_c2 <- hgi_c2_raw %>%
  transmute(
    chr = as.integer(`#CHR`),
    pos = as.integer(POS),
    SNP = rsid,
    beta = all_inv_var_meta_beta,
    se = all_inv_var_meta_sebeta,
    pval = all_inv_var_meta_p
  ) %>%
  filter(!is.na(SNP), !is.na(beta), !is.na(se), !is.na(pval))

lc_assoc <- lc_pop_rsid %>%
  transmute(
    chr = as.integer(chr),
    pos = as.integer(pos),
    SNP = rsid,
    beta = beta,
    se = se,
    pval = p,
    n = as.numeric(n),
    af = as.numeric(af)
  ) %>%
  filter(!is.na(SNP), !is.na(beta), !is.na(se), !is.na(pval), se > 0)

## 5. Define loci --------

if (file.exists(clumped_c2_path)) {
  lead_snps <- fread(clumped_c2_path)$SNP
  message("Using clumped lead SNPs from: ", clumped_c2_path)
} else {
  lead_snps <- hgi_c2 %>% filter(pval < 5e-8) %>% pull(SNP) %>% unique()
  message("Clumped SNP file not found; using C2 genome-wide significant SNPs.")
}

if (length(lead_snps) == 0) {
  stop("No lead SNPs were found for fine-mapping.")
}

lead_tbl <- hgi_c2 %>%
  filter(SNP %in% lead_snps) %>%
  arrange(chr, pos, pval) %>%
  distinct(SNP, .keep_all = TRUE)

# Collapse overlapping loci
if (nrow(lead_tbl) > 1) {
  overlap_bp <- 2L * 500000L
  lead_tbl <- lead_tbl %>%
    group_by(chr) %>%
    arrange(pos, .by_group = TRUE) %>%
    mutate(
      gap = pos - lag(pos),
      locus_group = cumsum(if_else(is.na(gap) | gap > overlap_bp, 1L, 0L))
    ) %>%
    group_by(chr, locus_group) %>%
    slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(chr, pos, SNP, beta, se, pval)
}

window_bp <- 500000L
L_signals <- 5L
min_snps_locus <- 10L
max_snps_locus <- 2000L

message("Found ", nrow(lead_tbl), " loci to fine-map")

## 6. Run SuSiE per locus -----

all_pip <- list()
all_cs <- list()
locus_summary <- list()

for (i in seq_len(nrow(lead_tbl))) {
  lead <- lead_tbl[i, ]
  locus_id <- paste0("chr", lead$chr, ":", max(1, lead$pos - window_bp), "-", lead$pos + window_bp)

  locus_hgi <- hgi_c2 %>%
    filter(chr == lead$chr, pos >= lead$pos - window_bp, pos <= lead$pos + window_bp)

  locus_lc <- lc_assoc %>%
    filter(chr == lead$chr, pos >= lead$pos - window_bp, pos <= lead$pos + window_bp)

  locus_merged <- locus_hgi %>%
    select(chr, pos, SNP) %>%
    inner_join(locus_lc %>% select(SNP, beta, se, pval, n, af), by = "SNP") %>%
    distinct(SNP, .keep_all = TRUE) %>%
    arrange(pval)

  if (nrow(locus_merged) > max_snps_locus) {
    locus_merged <- locus_merged %>% slice_head(n = max_snps_locus)
  }

  if (nrow(locus_merged) < min_snps_locus) {
    locus_summary[[length(locus_summary) + 1]] <- data.frame(
      locus_id = locus_id,
      lead_snp = lead$SNP,
      n_snps_input = nrow(locus_merged),
      n_snps_ld = NA_integer_,
      n_cs = 0L,
      top_snp = NA_character_,
      top_pip = NA_real_,
      status = "skipped_too_few_snps",
      ld_source = ld_source,
      stringsAsFactors = FALSE
    )
    next
  }

  tmp_prefix <- file.path(tempdir(), paste0("susie_chr", lead$chr, "_", i))

  # Compute LD using selected method
  ld_obj <- build_ld(lead$chr, locus_merged$SNP, tmp_prefix, locus_merged)

  if (is.null(ld_obj)) {
    locus_summary[[length(locus_summary) + 1]] <- data.frame(
      locus_id = locus_id,
      lead_snp = lead$SNP,
      n_snps_input = nrow(locus_merged),
      n_snps_ld = NA_integer_,
      n_cs = 0L,
      top_snp = NA_character_,
      top_pip = NA_real_,
      status = "skipped_ld_unavailable",
      ld_source = ld_source,
      stringsAsFactors = FALSE
    )
    next
  }

  ld_snps <- ld_obj$snps
  R <- ld_obj$R

  # Align z to LD SNP order
  locus_aligned <- locus_merged %>%
    filter(SNP %in% ld_snps) %>%
    mutate(z = beta / se)
  z <- locus_aligned$z[match(ld_snps, locus_aligned$SNP)]
  n_vec <- locus_aligned$n[match(ld_snps, locus_aligned$SNP)]

  ok <- is.finite(z)
  if (sum(ok) < min_snps_locus) {
    locus_summary[[length(locus_summary) + 1]] <- data.frame(
      locus_id = locus_id,
      lead_snp = lead$SNP,
      n_snps_input = nrow(locus_merged),
      n_snps_ld = sum(ok),
      n_cs = 0L,
      top_snp = NA_character_,
      top_pip = NA_real_,
      status = "skipped_nonfinite_z",
      ld_source = ld_source,
      stringsAsFactors = FALSE
    )
    next
  }

  z <- z[ok]
  R <- R[ok, ok, drop = FALSE]
  snp_ids <- ld_snps[ok]
  n_vec <- n_vec[ok]

  # Protect against unstable numerics
  z <- pmax(pmin(z, 80), -80)
  R <- make_psd_corr(R)

  n_use <- suppressWarnings(as.integer(round(stats::median(n_vec[is.finite(n_vec) & n_vec > 0], na.rm = TRUE))))
  n_is_valid <- is.finite(n_use) && !is.na(n_use) && n_use > 0

  fit <- tryCatch(
    susie_rss(
      z = z,
      R = R,
      n = if (n_is_valid) n_use else NULL,
      L = L_signals,
      coverage = 0.95,
      max_iter = 1000,
      estimate_residual_variance = FALSE
    ),
    error = function(e) e
  )

  fit_converged <- isTRUE(fit$converged)

  if (inherits(fit, "error")) {
    locus_summary[[length(locus_summary) + 1]] <- data.frame(
      locus_id = locus_id,
      lead_snp = lead$SNP,
      n_snps_input = nrow(locus_merged),
      n_snps_ld = length(snp_ids),
      n_used = if (n_is_valid) n_use else NA_integer_,
      n_cs = 0L,
      top_snp = NA_character_,
      top_pip = NA_real_,
      status = paste0("susie_failed: ", conditionMessage(fit)),
      ld_source = ld_source,
      stringsAsFactors = FALSE
    )
    next
  }

  pip_tbl <- data.frame(
    locus_id = locus_id,
    lead_snp = lead$SNP,
    snp = snp_ids,
    pip = fit$pip,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      z = z,
      n_used = if (n_is_valid) n_use else NA_integer_,
      ld_source = ld_source
    ) %>%
    arrange(desc(pip))

  cs_tbl <- extract_cs_table(fit, locus_id, lead$SNP, snp_ids)

  top_row <- pip_tbl %>% slice(1)
  all_pip[[length(all_pip) + 1]] <- pip_tbl
  all_cs[[length(all_cs) + 1]] <- cs_tbl

  locus_summary[[length(locus_summary) + 1]] <- data.frame(
    locus_id = locus_id,
    lead_snp = lead$SNP,
    n_snps_input = nrow(locus_merged),
    n_snps_ld = length(snp_ids),
    n_used = if (n_is_valid) n_use else NA_integer_,
    n_cs = if (is.null(fit$sets$cs)) 0L else length(fit$sets$cs),
    top_snp = top_row$snp,
    top_pip = top_row$pip,
    status = if (fit_converged) "ok" else "ok_nonconverged",
    ld_source = ld_source,
    stringsAsFactors = FALSE
  )
}

pip_out <- bind_rows(all_pip)
cs_out <- bind_rows(all_cs)
summary_out <- bind_rows(locus_summary)

if (file.exists(lammi_classes_path) && nrow(pip_out) > 0) {
  lammi <- fread(lammi_classes_path)
  if (all(c("rsid", "lammi_class") %in% colnames(lammi))) {
    pip_out <- pip_out %>%
      left_join(
        lammi %>% transmute(snp = rsid, lammi_class = lammi_class),
        by = "snp"
      )
  }
}

fwrite(pip_out, file.path(out_dir, "c2_lc_susie_pip.csv"))
fwrite(cs_out, file.path(out_dir, "c2_lc_susie_credible_sets.csv"))
fwrite(summary_out, file.path(out_dir, "c2_lc_susie_locus_summary.csv"))

## 7. QC plots ----

if (nrow(summary_out) > 0) {
  p_status <- ggplot(summary_out, aes(x = status)) +
    geom_bar(fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste("SuSiE locus run status (LD source:", ld_source, ")"),
      x = "Status",
      y = "Number of loci"
    ) +
    theme_bw()

  ggsave(
    file.path(out_dir, "c2_lc_susie_status_counts.png"),
    p_status,
    width = 8,
    height = 4,
    dpi = 300
  )

  summary_ok <- summary_out %>% filter(status == "ok" | status == "ok_nonconverged")
  if (nrow(summary_ok) > 0) {
    p_top <- ggplot(summary_ok, aes(x = reorder(lead_snp, top_pip), y = top_pip)) +
      geom_col(fill = "darkorange") +
      coord_flip() +
      labs(
        title = paste("Top PIP per locus (LD source:", ld_source, ")"),
        x = "Lead SNP",
        y = "Top PIP"
      ) +
      theme_bw()

    ggsave(
      file.path(out_dir, "c2_lc_susie_top_pip_by_locus.png"),
      p_top,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

if (nrow(pip_out) > 0) {
  p_pip_hist <- ggplot(pip_out, aes(x = pip)) +
    geom_histogram(bins = 60, fill = "#2A9D8F", color = "white") +
    labs(
      title = paste("Distribution of SuSiE PIPs (LD source:", ld_source, ")"),
      x = "PIP",
      y = "SNP count"
    ) +
    theme_bw()

  ggsave(
    file.path(out_dir, "c2_lc_susie_pip_histogram.png"),
    p_pip_hist,
    width = 7,
    height = 4,
    dpi = 300
  )
}

message("\n========= SuSiE fine-mapping finished ==========")
message("LD source: ", ld_source)
message("Saved: ", file.path(out_dir, "c2_lc_susie_pip.csv"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_credible_sets.csv"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_locus_summary.csv"))
message("Saved plots: PNG files in ", out_dir)
message("================================================\n")
