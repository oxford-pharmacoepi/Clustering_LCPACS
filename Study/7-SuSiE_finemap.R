############################################################
## SuSiE fine-mapping for C2 loci on LC GWAS
##
## Inputs:
## - GWAS/pcc_vsallfil.txt.gz (or .txt)
## - GWAS/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz
## - 1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.{bed,bim,fam}
##
## Outputs:
## - bayesResults/susie/c2_lc_susie_pip.csv
## - bayesResults/susie/c2_lc_susie_credible_sets.csv
## - bayesResults/susie/c2_lc_susie_locus_summary.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(susieR)
})

## 1. Paths ---------------------------------------------------------------

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

ldref_dir <- here::here("1000G_EUR_Phase3_plink")
clumped_c2_path <- here::here("Study/clumped_results/c2_clumped_snps.txt")

out_dir <- here::here("bayesResults/susie")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Optional join for downstream interpretation.
lammi_classes_path <- here::here("bayesResults/c2_lc_lammi_classes.csv")

# PLINK binary: repo local first, then PATH.
plink_local <- here::here("plink")
plink_bin <- if (file.exists(plink_local)) plink_local else "plink"

## 2. Helpers -------------------------------------------------------------

add_rsid_from_hgi <- function(gwas_df, hgi_ref) {
  has_n <- "N" %in% names(gwas_df)
  n_vec <- if (has_n) as.numeric(gwas_df$N) else rep(NA_real_, nrow(gwas_df))

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

  # Drop strand-ambiguous SNPs (A/T or C/G) to avoid direction ambiguity.
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
      n = .data$n
    ) %>%
    distinct(.data$rsid, .keep_all = TRUE)
}

make_psd_corr <- function(R, min_eig = 1e-6) {
  R <- (R + t(R)) / 2
  diag(R) <- 1

  eig <- eigen(R, symmetric = TRUE)
  eig$values[eig$values < min_eig] <- min_eig
  R_psd <- eig$vectors %*% (eig$values * t(eig$vectors))

  # Re-scale to correlation form.
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

build_plink_ld <- function(chr, snps, tmp_prefix) {
  if (length(snps) < 2) {
    return(NULL)
  }

  bfile <- file.path(ldref_dir, paste0("1000G.EUR.QC.", chr))
  bed_file <- paste0(bfile, ".bed")
  bim_file <- paste0(bfile, ".bim")

  if (!file.exists(bed_file) || !file.exists(bim_file)) {
    message("Skipping chr", chr, ": missing reference files in ", ldref_dir)
    return(NULL)
  }

  bim <- fread(bim_file, header = FALSE)
  colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

  # Keep SNPs in genotype order to align with PLINK LD output.
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
    bim = bim %>% filter(.data$SNP %in% ld_snps)
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

## 3. Inputs + harmonization ---------------------------------------------

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
    n = as.numeric(n)
  ) %>%
  filter(!is.na(SNP), !is.na(beta), !is.na(se), !is.na(pval), se > 0)

## 4. Define loci ---------------------------------------------------------

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

# Collapse overlapping lead windows so each locus is fine-mapped once.
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
max_snps_locus <- 5000L

## 5. Run SuSiE per locus -------------------------------------------------

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

  # Restrict to SNPs observed in both summary datasets.
  locus_merged <- locus_hgi %>%
    select(chr, pos, SNP) %>%
    inner_join(locus_lc %>% select(SNP, beta, se, pval, n), by = "SNP") %>%
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
      stringsAsFactors = FALSE
    )
    next
  }

  tmp_prefix <- file.path(tempdir(), paste0("susie_chr", lead$chr, "_", i))
  ld_obj <- build_plink_ld(lead$chr, locus_merged$SNP, tmp_prefix)

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
      stringsAsFactors = FALSE
    )
    next
  }

  ld_snps <- ld_obj$snps
  R <- ld_obj$R

  # Align z to LD SNP order.
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
      stringsAsFactors = FALSE
    )
    next
  }

  z <- z[ok]
  R <- R[ok, ok, drop = FALSE]
  snp_ids <- ld_snps[ok]
  n_vec <- n_vec[ok]

  # Protect against unstable numerics from LD/stat mismatch.
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
      n_used = if (n_is_valid) n_use else NA_integer_
    ) %>%
    left_join(ld_obj$bim %>% transmute(snp = SNP, bp = BP), by = "snp") %>%
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

## 6. QC plots ------------------------------------------------------------

if (nrow(summary_out) > 0) {
  p_status <- ggplot(summary_out, aes(x = status)) +
    geom_bar(fill = "steelblue") +
    coord_flip() +
    labs(
      title = "SuSiE locus run status",
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

  summary_ok <- summary_out %>% filter(status == "ok")
  if (nrow(summary_ok) > 0) {
    p_top <- ggplot(summary_ok, aes(x = reorder(lead_snp, top_pip), y = top_pip)) +
      geom_col(fill = "darkorange") +
      coord_flip() +
      labs(
        title = "Top PIP per locus",
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
      title = "Distribution of SuSiE PIPs",
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

message("SuSiE fine-mapping finished.")
message("Saved: ", file.path(out_dir, "c2_lc_susie_pip.csv"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_credible_sets.csv"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_locus_summary.csv"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_status_counts.png"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_top_pip_by_locus.png"))
message("Saved: ", file.path(out_dir, "c2_lc_susie_pip_histogram.png"))