############################################################
## MR helper functions used by 2.mr_run_all.R
############################################################

setup_parallel_with_progress <- function(parallel_jobs = 1,
                                         max_global_size_gb = 50) {
  library(future)
  library(progressr)

  is_interactive <- interactive()

  if (parallel_jobs > 1) {
    if (is_interactive) {
      handlers(list(
        handler_progress(
          format = ":spin [:bar] :percent | :current/:total | Elapsed: :elapsedfull | ETA: :eta",
          width = 100,
          complete = "=",
          incomplete = "-",
          clear = FALSE
        )
      ))
    } else {
      handlers(list(
        handler_progress(
          format = "[-] :current/:total (:percent) | Elapsed: :elapsedfull | ETA: :eta | :message",
          clear = FALSE
        )
      ))
    }

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }

    options(future.globals.maxSize = max_global_size_gb * 1024^3)
    plan(multicore, workers = parallel_jobs)
  } else {
    handlers("progress")
    plan(sequential)
  }

  invisible(NULL)
}

get_f <- function(dat, trait1, trait2, F_value = 10) {
  log <- is.na(dat$eaf.exposure)
  log <- unique(log)

  if (length(log) == 1) {
    if (log == TRUE) {
      print("Frequency is not found, F-value cannot be calculated.")
      return(dat)
    }
  }

  if (is.null(dat$beta.exposure[1]) == TRUE || is.na(dat$beta.exposure[1]) == TRUE) {
    print("BETA is not found, F-value cannot be calculated.")
    return(dat)
  }

  if (is.null(dat$se.exposure[1]) == TRUE || is.na(dat$se.exposure[1]) == TRUE) {
    print(paste0("SE is not found, F-value cannot be calculated.", trait1, trait2))
    return(dat)
  }

  if (is.null(dat$samplesize.exposure[1]) == TRUE || is.na(dat$samplesize.exposure[1]) == TRUE) {
    print("Sample size is not found, F-value cannot be calculated.")
    return(dat)
  }

  if ("FALSE" %in% log &&
      is.null(dat$beta.exposure[1]) == FALSE && is.na(dat$beta.exposure[1]) == FALSE &&
      is.null(dat$se.exposure[1]) == FALSE && is.na(dat$se.exposure[1]) == FALSE &&
      is.null(dat$samplesize.exposure[1]) == FALSE && is.na(dat$samplesize.exposure[1]) == FALSE) {
    R2 <- (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) /
      ((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) +
         (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$se.exposure^2) * dat$samplesize.exposure))
    F <- (dat$samplesize.exposure - 2) * R2 / (1 - R2)
    dat$R2 <- R2
    dat$F <- F
    dat <- subset(dat, F > F_value)
    return(dat)
  }
}

mr_single <- function(exp_dat, outcome_dat, trait1, trait2) {
  mr_result <- c(
    trait1 = trait1,
    trait2 = trait2,
    nsnp = 0,
    method = "none",
    beta = -99,
    se = -99,
    beta_low95 = -99,
    beta_up95 = -99,
    pval = -99
  )

  outcome_dat <- outcome_dat |>
    filter(SNP %in% exp_dat$SNP)

  if (nrow(outcome_dat) == 0) {
    return(mr_result)
  }

  suppressMessages(dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = outcome_dat,
    action = 2
  ))

  dat1 <- get_f(dat, trait1 = trait1, trait2 = trait2)

  if (nrow(dat1) == 0) {
    method <- "none"
  } else {
    dat2 <- dat1 |>
      dplyr::filter(mr_keep == TRUE)
    if (nrow(dat2) == 0) {
      method <- "none"
    } else if (nrow(dat2) == 1) {
      method <- "mr_wald_ratio"
    } else {
      method <- "mr_ivw"
    }
  }

  if (method %in% c("mr_wald_ratio", "mr_ivw")) {
    tryCatch({
      suppressMessages(tsmr <- mr(dat1, method_list = c(method)))
      mr_result["nsnp"] <- tsmr$nsnp[1]
      mr_result["method"] <- method
      mr_result["beta"] <- tsmr$b[1]
      mr_result["se"] <- tsmr$se[1]
      mr_result["beta_low95"] <- tsmr$b[1] - 1.96 * tsmr$se[1]
      mr_result["beta_up95"] <- tsmr$b[1] + 1.96 * tsmr$se[1]
      mr_result["pval"] <- tsmr$pval
    }, error = function(e) {
      message(paste0("MR analysis failed: ", trait1, trait2, " | Error: ", e$message))
      mr_result["nsnp"] <- -98
      mr_result["method"] <- "error"
      mr_result["beta"] <- -98
      mr_result["se"] <- -98
      mr_result["beta_low95"] <- -98
      mr_result["beta_up95"] <- -98
      mr_result["pval"] <- -98
    })
  }

  mr_result
}

############################################################
## Coloc helper functions
############################################################

run_coloc <- function(dataset1, dataset2, trait1_name, outcome_diag, lead_snp) {
  tryCatch({
    res <- coloc::coloc.abf(dataset1 = dataset1, dataset2 = dataset2)
    s <- res$summary
    snp_results <- res$results
    top_h4_snp <- snp_results$snp[which.max(snp_results$SNP.PP.H4)]

    data.frame(
      trait1 = trait1_name,
      trait2 = outcome_diag,
      lead_snp = lead_snp,
      nsnp = s["nsnps"],
      PP.H0 = s["PP.H0.abf"],
      PP.H1 = s["PP.H1.abf"],
      PP.H2 = s["PP.H2.abf"],
      PP.H3 = s["PP.H3.abf"],
      PP.H4 = s["PP.H4.abf"],
      top_H4_snp = top_h4_snp,
      row.names = NULL
    )
  }, error = function(e) {
    warning(sprintf("Coloc failed: %s -> %s [%s]: %s",
                    trait1_name, outcome_diag, lead_snp, e$message))
    NULL
  })
}

align_alleles <- function(dat1, dat2,
                          beta_col1 = "BETA", beta_col2 = "beta",
                          se_col1 = "SE", se_col2 = "sebeta",
                          af_col1 = "AF", af_col2 = "AF",
                          a1_col1 = "A1", a2_col1 = "A2",
                          a1_col2 = "A1", a2_col2 = "A2",
                          chr_col1 = "CHR", bp_col1 = "BP_hg37",
                          chr_col2 = "CHR", bp_col2 = "BP_hg37") {
  dat1 <- data.table::copy(data.table::as.data.table(dat1))
  dat2 <- data.table::copy(data.table::as.data.table(dat2))

  data.table::setnames(
    dat1,
    c(beta_col1, se_col1, af_col1, a1_col1, a2_col1, chr_col1, bp_col1),
    c("BETA", "SE", "AF", "A1", "A2", "CHR", "BP")
  )
  data.table::setnames(
    dat2,
    c(beta_col2, se_col2, af_col2, a1_col2, a2_col2, chr_col2, bp_col2),
    c("BETA", "SE", "AF", "A1", "A2", "CHR", "BP")
  )

  dat1[, CHR := as.character(CHR)]
  dat2[, CHR := as.character(CHR)]
  dat1[, snp_key := paste(CHR, BP, A1, A2, sep = ":")]
  dat2[, snp_key := paste(CHR, BP, A1, A2, sep = ":")]
  dat2[, snp_key_fl := paste(CHR, BP, A2, A1, sep = ":")]

  forward <- dat2[snp_key %in% dat1$snp_key]
  to_flip <- dat2[snp_key_fl %in% dat1$snp_key & !snp_key %in% dat1$snp_key]
  if (nrow(to_flip) > 0) {
    to_flip[, BETA := -BETA]
    to_flip[, AF := 1 - AF]
    to_flip[, A1_new := A2]
    to_flip[, A2 := A1]
    to_flip[, A1 := A1_new]
    to_flip[, A1_new := NULL]
    to_flip[, snp_key := snp_key_fl]
  }

  dat2_aligned <- data.table::rbindlist(list(forward, to_flip), fill = TRUE)
  dat2_aligned[, snp_key_fl := NULL]

  dedup <- function(dat) {
    dat[, MAF := pmin(AF, 1 - AF)]
    dat[order(-MAF)][!duplicated(snp_key)]
  }

  list(dat1 = dedup(dat1), dat2 = dedup(dat2_aligned))
}

run_one_coloc_task <- function(
  trait1_rawcol, trait1_name, file_class, lead_snp,
  outcome_diag, Nca, Nco,
  finn_with_rs, coloc_sub, data_list, window_size
) {
  exp_coloc_region <- coloc_sub
  if (is.null(exp_coloc_region)) return(NULL)

  exp_info <- exp_coloc_region[["INFO"]]
  exp_type <- exp_info[["type"]]
  exp_N_global <- as.numeric(exp_info[["N"]])
  exp_s <- if ("s" %in% names(exp_info)) exp_info[["s"]] else NA

  if (file_class %in% c("metab", "IDP", "panukb")) {
    exp_region <- exp_coloc_region[[lead_snp]]
    if (is.null(exp_region) || nrow(exp_region) == 0) return(NULL)

    lead_snps_dat <- data_list[[trait1_rawcol]]
    lead_row <- lead_snps_dat[SNP == lead_snp]
    if (nrow(lead_row) == 0) return(NULL)

    lead_chr <- lead_row$CHR[1]
    lead_bp <- lead_row$BP[1]

    finn_region <- finn_with_rs[
      CHR == lead_chr &
        BP_hg37 >= (lead_bp - window_size) &
        BP_hg37 <= (lead_bp + window_size)
    ]
  } else if (file_class == "pqtl") {
    exp_region <- exp_coloc_region[[trait1_name]]
    if (is.null(exp_region) || nrow(exp_region) == 0) return(NULL)

    exp_chr <- exp_region$CHR[1]
    exp_bp_min <- min(exp_region$BP_hg37)
    exp_bp_max <- max(exp_region$BP_hg37)

    finn_region <- finn_with_rs[
      CHR == exp_chr &
        BP_hg37 >= exp_bp_min &
        BP_hg37 <= exp_bp_max
    ]
  } else {
    return(NULL)
  }

  if (nrow(finn_region) == 0) return(NULL)

  aligned <- align_alleles(
    dat1 = exp_region,
    dat2 = finn_region,
    beta_col1 = "BETA", beta_col2 = "beta",
    se_col1 = "SE", se_col2 = "sebeta",
    af_col1 = "AF", af_col2 = "af_alt",
    a1_col1 = "A1", a1_col2 = "alt",
    a2_col1 = "A2", a2_col2 = "ref",
    chr_col1 = "CHR", chr_col2 = "CHR",
    bp_col1 = "BP_hg37", bp_col2 = "BP_hg37"
  )
  exp_aln <- aligned$dat1
  finn_aln <- aligned$dat2

  common_snps <- intersect(exp_aln$snp_key, finn_aln$snp_key)
  if (length(common_snps) < 5) return(NULL)

  exp_aln <- exp_aln[snp_key %in% common_snps]
  finn_aln <- finn_aln[snp_key %in% common_snps]
  data.table::setkey(exp_aln, snp_key)
  data.table::setkey(finn_aln, snp_key)

  dataset1 <- list(
    snp = exp_aln$snp_key,
    beta = exp_aln$BETA,
    varbeta = exp_aln$varBETA,
    MAF = exp_aln$MAF,
    type = exp_type,
    N = if (exp_N_global == -1) exp_aln$N else exp_N_global
  )
  if (exp_type == "cc" && !is.null(exp_s) && !is.na(exp_s)) {
    dataset1$s <- as.numeric(exp_s)
  }

  dataset2 <- list(
    snp = finn_aln$snp_key,
    beta = finn_aln$BETA,
    varbeta = finn_aln$SE^2,
    MAF = finn_aln$MAF,
    type = "cc",
    N = Nca + Nco,
    s = Nca / (Nca + Nco)
  )

  invisible(capture.output(
    res_coloc <- run_coloc(dataset1, dataset2, trait1_name, outcome_diag, lead_snp)
  ))
  res_coloc
}

infer_coloc_file_class <- function(df) {
  if (!"file_class" %in% names(df)) {
    df$file_class <- NA_character_
  }

  df |>
    dplyr::mutate(
      file_class = dplyr::case_when(
        !is.na(.data$file_class) ~ .data$file_class,
        .data$exp_cate == "Proteome" ~ "pqtl",
        .data$exp_cate == "Metabolome" ~ "metab",
        .data$exp_cate %in% c(
          "IDP", "Brain imaging", "Brain image",
          "Imaging-derived phenotypes",
          "Magnetic resonance imaging phenotypes"
        ) ~ "IDP",
        grepl("^[0-9]+_p[0-9]+_i[0-9]+$", .data$trait1_rawcol) ~ "IDP",
        grepl("^p[0-9]+_i[0-9]+$", .data$trait1) ~ "IDP",
        TRUE ~ "panukb"
      )
    )
}

add_coloc_outcome_diag <- function(df) {
  if ("outcome_diag" %in% names(df)) {
    return(df)
  }

  outcome_col <- dplyr::case_when(
    "trait2" %in% names(df) ~ "trait2",
    "disease_code" %in% names(df) ~ "disease_code",
    "disease_code.x" %in% names(df) ~ "disease_code.x",
    "disease_code.y" %in% names(df) ~ "disease_code.y",
    TRUE ~ NA_character_
  )

  if (is.na(outcome_col)) {
    stop("MR result has no outcome column: expected trait2 or disease_code.")
  }

  df$outcome_diag <- df[[outcome_col]]
  df
}

read_coloc_mr_result <- function(mr_file) {
  if (!file.exists(mr_file)) {
    stop("MR result file not found: ", mr_file)
  }

  ext <- tolower(tools::file_ext(mr_file))
  result <- switch(
    ext,
    "csv" = data.table::fread(mr_file) |> as.data.frame(),
    "xlsx" = readxl::read_xlsx(mr_file),
    "xls" = readxl::read_xlsx(mr_file),
    "rds" = readRDS(mr_file),
    stop("Unsupported MR result file type: ", ext)
  )

  result <- infer_coloc_file_class(result)
  result <- add_coloc_outcome_diag(result)

  if ("pval" %in% names(result)) {
    result <- result |>
      dplyr::filter(is.na(.data$pval) | .data$pval >= 0)
  }

  result
}

get_coloc_exposure_table <- function(result, target_class) {
  required_cols <- c("trait1", "trait1_rawcol", "file_class")
  missing_cols <- setdiff(required_cols, names(result))
  if (length(missing_cols) > 0) {
    stop("MR result is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  result |>
    dplyr::distinct(.data$trait1, .keep_all = TRUE) |>
    dplyr::filter(.data$file_class %in% target_class)
}

print_coloc_exposure_summary <- function(result, result_in_loop, target_class) {
  cat("  Total MR rows:", nrow(result), "\n")
  cat("  Exposure classes:\n")
  print(table(result$file_class, useNA = "ifany"))
  cat("  Target class:", paste(target_class, collapse = ", "), "\n")
  cat("  Exposures to format:", nrow(result_in_loop), "\n")
}

get_leadsnp_dt <- function(data_list, trait1_rawcol) {
  if (!trait1_rawcol %in% names(data_list)) {
    warning(sprintf("trait1_rawcol not found in exp_data_list: %s", trait1_rawcol))
    return(NULL)
  }

  leadsnp_dt <- data.table::as.data.table(data_list[[trait1_rawcol]])
  if (nrow(leadsnp_dt) == 0) {
    warning(sprintf("No lead SNPs for %s", trait1_rawcol))
    return(NULL)
  }

  required_cols <- c("SNP", "CHR", "BP")
  missing_cols <- setdiff(required_cols, names(leadsnp_dt))
  if (length(missing_cols) > 0) {
    warning(sprintf(
      "Lead SNP table %s is missing columns: %s",
      trait1_rawcol, paste(missing_cols, collapse = ", ")
    ))
    return(NULL)
  }

  leadsnp_dt
}

save_coloc_region_list <- function(region_list, region_names, info, outpath) {
  region_list <- stats::setNames(region_list, region_names)
  region_list$INFO <- info
  qs::qsave(region_list, outpath)
  invisible(outpath)
}

get_panukb_sample_info <- function(header_names, manifest_row) {
  af_col_cc <- c("af_cases_EUR", "af_controls_EUR")
  af_col_quant <- "af_EUR"

  if (all(af_col_cc %in% header_names)) {
    n_case <- as.numeric(manifest_row$n_cases_EUR[1L])
    n_control <- as.numeric(manifest_row$n_controls_EUR[1L])
    n_total <- n_case + n_control

    return(list(
      af_cols = af_col_cc,
      trait_type = "cc",
      N = n_total,
      Nca = n_case,
      Nco = n_control,
      s = n_case / n_total
    ))
  }

  if (af_col_quant %in% header_names) {
    return(list(
      af_cols = af_col_quant,
      trait_type = "quant",
      N = as.numeric(manifest_row$n_cases_EUR[1L]),
      Nca = NA_real_,
      Nco = NA_real_,
      s = NA_real_
    ))
  }

  stop("Cannot find Pan-UKB EUR allele frequency columns.")
}

format_lead_window_regions <- function(dt, leadsnp_dt, window_size, format_region_fun) {
  cat("|| Begin looping ... Num:", nrow(leadsnp_dt), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  progressr::with_progress({
    p <- progressr::progressor(along = seq_len(nrow(leadsnp_dt)))

    future.apply::future_lapply(seq_len(nrow(leadsnp_dt)), function(j) {
      snp_id <- leadsnp_dt$SNP[j]
      snp_chr <- as.integer(leadsnp_dt$CHR[j])
      snp_bp <- as.integer(leadsnp_dt$BP[j])

      region_dt <- dt[
        CHR == snp_chr &
          BP_hg37 %between% c(snp_bp - window_size, snp_bp + window_size)
      ]

      region_dt <- format_region_fun(region_dt)
      p(sprintf("SNP %s | chr%s", snp_id, snp_chr))
      region_dt
    }, future.seed = TRUE, future.scheduling = 1)
  })
}

format_panukb_one <- function(filename0, trait1_rawcol, data_list, pan_ukb,
                              input_dir, output_dir, window_size) {
  outpath <- glue::glue("{output_dir}/panukb_{filename0}.qs")

  leadsnp_dt <- get_leadsnp_dt(data_list, trait1_rawcol)
  if (is.null(leadsnp_dt)) return(NULL)

  manifest_row <- pan_ukb[filename == paste0(filename0, ".tsv.bgz")]
  if (nrow(manifest_row) == 0) {
    warning(sprintf("Pan-UKB manifest row not found: %s.tsv.bgz", filename0))
    return(NULL)
  }

  f_path <- glue::glue("{input_dir}/{filename0}.tsv.bgz")
  if (!file.exists(f_path)) {
    warning(sprintf("Pan-UKB summary file not found: %s", f_path))
    return(NULL)
  }

  header <- data.table::fread(f_path, nrows = 0L)
  sample_info <- get_panukb_sample_info(names(header), manifest_row)

  required_cols <- c(
    "chr", "pos", "ref", "alt",
    "beta_EUR", "se_EUR",
    "neglog10_pval_EUR", "low_confidence_EUR",
    sample_info$af_cols
  )
  missing_cols <- setdiff(required_cols, names(header))
  if (length(missing_cols) > 0) {
    warning(sprintf("Missing columns for %s: %s", filename0, paste(missing_cols, collapse = ", ")))
    return(NULL)
  }

  dt <- data.table::fread(f_path, select = required_cols)
  data.table::setnames(
    dt,
    old = c("chr", "pos", "ref", "alt", "beta_EUR", "se_EUR"),
    new = c("CHR", "BP_hg37", "A2", "A1", "BETA", "SE")
  )
  dt[, `:=`(
    CHR = as.integer(CHR),
    BP_hg37 = as.integer(BP_hg37),
    P = pmax(10^(-neglog10_pval_EUR), .Machine$double.xmin),
    varBETA = SE * SE
  )]
  dt <- dt[!is.na(BETA) & !is.na(SE) & !is.na(CHR) & !is.na(BP_hg37)]
  data.table::setkey(dt, CHR, BP_hg37)

  format_region <- function(region_dt) {
    if (sample_info$trait_type == "cc") {
      region_dt[, `:=`(
        AF = (af_cases_EUR * sample_info$Nca + af_controls_EUR * sample_info$Nco) / sample_info$N,
        SNP = paste(CHR, BP_hg37, sep = ":")
      )]
    } else {
      region_dt[, `:=`(
        AF = af_EUR,
        SNP = paste(CHR, BP_hg37, sep = ":")
      )]
    }
    region_dt[, MAF := pmin(AF, 1 - AF)]
    region_dt[
      !is.na(AF) & !is.na(MAF),
      .(SNP, CHR, BP_hg37, A1, A2, BETA, SE, varBETA, P, AF, MAF)
    ]
  }

  region_list <- format_lead_window_regions(dt, leadsnp_dt, window_size, format_region)
  save_coloc_region_list(
    region_list,
    leadsnp_dt$SNP,
    info = c(type = sample_info$trait_type, N = sample_info$N, s = sample_info$s),
    outpath = outpath
  )
}

format_metab_one <- function(filename0, trait1_rawcol, data_list,
                             input_dir, output_dir, window_size) {
  outpath <- glue::glue("{output_dir}/metab_{filename0}.qs")

  leadsnp_dt <- get_leadsnp_dt(data_list, trait1_rawcol)
  if (is.null(leadsnp_dt)) return(NULL)

  f_path <- glue::glue("{input_dir}/total_{filename0}_GWAS.linear")
  if (!file.exists(f_path)) {
    warning(sprintf("Metabolomics summary file not found: %s", f_path))
    return(NULL)
  }

  dt <- data.table::fread(f_path)
  required_cols <- c("#CHROM", "POS", "ID", "AX", "A1", "A1_FREQ", "BETA", "SE", "P", "OBS_CT")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    warning(sprintf("Missing columns for %s: %s", filename0, paste(missing_cols, collapse = ", ")))
    return(NULL)
  }

  data.table::setnames(
    dt,
    old = c("#CHROM", "POS", "ID", "AX", "A1_FREQ", "OBS_CT"),
    new = c("CHR", "BP_hg37", "SNP", "A2", "AF", "N")
  )
  dt[, `:=`(
    CHR = as.integer(CHR),
    BP_hg37 = as.integer(BP_hg37),
    P = pmax(as.numeric(P), .Machine$double.xmin),
    varBETA = SE * SE
  )]
  dt <- dt[!is.na(BETA) & !is.na(SE) & !is.na(CHR) & !is.na(BP_hg37)]
  data.table::setkey(dt, CHR, BP_hg37)

  format_region <- function(region_dt) {
    region_dt[, MAF := pmin(AF, 1 - AF)]
    region_dt[
      !is.na(AF) & !is.na(MAF),
      .(SNP, CHR, BP_hg37, A1, A2, BETA, SE, varBETA, P, AF, MAF, N)
    ]
  }

  region_list <- format_lead_window_regions(dt, leadsnp_dt, window_size, format_region)
  save_coloc_region_list(
    region_list,
    leadsnp_dt$SNP,
    info = c(type = "quant", N = -1, s = NA_real_),
    outpath = outpath
  )
}

format_idp_one <- function(filename0, trait1_rawcol, data_list, idp_info, idp_snplist,
                           input_dir, output_dir, window_size) {
  outpath <- glue::glue("{output_dir}/IDP_{filename0}.qs")

  leadsnp_dt <- get_leadsnp_dt(data_list, trait1_rawcol)
  if (is.null(leadsnp_dt)) return(NULL)

  idp_fid <- stringr::str_extract(filename0, "(?<=p)[0-9]+")
  idp_idx <- which(idp_info$`UKB ID` == idp_fid)[1]
  if (is.na(idp_idx)) {
    warning(sprintf("IDP UKB ID not found in idp_info: %s", filename0))
    return(NULL)
  }

  idp_number <- stringr::str_pad(as.character(idp_info$Pheno[idp_idx]), width = 4, side = "left", pad = "0")
  f_path <- glue::glue("{input_dir}/{idp_number}.txt")
  if (!file.exists(f_path)) {
    warning(sprintf("IDP summary file not found: %s", f_path))
    return(NULL)
  }

  dt <- data.table::fread(f_path)
  dt <- dt[chr != "0X"]
  data.table::setkey(dt, chr, rsid, pos, a1, a2)
  dt <- idp_snplist[dt, nomatch = 0]
  dt[, `:=`(
    chr = as.numeric(chr),
    P = pmax(10^(-`pval(-log10)`), .Machine$double.xmin)
  )]

  data.table::setnames(
    dt,
    old = c("chr", "pos", "rsid", "a1", "a2", "af", "beta", "se"),
    new = c("CHR", "BP_hg37", "SNP", "A1", "A2", "AF", "BETA", "SE")
  )
  dt[, `:=`(
    CHR = as.integer(CHR),
    BP_hg37 = as.integer(BP_hg37),
    varBETA = SE * SE
  )]
  dt <- dt[!is.na(BETA) & !is.na(SE) & !is.na(CHR) & !is.na(BP_hg37)]
  data.table::setkey(dt, CHR, BP_hg37)

  format_region <- function(region_dt) {
    region_dt[, MAF := pmin(AF, 1 - AF)]
    region_dt[
      !is.na(AF) & !is.na(MAF),
      .(SNP, CHR, BP_hg37, A1, A2, BETA, SE, varBETA, P, AF, MAF)
    ]
  }

  region_list <- format_lead_window_regions(dt, leadsnp_dt, window_size, format_region)
  save_coloc_region_list(
    region_list,
    leadsnp_dt$SNP,
    info = c(type = "quant", N = 33224, s = NA_real_),
    outpath = outpath
  )
}

format_protein_one <- function(protein, protein_map, pgwas_dir, output_dir,
                               window_size, work_dir = here::here()) {
  outpath <- glue::glue("{output_dir}/pqtl_{protein}.qs")

  if (file.exists(outpath)) {
    cat("Skip existing:", outpath, "\n")
    return(outpath)
  }

  pm_rows <- protein_map[Assay == protein]
  if (nrow(pm_rows) == 0) {
    warning(sprintf("Protein %s not found in protein_map, skipping", protein))
    return(NULL)
  }

  pm_info <- pm_rows[1]
  chr_num <- pm_info$chr
  gene_start <- pm_info$gene_start
  gene_end <- pm_info$gene_end
  uniprot <- pm_info$UniProt
  oid <- pm_info$OlinkID
  panel <- pm_info$Panel

  dir_name <- glue::glue("{protein}_{uniprot}_{oid}_v1_{panel}")
  protein_dir <- glue::glue("{pgwas_dir}/{dir_name}")

  if (dir_name == "DFFA_O00273_OID20620_v1_Inflammation") {
    protein_dir <- glue::glue("{work_dir}/rawdata/{dir_name}")
  }

  if (!dir.exists(protein_dir)) {
    warning(sprintf("GWAS dir not found: %s", protein_dir))
    return(NULL)
  }

  chr_file <- glue::glue("{protein_dir}/discovery_chr{chr_num}_{protein}:{uniprot}:{oid}:v1:{panel}.gz")
  if (!file.exists(chr_file)) {
    warning(sprintf("GWAS file not found: %s", chr_file))
    return(NULL)
  }

  dt <- tryCatch(
    data.table::fread(chr_file),
    error = function(e) {
      warning(sprintf("Error reading %s: %s", chr_file, e$message))
      NULL
    }
  )
  if (is.null(dt)) return(NULL)

  region_start <- gene_start - window_size
  region_end <- gene_end + window_size
  dt_region <- dt[GENPOS >= region_start & GENPOS <= region_end]
  if (nrow(dt_region) == 0) {
    warning(sprintf("No SNPs in region for %s, skipping", protein))
    return(NULL)
  }

  dt_region[, c("CHR", "BP_hg37") := data.table::tstrsplit(ID, ":", keep = 1:2)]
  dt_region[, `:=`(
    CHR = as.integer(CHR),
    BP_hg37 = as.integer(BP_hg37)
  )]

  data.table::setnames(
    dt_region,
    old = c("ALLELE0", "ALLELE1", "A1FREQ", "LOG10P"),
    new = c("A2", "A1", "AF", "LOG10P"),
    skip_absent = TRUE
  )

  dt_region[, `:=`(
    SNP = ID,
    varBETA = SE * SE,
    MAF = pmin(AF, 1 - AF),
    P = pmax(10^(-LOG10P), .Machine$double.xmin)
  )]

  dt_region <- dt_region[
    !is.na(AF) & !is.na(MAF),
    .(SNP, CHR, BP_hg37, A1, A2, BETA, SE, varBETA, P, AF, MAF, N)
  ]

  trait1_result_list <- list()
  trait1_result_list[[protein]] <- dt_region
  trait1_result_list$INFO <- c(type = "quant", N = -1, s = NA_real_)

  qs::qsave(trait1_result_list, outpath)
  outpath
}

get_coloc_region_path <- function(coloc_dir, file_class, trait1_name) {
  file.path(coloc_dir, sprintf("%s_%s.qs", file_class, trait1_name))
}

preload_coloc_exposure_regions <- function(params_pairs_all, coloc_dir) {
  exposure_files <- params_pairs_all |>
    dplyr::select(file_class, trait1_name) |>
    dplyr::distinct()
  exposure_files$coloc_file <- mapply(
    get_coloc_region_path,
    coloc_dir = coloc_dir,
    file_class = exposure_files$file_class,
    trait1_name = exposure_files$trait1_name,
    USE.NAMES = FALSE
  )

  missing_files <- exposure_files$coloc_file[!file.exists(exposure_files$coloc_file)]
  if (length(missing_files) > 0) {
    warning(sprintf("Missing coloc exposure region files: %d", length(missing_files)))
    warning(paste(utils::head(missing_files, 10), collapse = "\n"))
  }

  exposure_files <- exposure_files[file.exists(exposure_files$coloc_file), , drop = FALSE]
  cache <- vector("list", nrow(exposure_files))
  names(cache) <- exposure_files$coloc_file

  cat(sprintf("|| Preloading exposure region files: %d\n", nrow(exposure_files)))
  for (i in seq_len(nrow(exposure_files))) {
    cache[[i]] <- qs::qread(exposure_files$coloc_file[i])
    if (i %% 25L == 0L || i == nrow(exposure_files)) {
      cat(sprintf("||   loaded %d/%d exposure files\n", i, nrow(exposure_files)))
      gc(verbose = FALSE)
    }
  }

  cache
}

build_coloc_params_for_outcome <- function(params_pairs_outcome,
                                           data_list,
                                           coloc_dir,
                                           exposure_region_cache = NULL) {
  if (nrow(params_pairs_outcome) == 0) {
    return(data.frame())
  }

  params_list <- vector("list", nrow(params_pairs_outcome))

  for (i in seq_len(nrow(params_pairs_outcome))) {
    trait1_rawcol <- params_pairs_outcome$trait1_rawcol[i]
    trait1_name <- params_pairs_outcome$trait1_name[i]
    file_class <- params_pairs_outcome$file_class[i]
    coloc_file <- get_coloc_region_path(coloc_dir, file_class, trait1_name)

    if (!file.exists(coloc_file)) {
      warning(sprintf("Coloc exposure region file not found: %s", coloc_file))
      next
    }

    if (!is.null(exposure_region_cache)) {
      exp_coloc_region <- exposure_region_cache[[coloc_file]]
      if (is.null(exp_coloc_region)) {
        warning(sprintf("Coloc exposure region not found in cache: %s", coloc_file))
        next
      }
    } else {
      exp_coloc_region <- qs::qread(coloc_file)
    }

    if (file_class %in% c("metab", "IDP", "panukb")) {
      lead_snps_dat <- data_list[[trait1_rawcol]]
      if (is.null(lead_snps_dat) || nrow(lead_snps_dat) == 0) {
        warning(sprintf("No lead SNPs in exp_data_list for %s", trait1_rawcol))
        next
      }
      lead_snps <- intersect(lead_snps_dat$SNP, setdiff(names(exp_coloc_region), "INFO"))
    } else if (file_class == "pqtl") {
      lead_snps <- trait1_name
    } else {
      next
    }

    if (length(lead_snps) == 0) {
      warning(sprintf("No coloc regions available for %s_%s", file_class, trait1_name))
      next
    }

    params_i <- data.frame(
      outcome_diag = params_pairs_outcome$outcome_diag[i],
      trait1_rawcol = trait1_rawcol,
      trait1_name = trait1_name,
      file_class = file_class,
      lead_snp = lead_snps,
      stringsAsFactors = FALSE
    )
    if (is.null(exposure_region_cache)) {
      params_i$coloc_sub <- rep(list(exp_coloc_region), nrow(params_i))
    } else {
      params_i$coloc_key <- coloc_file
    }
    params_list[[i]] <- params_i
  }

  dplyr::bind_rows(params_list)
}

read_finngen_outcome_for_coloc <- function(outcome_diag, diag_info, finn_sum_dir, finn_snplist) {
  diag_idx <- which(diag_info$NAME == outcome_diag)[1]
  if (is.na(diag_idx)) {
    warning(sprintf("Outcome %s not found in diag_info", outcome_diag))
    return(NULL)
  }

  outcome_file <- diag_info$gwas_summary[diag_idx]
  Nca <- diag_info$num_cases[diag_idx]
  Nco <- diag_info$num_controls[diag_idx]
  finn_path <- glue::glue("{finn_sum_dir}/finngen_R12_{outcome_file}")

  finn_raw <- tryCatch(
    data.table::fread(finn_path),
    error = function(e) {
      warning(sprintf("Error reading FinnGen data for %s: %s", outcome_diag, e$message))
      NULL
    }
  )
  if (is.null(finn_raw)) {
    return(NULL)
  }

  finn_raw[, join_key := paste(`#chrom`, pos, alt, ref, sep = ":")]
  finn_with_rs <- finn_snplist[finn_raw, on = "join_key", nomatch = 0]
  finn_with_rs[, `:=`(Nca = Nca, Nco = Nco)]
  data.table::setkey(finn_with_rs, CHR, BP_hg37)

  list(finn_with_rs = finn_with_rs, Nca = Nca, Nco = Nco)
}

run_coloc_by_outcome_streamed <- function(mr_file,
                                          target_classes = c("metab", "IDP", "panukb", "pqtl"),
                                          coloc_dir,
                                          result_dir,
                                          result_prefix,
                                          pval_threshold = Inf,
                                          parallel_jobs = 40L,
                                          window_size = 500000L,
                                          finn_sum_dir = "/public/home2/nodecw_group/Finngen_GWAS_summary/R12",
                                          data_list_file,
                                          finn_snplist_file,
                                          diag_info_file,
                                          finngen_manifest_file,
                                          skip_existing = TRUE,
                                          preload_exposure_regions = TRUE,
                                          preload_scope = c("all", "outcome"),
                                          max_tasks_per_outcome = 200000L) {
  preload_scope <- match.arg(preload_scope)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  by_outcome_dir <- file.path(result_dir, paste0(result_prefix, "_by_outcome"))
  dir.create(by_outcome_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n[1] Reading inputs...\n")
  mr_result <- read_coloc_mr_result(mr_file)
  if (is.finite(pval_threshold) && "pval" %in% names(mr_result)) {
    mr_result <- mr_result |>
      dplyr::filter(.data$pval < pval_threshold)
  }
  mr_result_sig <- mr_result |>
    dplyr::filter(.data$file_class %in% target_classes)

  data_list <- qs::qread(data_list_file)
  finn_snplist <- data.table::fread(finn_snplist_file)
  finn_snplist <- finn_snplist[, .(
    join_key = paste(CHR, BP_hg38, A1, A2, sep = ":"),
    CHR,
    BP_hg37
  )]

  diag_info <- readxl::read_xlsx(diag_info_file) |>
    dplyr::mutate(gwas_summary = dplyr::coalesce(.data$phenocode, .data$phenocode2)) |>
    dplyr::filter(!is.na(.data$gwas_summary))
  diag_info_num <- data.table::fread(finngen_manifest_file) |>
    dplyr::select(phenocode, num_cases, num_controls)
  diag_info <- diag_info |>
    dplyr::inner_join(diag_info_num, by = c("gwas_summary" = "phenocode"))

  params_pairs_all <- mr_result_sig |>
    dplyr::select(outcome_diag, trait1_rawcol, trait1, file_class) |>
    dplyr::distinct() |>
    dplyr::rename(trait1_name = trait1)

  cat("  MR rows:", nrow(mr_result), "\n")
  cat("  Selected MR rows:", nrow(mr_result_sig), "\n")
  cat("  Outcome-exposure pairs:", nrow(params_pairs_all), "\n")
  cat("  Classes:\n")
  print(table(mr_result_sig$file_class, useNA = "ifany"))

  if (nrow(params_pairs_all) == 0) {
    stop("No outcome-exposure pairs found for target classes.")
  }

  diag_vec <- unique(params_pairs_all$outcome_diag)
  exposure_region_cache <- NULL
  if (preload_exposure_regions && preload_scope == "all") {
    exposure_region_cache <- preload_coloc_exposure_regions(
      params_pairs_all = params_pairs_all,
      coloc_dir = coloc_dir
    )
  }

  setup_parallel_with_progress(parallel_jobs)

  start_time <- Sys.time()
  all_coloc_results <- list()

  for (i in seq_along(diag_vec)) {
    outcome_diag <- diag_vec[i]
    outcome_outfile <- file.path(by_outcome_dir, sprintf("%s_%s.rds", result_prefix, outcome_diag))

    cat("###### ", i, "/", length(diag_vec), outcome_diag, " ######\n")
    cat("|| Begin at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Ncores:", parallel_jobs, "\n")

    if (skip_existing && file.exists(outcome_outfile)) {
      cat("|| Existing outcome result found, loading:", outcome_outfile, "\n")
      all_coloc_results[[outcome_diag]] <- readRDS(outcome_outfile)
      next
    }

    params_pairs_outcome <- params_pairs_all |>
      dplyr::filter(.data$outcome_diag == .env$outcome_diag)

    outcome_data <- read_finngen_outcome_for_coloc(
      outcome_diag = outcome_diag,
      diag_info = diag_info,
      finn_sum_dir = finn_sum_dir,
      finn_snplist = finn_snplist
    )
    if (is.null(outcome_data)) next

    exposure_region_cache_outcome <- exposure_region_cache
    if (preload_exposure_regions && preload_scope == "outcome") {
      cat("|| Preloading exposure region files for this outcome...\n")
      exposure_region_cache_outcome <- preload_coloc_exposure_regions(
        params_pairs_all = params_pairs_outcome,
        coloc_dir = coloc_dir
      )
    }

    cat("|| Building coloc task table for this outcome...\n")
    params_outcome <- build_coloc_params_for_outcome(
      params_pairs_outcome = params_pairs_outcome,
      data_list = data_list,
      coloc_dir = coloc_dir,
      exposure_region_cache = exposure_region_cache_outcome
    )

    if (nrow(params_outcome) == 0) {
      warning(sprintf("No coloc tasks built for outcome %s", outcome_diag))
      next
    }
    if (!is.null(max_tasks_per_outcome) && nrow(params_outcome) > max_tasks_per_outcome) {
      stop(sprintf(
        "Too many coloc tasks for %s: %d. Check outcome filtering before running.",
        outcome_diag,
        nrow(params_outcome)
      ))
    }

    finn_with_rs <- outcome_data$finn_with_rs
    Nca <- outcome_data$Nca
    Nco <- outcome_data$Nco

    res_list <- progressr::with_progress({
      p <- progressr::progressor(steps = nrow(params_outcome))

      future.apply::future_lapply(
        seq_len(nrow(params_outcome)),
        function(k) {
          row <- params_outcome[k, ]
          coloc_sub <- if (!is.null(exposure_region_cache_outcome)) {
            exposure_region_cache_outcome[[row$coloc_key]]
          } else {
            row$coloc_sub[[1]]
          }

          res <- run_one_coloc_task(
            trait1_rawcol = row$trait1_rawcol,
            trait1_name = row$trait1_name,
            file_class = row$file_class,
            lead_snp = row$lead_snp,
            outcome_diag = outcome_diag,
            Nca = Nca,
            Nco = Nco,
            finn_with_rs = finn_with_rs,
            coloc_sub = coloc_sub,
            data_list = data_list,
            window_size = window_size
          )

          p(sprintf("%s [%s] %s", row$trait1_name, row$file_class, row$lead_snp))

          if (!is.null(res)) {
            res$trait1_rawcol <- row$trait1_rawcol
            res$file_class <- row$file_class
          }
          res
        },
        future.seed = TRUE,
        future.globals = FALSE
      )
    })

    outcome_res <- dplyr::bind_rows(res_list)
    saveRDS(outcome_res, outcome_outfile)
    all_coloc_results[[outcome_diag]] <- outcome_res

    rm(outcome_data, finn_with_rs, params_outcome, res_list, outcome_res, exposure_region_cache_outcome)
    gc(verbose = FALSE)

    cat("|| End at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("###############################\n\n")
  }

  coloc_all <- dplyr::bind_rows(all_coloc_results)
  final_file <- file.path(result_dir, sprintf("%s_all.rds", result_prefix))
  saveRDS(coloc_all, final_file)

  cat("\n========== Coloc Complete ==========\n")
  cat("Output:", final_file, "\n")
  cat("Rows:", nrow(coloc_all), "\n")
  cat("Total time:", format(Sys.time() - start_time), "\n")

  invisible(coloc_all)
}
