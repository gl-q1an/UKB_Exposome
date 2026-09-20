############################################################
## Mediation analysis helper functions
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
    if (future::supportsMulticore()) {
      plan(multicore, workers = parallel_jobs)
    } else {
      plan(multisession, workers = parallel_jobs)
    }
  } else {
    handlers("progress")
    plan(sequential)
  }

  invisible(NULL)
}

############################################################
## Result template
############################################################

make_empty_mediation_result <- function(ext_use,
                                        int_use,
                                        ill_use,
                                        num_all = 0,
                                        num_case = 0,
                                        error_msg = NA_character_,
                                        warning_msg = NA_character_) {
  data.table::data.table(
    external_exposure = ext_use,
    internal_exposure = int_use,
    diseases = ill_use,
    num_all = num_all,
    num_case = num_case,
    CFI = NA_real_,
    RMSEA = NA_real_,
    SRMR = NA_real_,
    a_estimate = NA_real_,
    a_pvalues = NA_real_,
    b_estimate = NA_real_,
    b_pvalues = NA_real_,
    c_estimate = NA_real_,
    c_pvalues = NA_real_,
    ab_estimate = NA_real_,
    ab_pvalues = NA_real_,
    total_estimate = NA_real_,
    total_pvalues = NA_real_,
    error_msg = error_msg,
    warning_msg = warning_msg
  )
}

############################################################
## Covariate preparation
############################################################

normalize_binary_covariates <- function(data_use) {
  data_use |>
    dplyr::mutate(
      Sex = dplyr::case_when(
        Sex == "Male" ~ 1,
        Sex == "Female" ~ 0,
        TRUE ~ suppressWarnings(as.numeric(as.character(Sex)))
      )
    )
}

############################################################
## Single mediation model
############################################################

run_single_mediation <- function(df_subset,
                                 ext_use,
                                 int_use,
                                 ill_use,
                                 covariates) {
  if (length(covariates) != 2L) {
    stop("covariates must contain Age and Sex inputs, in that order.")
  }

  required_cols <- unique(c("eid", ext_use, int_use, "status", covariates))
  missing_cols <- setdiff(required_cols, names(df_subset))
  if (length(missing_cols) > 0) {
    return(make_empty_mediation_result(
      ext_use = ext_use,
      int_use = int_use,
      ill_use = ill_use,
      error_msg = paste("Missing columns:", paste(missing_cols, collapse = ", "))
    ))
  }

  df_subset <- data.table::as.data.table(df_subset)
  cols_to_extract <- c("eid", ext_use, int_use, "status", covariates)
  data_use <- df_subset[, ..cols_to_extract]
  data.table::setnames(
    data_use,
    old = cols_to_extract,
    new = c("eid", "ext", "int", "status", "Age", "Sex")
  )

  data_use <- normalize_binary_covariates(data_use)
  data_use <- data_use |>
    dplyr::mutate(
      status = suppressWarnings(as.numeric(as.character(status))),
      ext = suppressWarnings(as.numeric(as.character(ext))),
      int = suppressWarnings(as.numeric(as.character(int))),
      Age = suppressWarnings(as.numeric(as.character(Age)))
    ) |>
    tidyr::drop_na(ext, int, status, Age, Sex) |>
    dplyr::filter(status %in% c(0, 1))

  num_all <- nrow(data_use)
  num_case <- sum(data_use$status == 1)

  if (nrow(data_use) < 100 || sum(data_use$status == 1, na.rm = TRUE) < 50) {
    return(make_empty_mediation_result(
      ext_use = ext_use,
      int_use = int_use,
      ill_use = ill_use,
      num_all = nrow(data_use),
      num_case = sum(data_use$status == 1, na.rm = TRUE),
      error_msg = "Insufficient complete samples or cases"
    ))
  }


  model <- paste(
    "int ~ a*ext + Age + Sex",
    "status ~ b*int + cp*ext + Age + Sex",
    "ab := a * b",
    "total := cp + ab",
    sep = "\n"
  )

  tryCatch({
    warnings_collected <- character(0)

    sem_result <- withCallingHandlers(
      lavaan::sem(model, data = data_use, estimator = "MLR"),
      warning = function(w) {
        warnings_collected <<- c(warnings_collected, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    fit_measures <- lavaan::fitMeasures(sem_result, c("cfi", "rmsea", "srmr"))
    pe <- lavaan::parameterEstimates(sem_result)

    get_est <- function(lbl) {
      val <- pe$est[pe$label == lbl]
      if (length(val) == 0) NA_real_ else val[1]
    }
    get_pval <- function(lbl) {
      val <- pe$pvalue[pe$label == lbl]
      if (length(val) == 0) NA_real_ else val[1]
    }

    data.table::data.table(
      external_exposure = ext_use,
      internal_exposure = int_use,
      diseases = ill_use,
      num_all = num_all,
      num_case = num_case,
      CFI = fit_measures["cfi"],
      RMSEA = fit_measures["rmsea"],
      SRMR = fit_measures["srmr"],
      a_estimate = get_est("a"),
      a_pvalues = get_pval("a"),
      b_estimate = get_est("b"),
      b_pvalues = get_pval("b"),
      c_estimate = get_est("cp"),
      c_pvalues = get_pval("cp"),
      ab_estimate = get_est("ab"),
      ab_pvalues = get_pval("ab"),
      total_estimate = get_est("total"),
      total_pvalues = get_pval("total"),
      error_msg = NA_character_,
      warning_msg = if (length(warnings_collected) == 0) {
        NA_character_
      } else {
        paste(warnings_collected, collapse = " | ")
      }
    )
  }, error = function(e) {
    make_empty_mediation_result(
      ext_use = ext_use,
      int_use = int_use,
      ill_use = ill_use,
      num_all = num_all,
      num_case = num_case,
      error_msg = e$message
    )
  })
}

############################################################
## Output formatting
############################################################

coerce_mediation_numeric_cols <- function(df) {
  numeric_cols <- c(
    "num_all", "num_case",
    "CFI", "RMSEA", "SRMR",
    "a_estimate", "a_pvalues",
    "b_estimate", "b_pvalues",
    "c_estimate", "c_pvalues",
    "ab_estimate", "ab_pvalues",
    "total_estimate", "total_pvalues"
  )

  for (col in intersect(numeric_cols, names(df))) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }

  df
}

############################################################
## Diagnosis status cache
############################################################

load_diag_status_cache <- function(diag_list, diseases) {
  cache <- vector("list", length(diseases))
  names(cache) <- diseases

  for (disease in diseases) {
    if (!disease %in% names(diag_list)) {
      cache[[disease]] <- NULL
      next
    }

    cache[[disease]] <- diag_list[[disease]] |>
      dplyr::select(eid, status)
  }

  cache
}

############################################################
## Batch mediation runner
############################################################

run_mediation_params <- function(params,
                                 data_all,
                                 diag_status_cache,
                                 covariates,
                                 progress_label = "mediation") {
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(params))

    future.apply::future_lapply(seq_len(nrow(params)), function(idx) {
      row <- params[idx, ]

      diag_df <- diag_status_cache[[row$diseases]]
      if (is.null(diag_df)) {
        p(sprintf("%s | missing disease", progress_label))
        return(make_empty_mediation_result(
          ext_use = row$external_exposure,
          int_use = row$internal_exposure,
          ill_use = row$diseases,
          error_msg = "Disease not found in diagnosis list"
        ))
      }

      cols_needed <- unique(c(
        "eid",
        row$external_exposure,
        row$internal_exposure,
        covariates
      ))

      df_subset <- data_all |>
        dplyr::select(dplyr::any_of(cols_needed)) |>
        dplyr::inner_join(diag_df, by = "eid")

      res <- run_single_mediation(
        df_subset = df_subset,
        ext_use = row$external_exposure,
        int_use = row$internal_exposure,
        ill_use = row$diseases,
        covariates = covariates
      )

      p(sprintf(
        "%s | %s -> %s -> %s",
        progress_label,
        row$external_exposure,
        row$internal_exposure,
        row$diseases
      ))
      res
    }, future.seed = TRUE, future.scheduling = 1)
  })
}

############################################################
## Save results
############################################################

save_mediation_result_rds <- function(result_dir, prefix, result) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(result_dir, paste0(prefix, ".rds"))
  saveRDS(result, out_file)
  out_file
}
