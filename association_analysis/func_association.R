####################################################
## Shared functions for the Cox analyses
####################################################

setup_parallel_with_progress <- function(parallel_jobs = 1L,
                                         max_global_size_gb = 50) {
  stopifnot(length(parallel_jobs) == 1L, is.finite(parallel_jobs),
            parallel_jobs >= 1L, parallel_jobs == as.integer(parallel_jobs))
  options(future.globals.maxSize = max_global_size_gb * 1024^3)
  progressr::handlers(progressr::handler_progress(
    format = "[-] :current/:total (:percent) | Elapsed: :elapsedfull | ETA: :eta",
    clear = FALSE
  ))
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  if (parallel_jobs == 1L) {
    future::plan(future::sequential)
  } else if (future::supportsMulticore()) {
    future::plan(future::multicore, workers = parallel_jobs)
  } else {
    future::plan(future::multisession, workers = parallel_jobs)
  }
}


####################################################
## Exposure-specific covariates
####################################################

build_exposure_params <- function(data, exp_info, protein = FALSE) {
  info <- as.data.frame(exp_info)
  if (!"clean_subcat" %in% names(info)) {
    info$clean_subcat <- if ("exp_subcate" %in% names(info)) {
      info$exp_subcate
    } else info$Level2
  }
  if (!"clean_subcat" %in% names(info)) stop("Missing exposure category")
  info <- info |>
    dplyr::filter(!is.na(exp_code), exp_code %in% names(data),
                  !is.na(Time_variableid))
  if (protein) {
    info <- info |> dplyr::filter(Time_variableid == "p53_i0")
    if ("Level2" %in% names(info)) {
      info <- info |> dplyr::filter(Level2 == "Proteomics")
    }
  } else {
    info <- info |> dplyr::filter(!is.na(coding_note))
  }
  info <- info |> dplyr::distinct(exp_code, .keep_all = TRUE)
  if (nrow(info) == 0L) stop("No eligible exposures")

  dplyr::bind_rows(lapply(seq_len(nrow(info)), function(i) {
    exposure <- info$exp_code[[i]]
    time_col <- info$Time_variableid[[i]]
    category <- info$clean_subcat[[i]]
    covars <- c("Sex", paste0("Age_", time_col))
    if (time_col == "p53_i0") covars <- c(covars, "AccCeni0")
    if (time_col == "p53_i2") covars <- c(covars, "AccCeni2")
    if (protein) {
      covars <- c(covars, "FastingTime", paste0(exposure, "_SampAge"))
    } else if (identical(category, "Blood biochemistry")) {
      covars <- c(covars, "FastingTime", "SamAge_BloodBioch")
    } else if (identical(category, "Blood count")) {
      covars <- c(covars, "FastingTime")
    } else if (identical(category, "Urine biochemistry")) {
      covars <- c(covars, "SamAge_UrineBioch")
    } else if (identical(category, "Blood metabolomics")) {
      covars <- c(covars, "FastingTime", "SamAge_BloodMetab", "MetabBatch")
    } else if (identical(category, "Brain MRI")) {
      covars <- c(covars, "TIV")
    }
    tibble::tibble(sing_colname = exposure, time_col = time_col,
                   covar_run = list(setdiff(unique(covars), exposure)))
  }))
}

get_cox_diag_names <- function(diag_list, disease_info) {
  out <- unique(c(intersect(disease_info$disease_code, names(diag_list)),
                  grep("^ZZZDEATH_", names(diag_list), value = TRUE)))
  if (!length(out)) stop("No eligible disease or mortality outcomes")
  out
}

expand_cox_pairs <- function(params, diag_names) {
  pairs <- as.data.frame(params)[
    rep(seq_len(nrow(params)), each = length(diag_names)), , drop = FALSE]
  pairs$sing_diag <- rep(diag_names, times = nrow(params))
  pairs
}

####################################################
## Cox
####################################################

run_single_cox <- function(df_all,
                           sing_colname,
                           sing_diag,
                           time_col,
                           diag_list,
                           covar_vec,
                           sex_vec) {
  case_date_col <- paste0("diag_time_", time_col)
  cox_result <- list(
    trait1 = sing_colname,
    trait2 = sing_diag,
    ncase = 0,
    ncontrol = 0,
    cox_hr = NA_real_,
    cox_se = NA_real_,
    cox_pval = NA_real_,
    formula = "",
    success = "n < 100"
  )

  df_diag <- dplyr::as_tibble(diag_list[[sing_diag]])
  required_diag_cols <- c("eid", "status", case_date_col)

  status_text <- as.character(df_diag$status)

  df_diag$status <- as.integer(status_text)

  df_run <- df_diag |>
    dplyr::select(dplyr::all_of(required_diag_cols)) |>
    dplyr::inner_join(df_all, by = "eid") |>
    tidyr::drop_na(dplyr::all_of(sing_colname))

  covar_vec_cox <- covar_vec
  if (sing_diag %in% sex_vec) {
    covar_vec_cox <- setdiff(covar_vec_cox, "Sex")
  }
  complete_case_cols <- unique(c(
    "status", case_date_col, sing_colname, covar_vec_cox
  ))

  df_cox <- df_run |>
    dplyr::filter(
      status %in% c(0L, 1L),
      is.finite(.data[[case_date_col]]),
      .data[[case_date_col]] > 0
    ) |>
    tidyr::drop_na(dplyr::all_of(complete_case_cols))
  numeric_cols <- complete_case_cols[vapply(df_cox[complete_case_cols], is.numeric, logical(1))]
  df_cox <- df_cox |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(numeric_cols), is.finite))

  n_control <- sum(df_cox$status == 0L)
  n_new_case <- sum(
    df_cox$status == 1L & df_cox[[case_date_col]] > 0
  )
  cox_result[c("ncase", "ncontrol")] <- c(n_new_case, n_control)

  if (n_new_case < 100L) {
    return(cox_result)
  }

  if (length(unique(df_cox[[sing_colname]])) < 2L) {
    cox_result$success <- "Exposure has < 2 values"
    return(cox_result)
  }

  covar_vec_cox <- covar_vec_cox[vapply(covar_vec_cox, function(x) {
    length(unique(df_cox[[x]][!is.na(df_cox[[x]])])) > 1L
  }, logical(1))]

  cox_formula_text <- paste0(
    "survival::Surv(", case_date_col, ", status) ~ ",
    paste(paste0("`", c(sing_colname, covar_vec_cox), "`"), collapse = " + ")
  )
  cox_result["formula"] <- cox_formula_text
  cox_fit <- survival::coxph(
    stats::as.formula(cox_formula_text),
    data = df_cox,
    id = eid
  )

  cox_coef <- summary(cox_fit)$coefficients
  rownames(cox_coef) <- gsub("`", "", rownames(cox_coef), fixed = TRUE)
  if (!sing_colname %in% rownames(cox_coef)) {
    cox_result[c("cox_hr", "cox_se", "cox_pval")] <- list(NA_real_, NA_real_, NA_real_)
    cox_result["success"] <- "Error: exposure coefficient is absent"
    return(cox_result)
  }

  estimates <- cox_coef[sing_colname, c("coef", "se(coef)", "Pr(>|z|)")]
  if (any(!is.finite(estimates)) || estimates[[2]] <= 0 ||
      !is.finite(exp(estimates[[1]])) || exp(estimates[[1]]) <= 0) {
    cox_result$success <- "Error: exposure coefficient is not estimable"
    return(cox_result)
  }

  cox_result["cox_hr"] <- exp(cox_coef[sing_colname, "coef"])
  cox_result["cox_se"] <- cox_coef[sing_colname, "se(coef)"]
  cox_result["cox_pval"] <- cox_coef[sing_colname, "Pr(>|z|)"]
  cox_result["success"] <- "1"

  cox_result
}

get_subgroup_defs <- function() {
  list(
    list(type = "sex", level = "Female", label = "sex_Female"),
    list(type = "sex", level = "Male", label = "sex_Male"),
    list(type = "age", level = "lt60", label = "age_lt60"),
    list(type = "age", level = "ge60", label = "age_ge60")
  )
}

subgroup_filter_col <- function(subgroup, time_col) {
  if (identical(subgroup$type, "sex")) {
    return("Sex")
  }
  if (identical(subgroup$type, "age")) {
    return(paste0("Age_", time_col))
  }
  stop("Unsupported subgroup type: ", subgroup$type)
}

adjust_params_for_subgroup <- function(params, subgroup) {
  params <- as.data.frame(params)
  params$covar_run <- Map(function(covar_run, time_col) {
    removed <- subgroup_filter_col(subgroup, time_col)
    setdiff(unique(covar_run), removed)
  }, params$covar_run, params$time_col)
  params
}

filter_for_subgroup <- function(df, subgroup, time_col) {
  filter_col <- subgroup_filter_col(subgroup, time_col)
  if (!filter_col %in% names(df)) {
    stop("Missing subgroup filter column: ", filter_col)
  }

  filter_value <- df[[filter_col]]
  keep <- !is.na(filter_value)
  if (identical(subgroup$type, "sex")) {
    keep <- keep & as.character(filter_value) == subgroup$level
  } else if (identical(subgroup$type, "age") && identical(subgroup$level, "lt60")) {
    keep <- keep & filter_value < 60
  } else if (identical(subgroup$type, "age") && identical(subgroup$level, "ge60")) {
    keep <- keep & filter_value >= 60
  } else {
    stop("Unsupported subgroup: ", subgroup$type, "/", subgroup$level)
  }

  df[keep, , drop = FALSE]
}

select_subgroup_model_data <- function(df, sing_colname, covar_vec, subgroup, time_col) {
  filter_col <- subgroup_filter_col(subgroup, time_col)
  cols_needed <- unique(c("eid", sing_colname, covar_vec, filter_col))
  missing_cols <- setdiff(cols_needed, names(df))
  if (length(missing_cols) > 0L) {
    stop("Missing model columns: ", paste(missing_cols, collapse = ", "))
  }

  selected <- dplyr::select(df, dplyr::all_of(cols_needed))
  filter_for_subgroup(selected, subgroup, time_col)
}

####################################################
## Parallel exposure-outcome analysis, without batch scheduling
####################################################

run_cox_pairs <- function(data, params, diag_list, sex_vec, subgroup = NULL) {
  required <- unique(c("eid", params$sing_colname,
                       unlist(params$covar_run, use.names = FALSE)))
  if (!is.null(subgroup)) {
    required <- unique(c(required, vapply(params$time_col, function(time_col) {
      subgroup_filter_col(subgroup, time_col)
    }, character(1))))
  }
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) stop("Missing model columns: ", paste(missing_cols, collapse = ", "))

  results <- progressr::with_progress({
    p <- progressr::progressor(nrow(params))
    future.apply::future_lapply(seq_len(nrow(params)), function(i) {
      row <- params[i, , drop = FALSE]
      exposure <- row$sing_colname[[1]]
      covars <- row$covar_run[[1]]
      time_col <- row$time_col[[1]]
      if (is.null(subgroup)) {
        model_data <- dplyr::select(data, dplyr::all_of(unique(c("eid", exposure, covars))))
      } else {
        model_data <- select_subgroup_model_data(data, exposure, covars, subgroup, time_col)
      }
      result <- run_single_cox(model_data, exposure, row$sing_diag[[1]],
                               time_col, diag_list, covars, sex_vec)
      result$time_col <- time_col
      p()
      tibble::as_tibble(result)
    }, future.seed = TRUE, future.scheduling = 1)
  })
  dplyr::bind_rows(results)
}