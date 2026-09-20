############################################################
## Shared functions: ERS, PRS, Cox interaction and weighted PAF
############################################################

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

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  values <- sort(table(x), decreasing = TRUE)
  as.numeric(names(values)[1L])
}


convert_ers_column <- function(x, colname) {
  x_char <- as.character(x)
  x_char[x_char == ""] <- NA_character_
  unique_values <- unique(x_char[!is.na(x_char)])

  if (length(unique_values) == 0L) {
    return(rep(NA_real_, length(x_char)))
  }
  if (all(grepl("^-?[0-9]+(\\.[0-9]+)?$", unique_values))) {
    return(as.numeric(x_char))
  }

  second_character <- substr(x_char, 2L, 2L)
  bad_values <- unique_values[!grepl("^.[0-9]", unique_values)]
  if (length(bad_values) > 0L) {
    stop(
      "Cannot convert ", colname, "; unexpected values: ",
      paste(head(bad_values, 10L), collapse = ", ")
    )
  }
  as.numeric(second_character)
}


convert_ers_exposure <- function(x, coding_note, colname) {
  if (identical(coding_note, "Continuous")) {
    return(suppressWarnings(as.numeric(x)))
  }
  if (coding_note %in% c("Binary", "Ordinal categorical")) {
    return(convert_ers_column(x, colname))
  }
  stop("Unsupported coding_note for ", colname, ": ", coding_note)
}


ers_impute_value <- function(x, coding_note) {
  if (identical(coding_note, "Continuous")) {
    return(stats::median(x, na.rm = TRUE))
  }
  get_mode(x)
}


impute_ers_with_train <- function(train_x,
                                  test_x,
                                  train_sex,
                                  test_sex,
                                  coding_note) {
  global_value <- ers_impute_value(train_x, coding_note)
  if (!is.finite(global_value)) {
    return(NULL)
  }

  sex_levels <- sort(unique(train_sex[!is.na(train_sex)]))
  by_sex <- setNames(rep(global_value, length(sex_levels)), as.character(sex_levels))
  for (sex_level in sex_levels) {
    value <- ers_impute_value(
      train_x[!is.na(train_sex) & train_sex == sex_level],
      coding_note
    )
    if (is.finite(value)) {
      by_sex[[as.character(sex_level)]] <- value
    }
  }

  fill_missing <- function(x, sex) {
    for (sex_level in names(by_sex)) {
      index <- is.na(x) & !is.na(sex) & as.character(sex) == sex_level
      x[index] <- by_sex[[sex_level]]
    }
    x[is.na(x)] <- global_value
    x
  }

  list(
    train = fill_missing(train_x, train_sex),
    test = fill_missing(test_x, test_sex),
    global_value = global_value,
    by_sex = by_sex
  )
}


remove_ers_collinear_vars <- function(data, id_col = "eid", threshold = 0.8) {
  exposure_data <- dplyr::select(data, -dplyr::all_of(id_col))
  if (ncol(exposure_data) <= 1L) {
    return(names(exposure_data))
  }

  correlation <- suppressWarnings(stats::cor(
    exposure_data,
    use = "pairwise.complete.obs"
  ))
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  absolute_correlation <- abs(correlation)
  scores <- rowSums(absolute_correlation, na.rm = TRUE) - 1
  variables_to_keep <- colnames(correlation)

  repeat {
    current_correlation <- absolute_correlation[
      variables_to_keep,
      variables_to_keep,
      drop = FALSE
    ]
    high_pairs <- which(
      current_correlation >= threshold & upper.tri(current_correlation),
      arr.ind = TRUE
    )
    if (nrow(high_pairs) == 0L) {
      break
    }

    variable_1 <- rownames(current_correlation)[high_pairs[1L, 1L]]
    variable_2 <- rownames(current_correlation)[high_pairs[1L, 2L]]
    variable_to_remove <- if (scores[[variable_1]] >= scores[[variable_2]]) {
      variable_1
    } else {
      variable_2
    }
    variables_to_keep <- setdiff(variables_to_keep, variable_to_remove)
  }

  variables_to_keep
}


############################################################
## Exposure groups used in the manuscript
############################################################

make_ers_groups <- function(exp_info, lifestyle_info) {
  external <- c("Health history", "Lifestyle factors",
                "Natural/occupational environment", "Socioeconomic factors")
  internal <- c("Clinical biomarkers", "Metabolome", "Physical measures")
  subtypes <- c("Smoking_alcohol", "Diet", "Sleep_others", "Physical activity")
  groups <- tibble::tibble(
    group_id = c("health_history", "lifestyle", "environment", "socioeconomic",
                 "in", "out", "smoking_alcohol", "diet", "sleep_others", "physical_activity"),
    category = c(external, "Intrinsic ERS", "Extrinsic ERS", subtypes),
    group_set = c(rep("external_category", 4), rep("inout", 2), rep("lifestyle_subtype", 4)),
    modifiable = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, rep(TRUE, 4))
  )
  exposures <- lapply(external, function(category) {
    exp_info$exp_code[exp_info$exp_cate == category]
  })
  exposures <- c(exposures, list(
    exp_info$exp_code[exp_info$exp_cate %in% internal],
    exp_info$exp_code[exp_info$exp_cate %in% external]
  ), lapply(subtypes, function(category) {
    lifestyle_info$exp_code[lifestyle_info$exp_cate == "Lifestyle factors" &
                             lifestyle_info$new_exp_cate == category &
                             lifestyle_info$Time_variableid == "p53_i0"]
  }))
  groups$exposures <- lapply(exposures, function(x) {
    intersect(unique(stats::na.omit(x)), exp_info$exp_code)
  })
  groups
}

score_name <- function(diagnosis) {
  paste0("score_", gsub("[^A-Za-z0-9_.-]+", "_", diagnosis))
}

parallel_map <- function(x, fun) {
  progressr::with_progress({
    p <- progressr::progressor(along = x)
    future.apply::future_lapply(x, function(value) {
      result <- fun(value)
      p()
      result
    }, future.seed = TRUE, future.scheduling = 1)
  })
}

model_covariates <- function(data, covars, diagnosis, sex_vec) {
  if (diagnosis %in% sex_vec) covars <- setdiff(covars, "Sex")
  covars[vapply(covars, function(x) {
    dplyr::n_distinct(data[[x]], na.rm = TRUE) > 1L
  }, logical(1))]
}

incident_data <- function(data, diagnosis, diag_list, time_col = "p53_i0") {

  time_name <- paste0("diag_time_", time_col)
  diagnosis_data <- as.data.frame(diag_list[[diagnosis]])
  diagnosis_data$status <- as.integer(as.character(diagnosis_data$status))
  diagnosis_data |>
    dplyr::select(eid, status, dplyr::all_of(time_name)) |>
    dplyr::inner_join(data, by = "eid") |>
    dplyr::filter(status %in% c(0L, 1L),
                  is.finite(.data[[time_name]]), .data[[time_name]] > 0)
}

cox_formula <- function(terms, time_col = "p53_i0") {
  stats::as.formula(paste0("survival::Surv(diag_time_", time_col,
                          ", status) ~ ", paste(terms, collapse = " + ")))
}

quote_name <- function(x) paste0("`", x, "`")

############################################################
## Train disease-specific ERS with backward Cox selection
############################################################

fit_ers_beta <- function(data, exposures, diagnosis, diag_list, covars, sex_vec,
                         p_threshold = 0.05) {
  result <- list(diag = diagnosis, ncase = 0L, ncontrol = 0L,
                 beta = setNames(numeric(), character()), success = "no exposures")
  if (!length(exposures)) return(result)
  data <- data |> dplyr::select(eid, dplyr::all_of(c(exposures, covars)))
  data <- incident_data(data, diagnosis, diag_list)
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  data <- data |>
    tidyr::drop_na(dplyr::all_of(c("status", "diag_time_p53_i0", exposures, covars)))
  result$ncase <- sum(data$status == 1L)
  result$ncontrol <- sum(data$status == 0L)
  result$success <- "fewer than 50 cases or controls"
  if (result$ncase < 50L || result$ncontrol < 50L) return(result)
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  current <- exposures[vapply(exposures, function(x) dplyr::n_distinct(data[[x]]) > 1L, logical(1))]

  while (length(current)) {
    model <- survival::coxph(cox_formula(quote_name(c(current, covars))), data = data, id = eid)
    coefficients <- summary(model)$coefficients
    rownames(coefficients) <- gsub("`", "", rownames(coefficients), fixed = TRUE)
    pvalues <- coefficients[current, "Pr(>|z|)"]
    names(pvalues) <- current
    if (all(is.finite(pvalues)) && all(pvalues < p_threshold)) {
      result$beta <- setNames(coefficients[current, "coef"], current)
      result$success <- "1"
      return(result)
    }
    pvalues[!is.finite(pvalues)] <- Inf
    current <- setdiff(current, names(which.max(pvalues)))
  }
  result$success <- "no exposure retained"
  result
}

calculate_ers <- function(test_data, beta_results) {
  out <- data.frame(eid = test_data$eid)
  for (result in beta_results) {
    column <- score_name(result$diag)
    if (result$success == "1") {
      x <- as.matrix(test_data[, names(result$beta), drop = FALSE])
      out[[column]] <- as.numeric(x %*% result$beta)
    } else {
      out[[column]] <- NA_real_
    }
  }
  out
}

############################################################
## PRS threshold selection in training data only
############################################################

select_prs_threshold <- function(diagnosis, prs_dir, train_eids, test_eids,
                                 diag_list, covar_data, covars, sex_vec, min_snp = 20L) {
  prefix <- file.path(prs_dir, paste0("prs_", diagnosis), paste0("prs_", diagnosis, "_detail"))
  scores <- as.data.frame(data.table::fread(paste0(prefix, ".all_score")))
  info <- data.table::fread(paste0(prefix, ".prsice"), colClasses = list(character = "Threshold"))
  thresholds <- setdiff(names(scores), c("FID", "IID"))
  snp_map <- setNames(as.integer(info$Num_SNP), paste0("Pt_", info$Threshold))

  data <- as.data.frame(diag_list[[diagnosis]]) |>
    dplyr::mutate(status = as.integer(as.character(status))) |>
    dplyr::filter(eid %in% train_eids, status %in% c(0L, 1L)) |>
    dplyr::select(eid, status) |>
    dplyr::inner_join(scores, by = c("eid" = "IID")) |>
    dplyr::inner_join(covar_data, by = "eid")
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  fits <- lapply(thresholds, function(threshold) {
    model_data <- data |> tidyr::drop_na(dplyr::all_of(c("status", threshold, covars)))
    result <- tibble::tibble(threshold = threshold, num_snp = unname(snp_map[threshold]),
                             z = NA_real_, pval = NA_real_)
    if (dplyr::n_distinct(model_data[[threshold]]) < 2L ||
        dplyr::n_distinct(model_data$status) < 2L) return(result)

    model <- stats::glm(stats::reformulate(quote_name(c(threshold, covars)), "status"),
                         data = model_data, family = stats::binomial())
    coefficients <- summary(model)$coefficients
    rownames(coefficients) <- gsub("`", "", rownames(coefficients), fixed = TRUE)
    result$z <- coefficients[threshold, "z value"]
    result$pval <- coefficients[threshold, "Pr(>|z|)"]
    result
  }) |> dplyr::bind_rows()
  candidates <- fits |> dplyr::filter(is.finite(z), num_snp >= min_snp)
  fallback <- nrow(candidates) == 0L
  if (fallback) candidates <- fits |> dplyr::filter(is.finite(z))
  if (nrow(candidates) == 0L) return(NULL)
  # Preserve the source rule: largest signed z; no absolute-value selection.
  best <- candidates |> dplyr::slice_max(z, n = 1L, with_ties = FALSE)
  best$diag <- diagnosis
  best$fewer_than_min_snp <- fallback
  list(selection = best,
       score = data.frame(eid = test_eids,
                          value = scores[[best$threshold]][match(test_eids, scores$IID)]))
}

############################################################
## Continuous ERS/PRS Cox models and multiplicative interaction
############################################################

standardize_scores <- function(data) {
  columns <- setdiff(names(data), "eid")
  for (column in columns) {
    value <- as.numeric(data[[column]])
    value[!is.finite(value)] <- NA_real_
    sigma <- stats::sd(value, na.rm = TRUE)
    data[[column]] <- if (is.finite(sigma) && sigma > 0) {
      (value - mean(value, na.rm = TRUE)) / sigma
    } else rep(NA_real_, length(value))
  }
  data
}

fit_score_cox <- function(data, diagnosis, diag_list, covars, sex_vec,
                          interaction = FALSE) {
  predictors <- if (interaction) c("ers", "prs") else "ers"
  data <- data |> dplyr::select(eid, dplyr::all_of(c(predictors, covars)))
  data <- incident_data(data, diagnosis, diag_list)
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  data <- data |> tidyr::drop_na(dplyr::all_of(c(predictors, covars)))
  ncase <- sum(data$status == 1L)
  ncontrol <- sum(data$status == 0L)
  if (ncase < 100L || ncontrol < 100L ||
      any(vapply(data[predictors], dplyr::n_distinct, integer(1)) < 2L)) return(NULL)
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  terms <- c(if (interaction) "ers * prs" else "ers", quote_name(covars))
  fit <- survival::coxph(cox_formula(terms), data = data, id = eid)
  coefficients <- summary(fit)$coefficients
  keep <- if (interaction) c("ers", "prs", "ers:prs") else "ers"
  tibble::tibble(diag = diagnosis, term = keep,
                 hr = exp(coefficients[keep, "coef"]),
                 se = coefficients[keep, "se(coef)"],
                 pval = coefficients[keep, "Pr(>|z|)"],
                 ncase = ncase, ncontrol = ncontrol,
                 formula = paste(deparse(stats::formula(fit)), collapse = " "))
}

############################################################
## Component HR/AF: intermediate quantities for weighted PAF
############################################################

component_af <- function(data, diagnosis, diag_list, covars, sex_vec,
                         top_fraction = 1 / 3) {
  result <- tibble::tibble(diag = diagnosis, source_hr = NA_real_,
                           source_hr_pval = NA_real_, source_paf = NA_real_, AF = 0,
                           high_threshold = NA_real_, exposed_proportion = NA_real_,
                           source_paf_n = 0L, n_model = 0L, success = "insufficient data")
  data <- data |> dplyr::select(eid, score, dplyr::all_of(covars))
  data <- incident_data(data, diagnosis, diag_list) |>
    dplyr::filter(is.finite(score))
  if (nrow(data) == 0L || dplyr::n_distinct(data$score) < 2L) return(result)
  threshold <- as.numeric(stats::quantile(data$score, 1 - top_fraction))
  data$high <- as.integer(data$score >= threshold)

  n_high <- sum(data$high == 1L)
  n_low <- sum(data$high == 0L)
  prevalence <- n_high / (n_high + n_low)
  if (n_high < 50L || n_low < 50L) return(result)
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  data <- data |> tidyr::drop_na(dplyr::all_of(covars))
  if (sum(data$status == 1L) < 50L || sum(data$status == 0L) < 50L ||
      dplyr::n_distinct(data$high) < 2L) return(result)
  covars <- model_covariates(data, covars, diagnosis, sex_vec)
  fit <- survival::coxph(cox_formula(c("high", quote_name(covars))), data = data, id = eid)
  coefficients <- summary(fit)$coefficients["high", ]
  hr <- exp(coefficients[["coef"]])
  pval <- coefficients[["Pr(>|z|)"]]

  source_af <- prevalence * (hr - 1) / (1 + prevalence * (hr - 1))
  result$source_hr <- hr
  result$source_hr_pval <- pval
  result$source_paf <- source_af
  result$AF <- if (is.finite(hr) && hr > 1 && is.finite(pval) &&
                   pval >= 0 && pval < 0.05 && is.finite(source_af)) source_af else 0
  result$high_threshold <- threshold
  result$exposed_proportion <- prevalence
  result$source_paf_n <- n_high + n_low
  result$n_model <- nrow(data)
  result$success <- if (is.finite(hr) && is.finite(pval)) "1" else "not estimable"
  result

}

############################################################
## Weighted PAF using the original weighting convention
############################################################

weighted_paf <- function(scores, component_results, top_fraction = 1 / 3) {

  categories <- setdiff(names(scores), "eid")
  binary <- list()
  thresholds <- numeric()
  for (category in categories) {
    value <- scores[[category]]
    value[!is.finite(value)] <- NA_real_
    if (sum(!is.na(value)) < 100L || dplyr::n_distinct(value, na.rm = TRUE) < 2L) next
    threshold <- as.numeric(stats::quantile(value, 1 - top_fraction, na.rm = TRUE))
    high <- as.integer(value >= threshold)
    if (stats::var(high, na.rm = TRUE) < 1e-10) next
    binary[[category]] <- high
    thresholds[[category]] <- threshold
  }
  if (length(binary) < 2L) return(NULL)
  binary <- stats::na.omit(as.data.frame(binary, check.names = FALSE))
  if (nrow(binary) < 100L) return(NULL)
  correlation <- psych::tetrachoric(binary, na.rm = TRUE)$rho
  eigen_result <- eigen(correlation, symmetric = TRUE)
  retained <- which(eigen_result$values > 1)
  if (!length(retained)) retained <- 1L

  communality <- rowSums(eigen_result$vectors[, retained, drop = FALSE]^2)
  names(communality) <- colnames(binary)
  result <- component_results |>
    dplyr::filter(category_key %in% names(communality)) |>
    dplyr::mutate(communality = unname(communality[category_key]),
                  weight_threshold = unname(thresholds[category_key]))
  if (nrow(result) < 2L) return(NULL)
  overall <- 1 - prod(1 - (1 - result$communality) * result$AF)
  result$weighted_PAF <- if (sum(result$AF) > 0) result$AF / sum(result$AF) * overall else 0
  result$overall_PAF <- overall
  result$n_samples <- nrow(binary)
  result$n_input_categories <- length(categories)
  result$n_valid_categories <- nrow(result)
  result$n_components <- length(retained)
  result
}
