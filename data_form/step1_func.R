library(glue)

################################################################
## Sort matched column names in numeric order.
################################################################
sort_ukb_cols <- function(x) {
  stopifnot(is.character(x))

  library(stringr)
  library(dplyr)

  tibble(col = x) %>%
    mutate(
      p = as.integer(str_extract(col, "(?<=^p)\\d+")),
      i = as.integer(str_extract(col, "(?<=_i)\\d+")),
      a = as.integer(str_extract(col, "(?<=_a)\\d+"))
    ) %>%
    arrange(p, i, a) %>%
    pull(col)
}

################################################################
## This function processes work-related data.
## Group occupations into broad categories.
################################################################
get_top_level_mapping <- function(df) {
  # Assume columns are coding, meaning, node_id, and parent_id.
  # Modify this if column names differ.

  # Create a node_id-to-meaning mapping.
  node_to_meaning <- setNames(df$meaning, df$node_id)

  # Create a node_id-to-parent_id mapping.
  node_to_parent <- setNames(df$parent_id, df$node_id)

  # Create a mapping from top-level categories to abbreviations.
  category_mapping <- c(
    "agriculture, horticulture, fishing, other work with animals (including managers)" = "agriculture",
    "armed forces, emergency services, security, health & safety (including managers)" = "armed_forces_emergency",
    "cleaning, caretaking, waste collection, pest control (including managers)" = "cleaning",
    "construction, building, demolition or maintenance (including managers)" = "construction",
    "education, school-related work (including managers)" = "education",
    "health (human or animal), residential/social/religious care, undertaking (including managers)" = "health",
    "mining, quarrying, energy production, water treatment (including managers)" = "mining",
    "office-based work: professional, managerial, administrative or general office/clerical" = "office_work",
    "personal services, travel/tourism, hospitality (including managers)" = "tourism",
    "routine factory-based manufacturing (including managers)" = "manufacturing",
    "science, research, engineering, computer technology (including managers)" = "science",
    "selling and shop work (retail/wholesale), storage and distribution (including managers)" = "sales",
    "skilled manual work (including managers)" = "skilled_trades",
    "sport, culture, arts, media, entertainment (including managers)" = "sport_media",
    "transport (road, rail, air, water), work with other mobile machinery (including managers)" = "transport"
  )

  # Recursively find the top-level node.
  find_top_level <- function(node_id) {
    parent_id <- node_to_parent[as.character(node_id)]

    # If parent_id is 0, the current node is top-level.
    if (!is.na(parent_id) && parent_id == 0) {
      return(node_id)
    }

    # If parent_id does not exist, return the current node.
    if (is.na(parent_id) || !as.character(parent_id) %in% names(node_to_parent)) {
      return(node_id)
    }

    # Recursively find the top-level parent node.
    return(find_top_level(parent_id))
  }

  # Find the top-level node for each node.
  result <- df %>%
    rowwise() %>%
    mutate(
      top_node_id = find_top_level(node_id),
      top_level_meaning = node_to_meaning[as.character(top_node_id)],
      category_short = category_mapping[top_level_meaning]
    ) %>%
    ungroup() %>%
    mutate(unique_sig = paste(coding,parent_id,sep='_'),
           check_col = paste(meaning,category_short,sep='_')) |>
    select(-coding,-node_id,-parent_id,-selectable,-top_node_id,-top_level_meaning,)

  return(result)
}

################################################################
## This function log-transforms data and replaces zeros with half of the smallest positive value.
## Used for blood-sample and metabolomics variables.
################################################################
lognmin <- function(x) {
  # Handle non-numeric input.
  if (!is.numeric(x)) {
    original_na <- is.na(x)
    x_numeric <- suppressWarnings(as.numeric(x))
    conversion_na <- is.na(x_numeric)
    new_na_count <- sum(conversion_na & !original_na)

    if (new_na_count > 0) {
      cat(paste("Non-numeric input produced", new_na_count, "NA values during numeric conversion\n"))
    } else {
      cat("Input was converted from non-numeric to numeric type\n")
    }
    x <- x_numeric
  }

  # Handle values below zero.
  if (any(x < 0, na.rm = TRUE)) {
    warning("Input contains values below zero; these values will be replaced with NA\n")
    x[x < 0] <- NA
  }

  # Return directly if all values are NA.
  if (all(is.na(x))) {
    cat("done!\n")
    return(x)
  }

  # Find the smallest positive value.
  min_positive <- min(x[x > 0], na.rm = TRUE)

  # Replace zeros with half of the smallest positive value.
  x[x == 0 & !is.na(x)] <- min_positive / 2

  # Take the log.
  result <- log(x)

  cat("done!\n")
  return(result)
}

################################################################
## General processing function.
################################################################

process_work_data_vectorized <- function(data,
                                         cutoff_year_col,
                                         mode = c("work_category", "gap", "work_exposure"),
                                         job_mapping = NULL,
                                         work_categories = NULL,
                                         exp_list = NULL) {

  mode <- match.arg(mode)

  # Create a copy.
  data <- copy(data)

  # Process cutoff_year.
  data[, cutoff_year_numeric := {
    cy <- get(cutoff_year_col)
    if (is.character(cy) || inherits(cy, "Date")) {
      as.numeric(format(as.Date(cy), "%Y"))
    } else {
      as.numeric(cy)
    }
  }]

  # Set column-name prefixes by mode.
  if (mode == "gap") {
    check_col <- "p22661"
    job_col_prefix <- NULL  # Processing note.
    start_col_prefix <- "p22663_a"
    end_col_prefix <- "p22664_a"
  } else {
    # work_category and work_exposure use the same prefix.
    check_col <- "p22599"
    job_col_prefix <- "p22601_a"
    start_col_prefix <- "p22602_a"
    end_col_prefix <- "p22603_a"
  }

  # Get and sort column names.
  get_sorted_cols <- function(prefix, data) {
    cols <- grep(paste0("^", prefix), names(data), value = TRUE)
    cols <- sort_ukb_cols(cols)
    return(cols)
  }

  job_cols <- get_sorted_cols(job_col_prefix, data)
  start_cols <- get_sorted_cols(start_col_prefix, data)
  end_cols <- get_sorted_cols(end_col_prefix, data)

  # Standardize column types.
  for (col in c(job_cols, start_cols, end_cols)) {
    if (col %in% names(data)) {
      data[, (col) := as.character(get(col))]
    }
  }

  # Reshape data.
  id_vars <- c("eid", check_col, "cutoff_year_numeric")

  if (mode == "gap") {
    measure_vars <- list(start_cols, end_cols)
    value_names <- c("start_year", "end_year")
  } else if (mode == "work_exposure") {
    # Add the exposure column.
    exp_cols_list <- lapply(exp_list, function(num) {
      get_sorted_cols(paste0("np", num, "_a"), data)
    })

    measure_vars <- c(list(job_cols, start_cols, end_cols), exp_cols_list)
    value_names <- c("work_content", "start_year", "end_year",
                     paste0("exp_", exp_list))
  } else {
    measure_vars <- list(job_cols, start_cols, end_cols)
    value_names <- c("work_content", "start_year", "end_year")
  }

  long_data <- melt(data,
                    id.vars = id_vars,
                    measure.vars = measure_vars,
                    variable.name = "job_index",
                    value.name = value_names)

  # Processing note.
  if (mode == "gap") {
    long_data <- long_data[!is.na(start_year) & start_year != "" & start_year != "NA"]
  } else {
    long_data <- long_data[!is.na(work_content) & work_content != "" & work_content != "NA"]
  }

  # Processing note.
  if (mode != "gap") {
    long_data[, unique_sig := {
      sapply(work_content, function(x) {
        match_result <- regmatches(x, regexpr("\\{([^}]+)\\}", x))
        if (length(match_result) > 0) {
          gsub("\\{|\\}", "", match_result)
        } else {
          NA_character_
        }
      })
    }]
  }

  # Process years.
  long_data[, start_year := suppressWarnings(as.numeric(start_year))]

  long_data[, end_year := {
    ey <- as.character(end_year)
    ey <- trimws(end_year)
    ifelse(is.na(ey) | ey == "" | ey == "NA",
           as.character(cutoff_year_numeric),
           ifelse(ey == "Ongoing when data entered",
                  as.character(cutoff_year_numeric),
                  ey))
  }]

  long_data[, end_year := suppressWarnings(as.numeric(end_year))]
  long_data[, cutoff_year_numeric := suppressWarnings(as.numeric(cutoff_year_numeric))]
  long_data[, end_year := pmin(end_year, cutoff_year_numeric, na.rm = TRUE)]

  # Filter and calculate duration.
  long_data <- long_data[!is.na(start_year) & start_year <= cutoff_year_numeric]
  long_data[, duration := pmax(0, end_year - start_year, na.rm = TRUE)]

  # Process differently by mode.
  if (mode == "work_category") {
    # Match work categories.
    if (!is.data.table(job_mapping)) {
      job_mapping <- as.data.table(job_mapping)
    }

    long_data[!is.na(unique_sig),
              work_category := job_mapping$work_category[match(unique_sig, job_mapping$unique_sig)]]
    long_data[is.na(unique_sig),
              work_category := job_mapping$work_category[match(work_content, job_mapping$work_content)]]

    # Aggregate.
    result <- long_data[!is.na(work_category),
                        .(total_years = sum(duration, na.rm = TRUE)),
                        by = .(eid, work_category)]

    # Convert to wide format.
    result_wide <- dcast(result, eid ~ work_category,
                         value.var = "total_years",
                         fill = 0)

  } else if (mode == "gap") {
    # Processing note.
    result <- long_data[, .(gap_time = sum(duration, na.rm = TRUE)), by = eid]
    result_wide <- result

  } else if (mode == "work_exposure") {
    # Processing note.
    exp_cols <- paste0("exp_", exp_list)

    # Convert exposure column to numeric.
    for (col in exp_cols) {
      long_data[, (col) := as.numeric(get(col))]
    }

    # Calculate weighted values.
    for (i in seq_along(exp_list)) {
      col <- exp_cols[i]
      out_col <- paste0("tidy_p", exp_list[i])
      long_data[, (out_col) := duration * get(col)]
    }

    # Aggregate.
    out_cols <- paste0("tidy_p", exp_list)
    result <- long_data[, lapply(.SD, function(x) sum(x, na.rm = TRUE)),
                        .SDcols = out_cols,
                        by = eid]
    result_wide <- result
  }

  # Handle cases where check_col is NA.
  all_eids <- data.table(eid = data$eid)
  # all_eids <- data.table(eid = data$eid,
  #                       has_check = !is.na(data[[check_col]]))
  result_wide <- merge(all_eids, result_wide, by = "eid", all.x = TRUE)

  return(result_wide)
}

###############################################################
## Parallel processing function.
##############################################################
setup_parallel_with_progress <- function(parallel_jobs = 1,
                                         max_global_size_gb = 50) {

  library(future)
  library(progressr)

  ## ========= Progress bar style. =========
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

    ## ========= BLAS / OMP =========
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }

    ## ========= future settings. =========
    options(future.globals.maxSize = max_global_size_gb * 1024^3)
    plan(multicore, workers = parallel_jobs)
  } else {
    handlers("progress")
    plan(sequential)
  }

  invisible(NULL)
}

###############################################################
## Diagnosis and death-outcome helper functions.
##############################################################
# For each diagnosis CSV file, calculate diagnosis-time information at each time point.
process_single_diag <- function(diag_csv_file, df_time, df_death, output_folder, case_thresh=100) {

  df <- fread(diag_csv_file,header = TRUE)

  df <- df |> filter(status != -1)

  df <- df |> mutate(eid = as.integer(eid))
  df_time <- df_time |> mutate(eid = as.integer(eid))

  df <- df |> left_join(df_time, by = "eid")
  df <- df |> left_join(df_death, by = "eid")

  time_columns <- setdiff(names(df_time), "eid")

  control_date <- as.Date("2023-03-31")

  for (time_col in time_columns) {
    new_col_name <- paste0("diag_time_", time_col)

    output_file_path <- file.path(glue("{output_folder}/each_diag_time_{time_col}"),
                                  basename(diag_csv_file))

    #if (file.exists(output_file)) {
    #  return(invisible(output_file))
    #}

    df_tmp <- df |>
      select(eid,status,date,all_of(time_col),death_date) |>
      drop_na(all_of(time_col))

    case_count <- sum(df_tmp$status == 1, na.rm = TRUE)

    if (case_count < case_thresh) {
        next
    }

    # censor_time: use the death date if it is available and earlier than control_date; otherwise use control_date.
    # Use if_else to preserve Date type and avoid ifelse converting Date values to numeric values.
    df_tmp$censor_time <- dplyr::if_else(
      !is.na(df_tmp$death_date) & df_tmp$death_date < control_date,
      as.Date(df_tmp$death_date),
      control_date
    )

    df_tmp[[new_col_name]] <- ifelse(
      df_tmp$status == 1,
      as.numeric(difftime(as.Date(df_tmp$date), as.Date(df_tmp[[time_col]]), units = "days")),
      as.numeric(difftime(df_tmp$censor_time, as.Date(df_tmp[[time_col]]), units = "days"))
    )

    df_tmp <- df_tmp |>
      select(eid,status,all_of(new_col_name))

    fwrite(df_tmp, output_file_path)
  }

  return(output_file_path)
}

# Convert all time columns to Date type.
# fread may convert date strings to days since epoch, usually around 14000-20000.
#   It may also parse them as Unix seconds, usually around 1.5e9, so convert based on value size.
safe_as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.character(x)) return(as.Date(x))
  if (is.numeric(x)) {
    max_val <- max(abs(x), na.rm = TRUE)
    if (max_val > 1e6) {
      # Unix timestamp in seconds: convert to POSIXct first, then extract Date.
      return(as.Date(as.POSIXct(x, origin = "1970-01-01")))
    } else {
      # Days since epoch.
      return(as.Date(x, origin = "1970-01-01"))
    }
  }
  return(as.Date(x))
}

# ========= ICD10 classification function =========
classify_icd10 <- function(code) {
  if (is.na(code) || code == "") return(NA_character_)
  first_char <- substr(code, 1, 1)
  second_char <- substr(code, 2, 2)

  if (first_char %in% c("A", "B")) return("Infectious")
  if (first_char == "C") return("Neoplasms")
  if (first_char == "D") {
    if (second_char %in% as.character(0:4)) return("Neoplasms")
    if (second_char %in% as.character(5:8)) return("Blood")
    return("Other")
  }
  if (first_char == "E") return("Endocrine")
  if (first_char == "F") return("Mental")
  if (first_char == "G") return("Nervous")
  if (first_char == "H") return("Ophthalmic_otic")
  if (first_char == "I") return("Circulatory")
  if (first_char == "J") return("Respiratory")
  if (first_char == "K") return("Digestive")
  if (first_char == "L") return("Dermatological")
  if (first_char == "M") return("Musculoskeletal")
  if (first_char == "N") return("Genitourinary")
  return("Other")
}


