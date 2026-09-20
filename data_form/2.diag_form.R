suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(lubridate)
  library(progressr)
  library(future)
  library(future.apply)
  library(readxl)
  library(here)
  library(glue)
  library(qs)
})

source(here("data_form/step1_func.R"))

################################################################################
# Diagnosis follow-up time files
################################################################################

# Create the folder used to store diagnosis files.
output_dir <- here("processed_data/diag_all_time")

df_time <- fread('processed_data/ukb_time.tsv')

# Read death dates as censoring dates for controls; use death date when available, otherwise use control_date 2023-03-31.
df_death <- readRDS(here('processed_data/ukb_death_info.RDS'))
df_death <- df_death |>
  mutate(eid = as.integer(eid),
         death_date = as.Date(p40000_i0)) |>
  select(eid, death_date)

time_columns <- setdiff(names(df_time), "eid")
for (time_col in time_columns) {
  dir.create(here(glue("{output_dir}/each_diag_time_{time_col}")),recursive = TRUE, showWarnings = FALSE)
}

# ========= Set variables =========
diag_dir <- here("ukb_diag_by_finn251226/each_diag")

diag_info <- read_xlsx(here("rawdata/ukb_disease_525.xlsx"))

diag_vec <- paste0(diag_dir,'/',diag_info$disease_code,'.csv')

samp_thresh <- 100

df_time <- df_time %>%
  mutate(across(-1, as.Date))

# ========= Set up parallel processing =========
parallel_jobs <- 20 # Adjust according to the number of CPU cores.
setup_parallel_with_progress(parallel_jobs = parallel_jobs, max_global_size_gb = 50)

# ========= Run parallel processing =========
with_progress({
  p <- progressor(along =  diag_vec)

  results <- future_lapply(diag_vec,function(diag_file) {

      res <- process_single_diag(
        diag_csv_file = diag_file,
        df_time       = df_time,
        df_death      = df_death,
        output_folder = output_dir,
        case_thresh   = samp_thresh
      )

      p(sprintf("%s", basename(diag_file)))
      res
    },future.seed = TRUE
  )
})

# One issue still needs attention after this step.
# Some controls can also have negative diag_time values.
# This means the assessment time for that item is after the censoring time.
# These participants can still be treated as controls in logistic models.
# However, they should not be included in Cox model calculations.

################################################################################
# Death outcome files
################################################################################
control_date <- as.Date("2023-03-31")

# ========= Read data =========
df_death <- readRDS(here("processed_data/ukb_death_info.RDS"))
df_death <- df_death[, .(
  eid = as.integer(eid),
  death_date = as.Date(p40000_i0),
  icd10 = as.character(p40001_i0)
)]

df_time <- fread(here("processed_data/ukb_time.tsv"))
time_columns <- setdiff(names(df_time), "eid")

df_time[, (time_columns) := lapply(.SD, safe_as_date), .SDcols = time_columns]

cat("[-] Time columns:", length(time_columns), "\n")
cat("[-] Death records:", nrow(df_death),
    "(died:", sum(!is.na(df_death$death_date)),
    "alive:", sum(is.na(df_death$death_date)), ")\n")

# ========= Create output directories =========
output_dir <- here("processed_data/diag_all_time")
for (time_col in time_columns) {
  dir.create(glue("{output_dir}/each_diag_time_{time_col}"),
             recursive = TRUE, showWarnings = FALSE)
}

df_death[, cause_category := sapply(icd10, classify_icd10)]

cat("[-] ICD10 category distribution (among deaths):\n")
print(table(df_death$cause_category, useNA = "ifany"))

cause_categories <- c("Infectious", "Neoplasms", "Blood", "Endocrine",
                      "Mental", "Nervous", "Ophthalmic_otic", "Circulatory",
                      "Respiratory", "Digestive", "Dermatological",
                      "Musculoskeletal", "Genitourinary", "Other")

# ========= Process each time point =========
for (time_col in time_columns) {
  cat("\n[*] Processing:", time_col, "\n")

  # Extract data for this time point.
  df_tmp <- df_time[, .(eid = as.integer(eid), baseline_time = get(time_col))]
  df_tmp <- df_tmp[!is.na(baseline_time)]

  # Left-join death information.
  df_tmp <- merge(df_tmp, df_death, by = "eid", all.x = TRUE)

  # Calculate all-cause death status: death date is non-missing and before the censoring date.
  df_tmp[, allcause_status := ifelse(!is.na(death_date) & death_date <= control_date, 1L, 0L)]

  # Calculate the censoring endpoint: cases use death date and controls use censor date.
  df_tmp[, censor_end := as.Date(ifelse(allcause_status == 1,
                                        as.character(death_date),
                                        as.character(control_date)))]

  # Calculate diag_time in days.
  df_tmp[, diag_time := as.numeric(difftime(censor_end, baseline_time, units = "days"))]

  # ===== All-cause death =====
  out_allcause <- df_tmp[, .(eid, status = allcause_status, diag_time)]
  new_col_name <- paste0("diag_time_", time_col)
  setnames(out_allcause, "diag_time", new_col_name)

  fwrite(out_allcause, file.path(output_dir, glue("each_diag_time_{time_col}"), "ZZZDEATH_AllCause.csv"))
  cat("  [AllCause] cases:", sum(out_allcause$status == 1),
      "controls:", sum(out_allcause$status == 0), "\n")

  # ===== Cause-specific death by broad categories =====
  for (cause in cause_categories) {
    # tmp_status logic:
    #   Start from allcause_status.
    #   If death occurred from any cause, first mark it as -1 for exclusion.
    #   Then recover deaths from the current cause and mark them as 1.
    df_tmp[, tmp_status := allcause_status]
    df_tmp[allcause_status == 1, tmp_status := -1L]
    df_tmp[allcause_status == 1 & cause_category == cause, tmp_status := 1L]

    # Exclude participants who died from other causes.
    df_cause <- df_tmp[tmp_status != -1, .(eid, status = tmp_status, diag_time)]
    setnames(df_cause, "diag_time", new_col_name)

    case_n <- sum(df_cause$status == 1)
    ctrl_n <- sum(df_cause$status == 0)
    excl_n <- sum(df_tmp$tmp_status == -1)

    cat(sprintf("  [%s] cases: %d, controls: %d, excluded(competing death): %d\n",
                cause, case_n, ctrl_n, excl_n))

    fwrite(df_cause, file.path(output_dir, glue("each_diag_time_{time_col}"),
                               glue("ZZZDEATH_{cause}.csv")))
  }
}

cat("\n[OK] Death outcome files generated successfully!\n")

# ========= Generate death-count summary table =========
cat("\n[*] Generating death summary table...\n")

death_dir <- here("processed_data/death_all_time")
time_dirs <- list.files(death_dir, pattern = "^each_diag_time_", full.names = TRUE)
time_dirs <- sort(time_dirs)
time_cols <- gsub("^each_diag_time_", "", basename(time_dirs))

cat_info <- data.table(
  category = c("AllCause", "Infectious", "Neoplasms", "Blood", "Endocrine",
               "Mental", "Nervous", "Ophthalmic_otic", "Circulatory",
               "Respiratory", "Digestive", "Dermatological",
               "Musculoskeletal", "Genitourinary", "Other"),
  icd10_code = c("All", "A*, B*", "C*, D0-D4*", "D5-D8*", "E*", "F*", "G*", "H*",
                 "I*", "J*", "K*", "L*", "M*", "N*",
                 "U*, V*, W*, X*, Y*, D9*, O*, P*, Q*, R*")
)

result <- copy(cat_info)
setnames(result, "category", "Death Category")
setnames(result, "icd10_code", "ICD10 Code")

for (i in seq_along(time_dirs)) {
  tc <- time_cols[i]
  dir_path <- time_dirs[i]

  counts <- sapply(cat_info$category, function(cat) {
    fpath <- file.path(dir_path, paste0(cat, ".csv"))
    if (file.exists(fpath)) {
      df <- fread(fpath)
      sum(df$status == 1)
    } else {
      NA_integer_
    }
  })

  result[[tc]] <- counts
}

cat("\n=== Death Outcome Summary (Censor:", as.character(control_date), ") ===\n")
cat(sprintf("Total N (first time col): %d\n\n", nrow(fread(file.path(time_dirs[1], "AllCause.csv")))))
print(result, nrows = Inf)

fwrite(result, file.path(death_dir, "death_summary.tsv"), sep = "\t")
cat("\n[OK] Summary saved to processed_data/death_summary.tsv\n")

################################################################################
# Combine diagnosis files
################################################################################
# Set the folder path.

all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

for (folder_path in all_dirs) {

    # Get all CSV files in the folder.
    csv_files <- list.files(path = folder_path, 
                            pattern = "\\.csv$", 
                            full.names = TRUE)

    # Create an empty list to store data.
    data_list <- list()

    # Read all CSV files.
    for (i in seq_along(csv_files)) {
        file <- csv_files[i]
        # Get the base name and remove the .csv extension.
        name <- tools::file_path_sans_ext(basename(file))
        
        # Read the CSV file and store it in the list.
        data_list[[name]] <- fread(file)
        
        # Print progress.
        cat(i,'/', length(csv_files),",read:", name, "\n")
    }

    # Save in qs format.
    qsave(data_list, paste0(folder_path,"/all_diag.qs"))
}