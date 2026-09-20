############################################################
## Mediation analysis using structural equation models
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(here)
  library(progressr)
  library(future)
  library(future.apply)
  library(glue)
  library(qs)
  library(lavaan)
})

source(here("mediation_analysis/func_mediation.R"))

start_time <- Sys.time()

############################################################
## Configure input and output paths
############################################################

external_data_path <- here("processed_data/ukb_all_external_data.qs")
internal_data_path <- here("processed_data/ukb_all_internal_data.qs")
covar_data_path <- here("processed_data/ukb_all_covar_common.qs")
diag_path <- here("processed_data/diag_all_time/each_diag_time_p53_i0/all_diag.qs")
mediation_pairs_path <- here("processed_data/mediation_analysis/mediation_pairs.tsv")

result_dir <- here("results/mediation_analysis")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(mediation_pairs_path), recursive = TRUE, showWarnings = FALSE)

parallel_jobs <- 10L
covariates <- c("Age_p53_i0", "Sex")

############################################################
## Load mediation pairs
############################################################

mediation_pairs <- data.table::fread(mediation_pairs_path)
required_pair_cols <- c("external_exposure", "internal_exposure", "diseases")
missing_pair_cols <- setdiff(required_pair_cols, names(mediation_pairs))

mediation_pairs <- mediation_pairs |>
  distinct(
    .data$external_exposure,
    .data$internal_exposure,
    .data$diseases,
    .keep_all = TRUE
  )

cat("Mediation pairs:", nrow(mediation_pairs), "\n")

############################################################
## Prepare Discover-only mediation data
############################################################

external_cols <- unique(mediation_pairs$external_exposure)
internal_cols <- unique(mediation_pairs$internal_exposure)

# Use the full Discover cohort, with the same definition as the ERS analysis.
# Ethnicity_Ca is used only for cohort selection, not as a model covariate.
covar_data <- qread(covar_data_path) |>
  filter(as.character(.data$Ethnicity_Ca) == "1") |>
  select(eid, all_of(covariates))

external_data <- qread(external_data_path) |>
  select(eid, any_of(external_cols))

internal_data <- qread(internal_data_path) |>
  select(eid, any_of(internal_cols))

data_all <- external_data |>
  inner_join(internal_data, by = "eid") |>
  inner_join(covar_data, by = "eid")

rm(external_data, internal_data, covar_data)
gc()

# Cache disease status tables so each SEM model only joins the required outcome.
diag_list <- qread(diag_path)
diag_status_cache <- load_diag_status_cache(diag_list, unique(mediation_pairs$diseases))
rm(diag_list)
gc()

setup_parallel_with_progress(parallel_jobs)

############################################################
## Run mediation models in the Discover cohort
############################################################

cat("\n==== Discover-only mediation models ====\n")
mediation_list <- run_mediation_params(
  params = mediation_pairs,
  data_all = data_all,
  diag_status_cache = diag_status_cache,
  covariates = covariates,
  progress_label = "Discover"
)

future::plan(future::sequential)

mediation_discover <- data.table::rbindlist(mediation_list, fill = TRUE) |>
  coerce_mediation_numeric_cols()

############################################################
## Save mediation results
############################################################

result_file <- save_mediation_result_rds(
  result_dir = result_dir,
  prefix = "mediation_sem_discover",
  result = mediation_discover
)

data.table::fwrite(mediation_discover, file.path(result_dir, "mediation_sem_discover.tsv"), sep = "\t")

cat("\n========== Mediation Analysis Complete ==========\n")
cat("Result file:", result_file, "\n")
cat("Rows:", nrow(mediation_discover), "\n")
cat("Total time:", format(Sys.time() - start_time), "\n")
cat("================================================\n")
