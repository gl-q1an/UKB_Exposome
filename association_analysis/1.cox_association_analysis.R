####################################################
## Cox association analysis: discovery and replication
####################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(survival)
  library(here)
  library(progressr)
  library(future)
  library(future.apply)
  library(qs)
})

source(here("association_analysis/func_association.R"))

####################################################
## Configure paths
####################################################

common_discover_path <- here("processed_data/ukb_all_common_discover_scaled.qs")
protein_discover_path <- here("processed_data/ukb_all_protein_discover_scaled.qs")
exp_info_path <- here("rawdata/exposome_list_730.csv")
disease_info_path <- here("rawdata/disease_list.xlsx")
sex_specific_diag_path <- here("processed_data/sex_specific_diseases.rds")
diag_dir <- here("processed_data/diag_all_time")
parallel_jobs <- 10L

common_replication_path <- here("processed_data/ukb_all_common_replication_scaled.qs")
protein_replication_path <- here("processed_data/ukb_all_protein_replication_scaled.qs")
result_dir <- here("results/association")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

####################################################
## Read shared data and define the two exposure sets
####################################################

exp_info <- fread(exp_info_path, na.strings = c("", "NA"))
disease_info <- read_xlsx(disease_info_path) |>
  filter(!is.na(disease_code)) |>
  distinct(disease_code, .keep_all = TRUE)
sex_vec <- readRDS(sex_specific_diag_path)

data_sets <- list(
  common = list(dis = common_discover_path, rep = common_replication_path,
                diag_dis = "all_diag_3yr.qs", diag_rep = "all_diag.qs"),
  protein = list(dis = protein_discover_path, rep = protein_replication_path,
                 diag_dis = "all_diag_protein_3yr.qs", diag_rep = "all_diag_protein.qs")
)

####################################################
## Read data, build Cox pairs, and run in parallel
####################################################

setup_parallel_with_progress(parallel_jobs)
start_time <- Sys.time()

for (data_type in names(data_sets)) {
  inputs <- data_sets[[data_type]]
  discover_data <- as_tibble(qread(inputs$dis))
  replication_data <- as_tibble(qread(inputs$rep))

  params <- build_exposure_params(discover_data, exp_info,
                                   protein = data_type == "protein")
  results_dis <- results_rep <- list()

  # Time-point grouping ensures every exposure uses the correct origin date.
  for (time_col in unique(params$time_col)) {
    params_time <- params[params$time_col == time_col, ]
    time_dir <- file.path(diag_dir, paste0("each_diag_time_", time_col))
    diag_dis <- qread(file.path(time_dir, inputs$diag_dis))
    diag_rep <- qread(file.path(time_dir, inputs$diag_rep))
    diag_names <- get_cox_diag_names(diag_dis, disease_info)
    pairs <- expand_cox_pairs(params_time, diag_names)
    message(data_type, " | ", time_col, " | Cox pairs: ", nrow(pairs))

    results_dis[[time_col]] <- run_cox_pairs(discover_data, pairs, diag_dis, sex_vec)
    results_rep[[time_col]] <- run_cox_pairs(replication_data, pairs, diag_rep, sex_vec)
    rm(diag_dis, diag_rep)
    gc(verbose = FALSE)
  }

  saveRDS(list(dis = bind_rows(results_dis), rep = bind_rows(results_rep)),
          file.path(result_dir, paste0(data_type, "_cox.rds")))
  rm(discover_data, replication_data, results_dis, results_rep)
  gc(verbose = FALSE)
}
future::plan(future::sequential)

cat("Cox analysis finished in", format(Sys.time() - start_time), "\n")
# Retain all rows for auditing. Use success == "1" and valid cox_pval
# for downstream inference; review warnings and choose multiplicity correction
# for the intended family of hypotheses before declaring associations.
