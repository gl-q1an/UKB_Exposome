####################################################
## Subgroup Cox association analysis (discovery only)
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

result_dir <- here("results/association_730_subgroup")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
subgroup_defs <- get_subgroup_defs()

####################################################
## Read shared data and define the two exposure sets
####################################################

exp_info <- fread(exp_info_path, na.strings = c("", "NA"))
disease_info <- read_xlsx(disease_info_path) |>
  filter(!is.na(disease_code)) |>
  distinct(disease_code, .keep_all = TRUE)
sex_vec <- readRDS(sex_specific_diag_path)

data_sets <- list(
  common = list(dis = common_discover_path, diag = "all_diag_3yr.qs"),
  protein = list(dis = protein_discover_path, diag = "all_diag_protein_3yr.qs")
)

####################################################
## Read data, form subgroups, and run Cox models in parallel
####################################################

setup_parallel_with_progress(parallel_jobs)
start_time <- Sys.time()

for (data_type in names(data_sets)) {
  inputs <- data_sets[[data_type]]
  discover_data <- as_tibble(qread(inputs$dis))
  params <- build_exposure_params(discover_data, exp_info,
                                   protein = data_type == "protein")
  subgroup_results <- list()

  for (time_col in unique(params$time_col)) {
    params_time <- params[params$time_col == time_col, ]
    diag_list <- qread(file.path(diag_dir, paste0("each_diag_time_", time_col), inputs$diag))
    diag_names <- get_cox_diag_names(diag_list, disease_info)

    for (subgroup in subgroup_defs) {
      subgroup_params <- adjust_params_for_subgroup(params_time, subgroup)
      pairs <- expand_cox_pairs(subgroup_params, diag_names)
      message(data_type, " | ", time_col, " | ", subgroup$label)
      result <- run_cox_pairs(discover_data, pairs, diag_list, sex_vec, subgroup)
      result$subgroup_type <- subgroup$type
      result$subgroup_level <- subgroup$level
      result$subgroup_label <- subgroup$label
      subgroup_results[[paste(time_col, subgroup$label, sep = "_")]] <- result
    }
    rm(diag_list)
    gc(verbose = FALSE)
  }

  saveRDS(list(dis = bind_rows(subgroup_results)),
          file.path(result_dir, paste0(data_type, "_subgroup_cox.rds")))
  rm(discover_data, subgroup_results)
  gc(verbose = FALSE)
}
future::plan(future::sequential)

cat("Subgroup Cox analysis finished in", format(Sys.time() - start_time), "\n")
# The >= 100 incident-case threshold is applied separately within each subgroup
# after complete-case filtering. Subgroup estimates are not interaction tests.
