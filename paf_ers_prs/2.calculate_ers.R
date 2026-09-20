############################################################
## ERS: external categories, internal/external totals, lifestyle subtypes
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(survival)
  library(qs)
  library(here)
  library(future)
  library(future.apply)
  library(progressr)
})

source(here("paf_ers_prs/func_ers_prs.R"))

############################################################
## Configure input and output paths
############################################################

prepared_path <- here("processed_data/ers/prepared_exposures.qs")
covar_path <- here("processed_data/common_covariates.qs")
lifestyle_info_path <- here("rawdata/modifiable_exposure_categories.xlsx")
disease_info_path <- here("rawdata/disease_list.xlsx")
diag_path <- here("processed_data/diag_all_time/each_diag_time_p53_i0/all_diag_3yr.qs")
sex_vec_path <- here("processed_data/sex_specific_diseases.rds")
score_dir <- here("processed_data/ers/scores")
dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
parallel_jobs <- 10L
covars <- c("Age_p53_i0", "Sex", "AccCeni0")

############################################################
## Read training/test data and define the selected ERS groups
############################################################

prepared <- qread(prepared_path)
covar_data <- as_tibble(qread(covar_path)) |> select(eid, all_of(covars))
train_data <- as_tibble(prepared$train) |> inner_join(covar_data, by = "eid")
test_data <- as.data.frame(prepared$test)
lifestyle_info <- read_xlsx(lifestyle_info_path)
groups <- make_ers_groups(prepared$exp_info, lifestyle_info)

saveRDS(groups, file.path(score_dir, "group_manifest.rds"))
diag_list <- qread(diag_path)
disease_info <- read_xlsx(disease_info_path)
diag_names <- unique(c(intersect(na.omit(disease_info$disease_code), names(diag_list)),
                       grep("^ZZZDEATH_", names(diag_list), value = TRUE)))
sex_vec <- readRDS(sex_vec_path)

############################################################
## Correlation pruning, backward Cox selection and held-out ERS
############################################################

setup_parallel_with_progress(parallel_jobs)
for (i in seq_len(nrow(groups))) {
  group <- groups[i, ]
  exposures <- remove_ers_collinear_vars(
    train_data |> select(eid, all_of(group$exposures[[1]])), threshold = 0.8)
  message(group$category, " | retained before Cox: ", length(exposures))
  beta_results <- parallel_map(diag_names, function(diagnosis) {
    fit_ers_beta(train_data, exposures, diagnosis, diag_list, covars, sex_vec,
                 p_threshold = 0.05)
  })
  scores <- calculate_ers(test_data, beta_results)
  saveRDS(list(group_id = group$group_id, category = group$category,
               group_set = group$group_set, modifiable = group$modifiable,
               diag = diag_names, beta = beta_results, scores = scores),
          file.path(score_dir, paste0("ers_", group$group_id, ".rds")))
}
future::plan(future::sequential)
