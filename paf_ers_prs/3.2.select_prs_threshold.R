############################################################
## Select PRS thresholds in training data; score held-out participants
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
## Configure paths
############################################################

prs_dir <- here("processed_data/PRS_results")
split_path <- here("processed_data/ers/train_test_eids.rds")
covar_path <- here("processed_data/common_covariates.qs")
pca_path <- here("processed_data/genetic_pcs.tsv")
diag_path <- here("processed_data/diag_all_time/each_diag_time_p53_i0/all_diag.qs")
sex_vec_path <- here("processed_data/sex_specific_diseases.rds")
out_dir <- here("processed_data/prs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
parallel_jobs <- 10L
pc_vec <- paste0("p22009_a", 1:5)
covars <- c("Age_p53_i0", "Sex", "AccCeni0", pc_vec)

############################################################
## Read data and available disease-specific PRS files
############################################################

split <- readRDS(split_path)
covar_data <- as_tibble(qread(covar_path)) |>
  select(eid, Age_p53_i0, Sex, AccCeni0) |>
  inner_join(as_tibble(fread(pca_path)) |> select(eid, all_of(pc_vec)), by = "eid")
diag_list <- qread(diag_path)
sex_vec <- readRDS(sex_vec_path)
trait_dirs <- list.dirs(prs_dir, full.names = FALSE, recursive = FALSE)
diag_names <- intersect(sub("^prs_", "", trait_dirs[grepl("^prs_", trait_dirs)]), names(diag_list))

############################################################
## Select thresholds in parallel and build the test score matrix
############################################################

setup_parallel_with_progress(parallel_jobs)
results <- parallel_map(diag_names, function(diagnosis) {
  select_prs_threshold(diagnosis, prs_dir, split$train, split$test,
                        diag_list, covar_data, covars, sex_vec, min_snp = 20L)
})
future::plan(future::sequential)
results <- Filter(Negate(is.null), results)
selection <- bind_rows(lapply(results, `[[`, "selection"))
scores <- data.frame(eid = split$test)
for (result in results) {
  diagnosis <- result$selection$diag
  scores[[score_name(diagnosis)]] <- result$score$value
}
saveRDS(list(selection = selection, scores = scores), file.path(out_dir, "selected_prs.rds"))
