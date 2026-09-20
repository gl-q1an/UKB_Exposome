############################################################
## Cox main effects and multiplicative ERS x PRS interaction
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
## Configure paths and covariates
############################################################

score_dir <- here("processed_data/ers/scores")
prs_path <- here("processed_data/prs/selected_prs.rds")
covar_path <- here("processed_data/common_covariates.qs")
pca_path <- here("processed_data/genetic_pcs.tsv")
diag_path <- here("processed_data/diag_all_time/each_diag_time_p53_i0/all_diag_3yr.qs")
sex_vec_path <- here("processed_data/sex_specific_diseases.rds")
out_dir <- here("results/ers_prs_cox")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
parallel_jobs <- 10L
base_covars <- c("Age_p53_i0", "Sex", "AccCeni0")
pc_vec <- paste0("p22009_a", 1:5)
interaction_covars <- c(base_covars, pc_vec)

############################################################
## Read held-out scores, covariates and incident outcomes
############################################################

groups <- readRDS(file.path(score_dir, "group_manifest.rds"))
prs <- readRDS(prs_path)
covar_data <- as_tibble(qread(covar_path)) |> select(eid, all_of(base_covars))
interaction_covar <- covar_data |>
  inner_join(as_tibble(fread(pca_path)) |> select(eid, all_of(pc_vec)), by = "eid")
prs_scores <- as.data.frame(prs$scores)
prs_scores <- prs_scores |> filter(eid %in% interaction_covar$eid)
prs_scores <- prs_scores[rowSums(!is.na(prs_scores[setdiff(names(prs_scores), "eid")])) > 0, ]
prs_scores <- standardize_scores(prs_scores)
diag_list <- qread(diag_path)
sex_vec <- readRDS(sex_vec_path)
main_results <- interaction_results <- list()

############################################################
## Fit Cox models in parallel
############################################################

setup_parallel_with_progress(parallel_jobs)
for (i in seq_len(nrow(groups))) {
  group <- groups[i, ]
  ers <- readRDS(file.path(score_dir, paste0("ers_", group$group_id, ".rds")))
  ers_scores <- standardize_scores(ers$scores)
  results <- parallel_map(ers$diag, function(diagnosis) {
    column <- score_name(diagnosis)
    ers_data <- ers_scores |> select(eid, ers = all_of(column))
    main <- fit_score_cox(inner_join(ers_data, covar_data, by = "eid"),
                          diagnosis, diag_list, base_covars, sex_vec)
    interaction <- NULL
    if (column %in% names(prs_scores)) {
      model_data <- ers_data |>
        inner_join(prs_scores |> select(eid, prs = all_of(column)), by = "eid") |>
        inner_join(interaction_covar, by = "eid")
      interaction <- fit_score_cox(model_data, diagnosis, diag_list,
                                   interaction_covars, sex_vec, interaction = TRUE)
    }
    list(main = main, interaction = interaction)
  })
  main_results[[group$group_id]] <- bind_rows(lapply(results, `[[`, "main")) |>
    mutate(group_id = group$group_id, category = group$category, modifiable = group$modifiable)
  interaction_results[[group$group_id]] <- bind_rows(lapply(results, `[[`, "interaction")) |>
    mutate(group_id = group$group_id, category = group$category, modifiable = group$modifiable)
}
future::plan(future::sequential)
main_result <- bind_rows(main_results)
interaction_result <- bind_rows(interaction_results)
saveRDS(list(main = main_result, interaction = interaction_result),
        file.path(out_dir, "cox_results.rds"))

interaction_summary <- interaction_result |>
  filter(modifiable, term == "ers:prs", is.finite(pval), pval >= 0, pval <= 1) |>
  mutate(significant = pval < 0.05) |>
  count(category, significant, name = "n_diseases")
fwrite(interaction_summary, file.path(out_dir, "modifiable_interaction_counts.tsv"), sep = "\t")
