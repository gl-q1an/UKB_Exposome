############################################################
## Weighted PAF: six modifiable categories and intrinsic/extrinsic ERS + PRS
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
  library(psych)
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
out_dir <- here("results/weighted_paf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
parallel_jobs <- 10L
base_covars <- c("Age_p53_i0", "Sex", "AccCeni0")
pc_vec <- paste0("p22009_a", 1:5)
prs_covars <- c(base_covars, pc_vec)

groups <- readRDS(file.path(score_dir, "group_manifest.rds")) |>
  filter(modifiable | group_set == "inout")
ers_inputs <- setNames(lapply(groups$group_id, function(group_id) {
  readRDS(file.path(score_dir, paste0("ers_", group_id, ".rds")))
}), groups$group_id)
prs <- readRDS(prs_path)
score_tables <- lapply(ers_inputs, `[[`, "scores")
prs_scores <- as.data.frame(prs$scores)

prs_scores <- prs_scores[rowSums(!is.na(prs_scores[setdiff(names(prs_scores), "eid")])) > 0, ]
score_tables$prs <- prs_scores
covar_data <- as_tibble(qread(covar_path)) |> select(eid, all_of(base_covars))
prs_covar <- covar_data |>
  inner_join(as_tibble(fread(pca_path)) |> select(eid, all_of(pc_vec)), by = "eid")
diag_list <- qread(diag_path)
sex_vec <- readRDS(sex_vec_path)

setup_parallel_with_progress(parallel_jobs)
component_results <- list()
for (key in names(score_tables)) {
  scores <- score_tables[[key]]
  diagnoses <- if (key == "prs") prs$selection$diag else ers_inputs[[key]]$diag
  diagnoses <- intersect(diagnoses, names(diag_list))
  covars <- if (key == "prs") prs_covars else base_covars
  covariates <- if (key == "prs") prs_covar else covar_data
  category <- if (key == "prs") "PRS" else groups$category[match(key, groups$group_id)]
  components <- parallel_map(diagnoses, function(diagnosis) {
    model_data <- scores |>
      select(eid, score = all_of(score_name(diagnosis))) |>
      inner_join(covariates, by = "eid")
    component_af(model_data, diagnosis, diag_list, covars, sex_vec)
  })
  component_results[[key]] <- bind_rows(components) |>
    mutate(category_key = key, category = category,
           source_population = if (key == "prs") "test_with_prs_and_pcs" else "test")
}
component_table <- bind_rows(component_results)

component_sets <- list(
  modifiable = groups$group_id[groups$modifiable],
  inout_prs = c("in", "out", "prs")
)
weighted_results <- list()
for (set_name in names(component_sets)) {
  keys <- component_sets[[set_name]]
  reference_eids <- sort(Reduce(intersect, lapply(score_tables[keys], `[[`, "eid")))
  diagnosis_sets <- lapply(component_results[keys], function(x) x$diag)
  diagnoses <- sort(Reduce(intersect, diagnosis_sets))
  results <- parallel_map(diagnoses, function(diagnosis) {
    aligned_scores <- data.frame(eid = reference_eids)
    for (key in keys) {
      table <- score_tables[[key]]
      aligned_scores[[key]] <- table[[score_name(diagnosis)]][match(reference_eids, table$eid)]
    }
    components <- component_table |>
      filter(diag == diagnosis, category_key %in% keys)
    weighted_paf(aligned_scores, components)
  })
  weighted_results[[set_name]] <- bind_rows(results) |>
    mutate(group_set = set_name)
}
future::plan(future::sequential)

############################################################
## Save weighted PAF outputs and their calculation inputs
############################################################

saveRDS(list(modifiable = weighted_results$modifiable,
             inout_prs = weighted_results$inout_prs,
             component_inputs = component_table),
        file.path(out_dir, "weighted_paf_results.rds"))
fwrite(weighted_results$modifiable,
       file.path(out_dir, "modifiable_weighted_paf.tsv"), sep = "\t")
fwrite(weighted_results$inout_prs,
       file.path(out_dir, "inout_prs_weighted_paf.tsv"), sep = "\t")
