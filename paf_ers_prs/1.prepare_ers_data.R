############################################################
## Prepare training/test exposures for ERS
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
## Configure paths and analysis settings
############################################################

common_data_path <- here("processed_data/common_exposures.qs")
covar_path <- here("processed_data/common_covariates.qs")
exp_info_path <- here("rawdata/exposome_list.csv")
out_dir <- here("processed_data/ers")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seed <- 830L
train_fraction <- 1 / 3
missing_threshold <- 1 / 3

############################################################
## Read baseline exposures and split the discovery cohort
############################################################

covar_data <- as_tibble(qread(covar_path))
exp_info <- as_tibble(fread(exp_info_path))
exposure_data <- as_tibble(qread(common_data_path))

discover_covar <- covar_data |>
  filter(as.character(Ethnicity_Ca) == "1") |>
  select(eid, Sex) |>
  arrange(eid)
set.seed(seed)
train_eids <- sort(sample(discover_covar$eid, floor(nrow(discover_covar) * train_fraction)))
test_eids <- sort(setdiff(discover_covar$eid, train_eids))
saveRDS(list(train = train_eids, test = test_eids), file.path(out_dir, "train_test_eids.rds"))

exp_info <- exp_info |>
  filter(Time_variableid == "p53_i0", !is.na(exp_cate), exp_cate != "Proteome",
         exp_code %in% names(exposure_data))
raw_data <- discover_covar |>
  left_join(exposure_data |> select(eid, all_of(exp_info$exp_code)), by = "eid")
train_raw <- raw_data |> filter(eid %in% train_eids)
test_raw <- raw_data |> filter(eid %in% test_eids)

############################################################
## Fit preprocessing in training data, then apply it to test data
############################################################

missing_rate <- vapply(train_raw[exp_info$exp_code], function(x) mean(is.na(x)), numeric(1))
kept <- names(missing_rate)[missing_rate <= missing_threshold]
train_data <- data.frame(eid = train_raw$eid)
test_data <- data.frame(eid = test_raw$eid)
preprocess <- list()

for (exposure in kept) {
  coding_note <- exp_info$coding_note[match(exposure, exp_info$exp_code)]
  train_x <- convert_ers_exposure(train_raw[[exposure]], coding_note, exposure)
  test_x <- convert_ers_exposure(test_raw[[exposure]], coding_note, exposure)
  imputed <- impute_ers_with_train(train_x, test_x, train_raw$Sex, test_raw$Sex, coding_note)
  if (is.null(imputed)) next
  mu <- mean(imputed$train)
  sigma <- sd(imputed$train)
  if (!is.finite(sigma) || sigma == 0) next
  train_data[[exposure]] <- (imputed$train - mu) / sigma
  test_data[[exposure]] <- (imputed$test - mu) / sigma
  preprocess[[exposure]] <- list(mean = mu, sd = sigma, by_sex = imputed$by_sex,
                                 global_value = imputed$global_value)
}
exp_info <- exp_info |> filter(exp_code %in% names(train_data))
qsave(list(train = train_data, test = test_data, exp_info = exp_info,
           preprocess = preprocess, train_missing_rate = missing_rate),
      file.path(out_dir, "prepared_exposures.qs"))
