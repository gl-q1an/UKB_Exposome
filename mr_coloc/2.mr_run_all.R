############################################################
## Run MR analyses for all clumped exposures and FinnGen outcomes
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(future)
  library(future.apply)
  library(progressr)
  library(ieugwasr)
  library(TwoSampleMR)
  library(readxl)
  library(glue)
  library(here)
  library(qs)
})

source(here("mr_coloc/step4_func.R"))

start_time <- Sys.time()

############################################################
## Load FinnGen outcome metadata
############################################################

diag_info <- read_xlsx(here("rawdata/Finndisease_MR_0335.xlsx"))
diag_info <- diag_info |>
  mutate(gwas_summary = coalesce(phenocode, phenocode2)) |>
  filter(!is.na(gwas_summary))

diag_info <- diag_info |>
  inner_join(diag_info_num, by = c("gwas_summary" = "phenocode"))

############################################################
## Load all clumped exposure instruments
############################################################

exp_all_dir <- here("processed_data/medelian_analysis/exp_sum_clump/all_1kg")
dir.create(exp_all_dir, showWarnings = FALSE, recursive = TRUE)
qs_file <- file.path(exp_all_dir, "all_clumped.qs")

base_dir <- here("processed_data/medelian_analysis/exp_sum_clump")
sub_dirs <- c("brainimage_1kg", "metab_1kg", "panukb_1kg")

if (!file.exists(qs_file)) {
  data_list <- list()
  counter <- 1

  for (sub in sub_dirs) {
    exp_dir <- file.path(base_dir, sub)

    tsv_files <- list.files(
      path = exp_dir,
      pattern = "\\.tsv$",
      full.names = FALSE
    )

    for (i in seq_along(tsv_files)) {
      file <- file.path(exp_dir, tsv_files[i])
      name <- tools::file_path_sans_ext(basename(file))
      pheno <- str_remove(name, ".EURsig_std_clumped")
      list_name <- paste(counter, pheno, sep = "_")

      data_list[[list_name]] <- fread(file)

      cat(counter, "Read exposure file:", file, "\n")
      counter <- counter + 1
    }
  }

  qsave(data_list, qs_file)
} else {
  message("Exposure list file already exists, reading: ", qs_file)
}

data_list_common <- qread(qs_file)
data_list_protein <- qread(here("processed_data/medelian_analysis/exp_sum_clump/protein/all_clumped.qs"))

common_names <- intersect(names(data_list_common), names(data_list_protein))

if (length(common_names) > 0) {
  message("Duplicated exposure names found:")
  print(common_names)
} else {
  message("No duplicated exposure names; merging common exposures and protein pQTL instruments.")
  data_list <- c(data_list_common, data_list_protein)
}

qsave(data_list, here("processed_data/medelian_analysis/exp_data_list.qs"))

snp_vec <- map(data_list, "SNP") |>
  unlist() |>
  unique()

############################################################
## Run MR for every FinnGen outcome against all exposure lists
############################################################

parallel_jobs <- 80
setup_parallel_with_progress(parallel_jobs)
all_mr_results <- list()
sched_strategy <- structure(1, ordering = "random")

for (i in seq_along(diag_info$NAME)) {
  outcome_diag <- diag_info$NAME[i]
  outcome_file <- diag_info$gwas_summary[i]

  cat("###### ", i, "/", nrow(diag_info), outcome_diag, " ######\n")
  cat("|| Begin at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " Ncores: ", parallel_jobs, "\n")

  tryCatch({
    out_dat <- fread(glue("{finn_sum_dir}/finngen_R12_{outcome_file}")) |>
      mutate(
        Nca = diag_info$num_cases[i],
        Nco = diag_info$num_controls[i],
        Phenotype = outcome_diag
      ) |>
      filter(rsids %in% snp_vec)

    if (nrow(out_dat) == 0) {
      cat("No SNPs found for ", outcome_diag, ", skipping...\n")
      next
    }

    outcome_dat <- format_data(
      as.data.frame(out_dat),
      type = "outcome",
      snp_col = "rsids",
      phenotype_col = "Phenotype",
      beta_col = "beta",
      se_col = "sebeta",
      eaf_col = "af_alt",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval",
      ncase_col = "Nca",
      ncontrol_col = "Nco",
      min_pval = 1e-200,
      log_pval = FALSE,
      chr_col = "#chrom",
      pos_col = "pos"
    )
  }, error = function(e) {
    cat("Error reading outcome data for ", outcome_diag, ": ", e$message, "\n")
    next
  })

  exposure_names <- names(data_list)

  res_list <- with_progress({
    p <- progressor(along = exposure_names)

    future_lapply(exposure_names, function(name) {
      res <- tryCatch({
        dat <- data_list[[name]]
        suppressMessages({
          if ("N" %in% names(dat)) {
            exp_dat <- format_data(
              as.data.frame(dat),
              type = "exposure",
              snp_col = "SNP",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "A1",
              other_allele_col = "A2",
              pval_col = "P",
              samplesize_col = "N",
              min_pval = 1e-200,
              log_pval = FALSE,
              chr_col = "CHR",
              pos_col = "BP"
            )
          } else if (all(c("Nca", "Nco") %in% names(dat))) {
            exp_dat <- format_data(
              as.data.frame(dat),
              type = "exposure",
              snp_col = "SNP",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "A1",
              other_allele_col = "A2",
              pval_col = "P",
              ncase_col = "Nca",
              ncontrol_col = "Nco",
              min_pval = 1e-200,
              log_pval = FALSE,
              chr_col = "CHR",
              pos_col = "BP"
            )
          }
        })

        out <- mr_single(exp_dat, outcome_dat, name, outcome_diag)

        if (is.null(out) || length(out) == 0) {
          out <- c(
            trait1 = name,
            trait2 = outcome_diag,
            nsnp = 0,
            method = "none",
            beta = -90,
            se = -90,
            beta_low95 = -90,
            beta_up95 = -90,
            pval = -90
          )
        }
        out
      }, error = function(e) {
        warning(sprintf("MR failed: %s (%s)", name, e$message))
        c(
          trait1 = name,
          trait2 = outcome_diag,
          nsnp = 0,
          method = "none",
          beta = -90,
          se = -90,
          beta_low95 = -90,
          beta_up95 = -90,
          pval = -90
        )
      })

      p(sprintf("%s -> %s", name, outcome_diag))
      res
    }, future.seed = TRUE, future.scheduling = sched_strategy)
  })

  mr_results_df <- do.call(rbind, lapply(res_list, function(x) {
    as.data.frame(t(x), stringsAsFactors = FALSE)
  }))

  numeric_cols <- c("nsnp", "beta", "se", "beta_low95", "beta_up95", "pval")
  mr_results_df[numeric_cols] <- lapply(mr_results_df[numeric_cols], as.numeric)

  all_mr_results[[outcome_diag]] <- mr_results_df

  cat("|| End at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("###############################\n\n")
}

############################################################
## Save MR results
############################################################

dir.create(here("results/mr"), showWarnings = FALSE, recursive = TRUE)
save(all_mr_results, file = here("tmp_files/final_mr_results_all_1kg.RData"))

final_results <- bind_rows(all_mr_results)
fwrite(final_results, here("results/mr/mr_results_all_1kg.csv"))

cat("\n========== MR Analysis Complete ==========\n")
cat("Total time:", format(Sys.time() - start_time), "\n")
cat("Total outcomes analyzed:", length(all_mr_results), "\n")
cat("Total MR tests:", nrow(final_results), "\n")