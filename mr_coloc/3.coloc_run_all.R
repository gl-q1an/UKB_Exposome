############################################################
## Prepare exposure regions and run coloc analyses
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(stringr)
  library(glue)
  library(here)
  library(qs)
  library(coloc)
  library(future)
  library(future.apply)
  library(progressr)
})

source(here("mr_coloc/step4_func.R"))

start_time <- Sys.time()

############################################################
## Configure shared input and output paths
############################################################

mr_file <- here("results/mr/MR_res_significant.xlsx")
data_list_file <- here("processed_data/medelian_analysis/exp_data_list.qs")
coloc_dir <- here("processed_data/coloc_exp_region")
result_dir <- here("results/coloc")
dir.create(coloc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

window_size <- 500000L
parallel_jobs <- 20L
coloc_parallel_jobs <- 80L

############################################################
## Load MR results and exposure instrument list
############################################################

cat("\n[1] Reading MR results and exposure instruments...\n")

mr_result <- read_coloc_mr_result(mr_file)
data_list <- qs::qread(data_list_file)

############################################################
## Format Pan-UKB exposure regions
############################################################

cat("\n[2] Formatting Pan-UKB exposure regions...\n")

panukb_loop <- get_coloc_exposure_table(mr_result, "panukb")
print_coloc_exposure_summary(mr_result, panukb_loop, "panukb")

setup_parallel_with_progress(parallel_jobs)
for (i in seq_len(nrow(panukb_loop))) {
  filename0 <- panukb_loop$trait1[i]
  trait1_rawcol <- panukb_loop$trait1_rawcol[i]

  cat("###### Pan-UKB ", i, "/", nrow(panukb_loop), " ", filename0, " ######\n")
  outpath <- format_panukb_one(
    filename0 = filename0,
    trait1_rawcol = trait1_rawcol,
    data_list = data_list,
    pan_ukb = pan_ukb,
    input_dir = panukb_input_dir,
    output_dir = coloc_dir,
    window_size = window_size
  )
  if (!is.null(outpath)) cat("Saved:", outpath, "\n")
}

############################################################
## Format metabolomics exposure regions
############################################################

cat("\n[3] Formatting metabolomics exposure regions...\n")

metab_loop <- get_coloc_exposure_table(mr_result, "metab")
print_coloc_exposure_summary(mr_result, metab_loop, "metab")

setup_parallel_with_progress(parallel_jobs)
for (i in seq_len(nrow(metab_loop))) {
  filename0 <- metab_loop$trait1[i]
  trait1_rawcol <- metab_loop$trait1_rawcol[i]

  cat("###### Metabolomics ", i, "/", nrow(metab_loop), " ", filename0, " ######\n")
  outpath <- format_metab_one(
    filename0 = filename0,
    trait1_rawcol = trait1_rawcol,
    data_list = data_list,
    input_dir = metab_input_dir,
    output_dir = coloc_dir,
    window_size = window_size
  )
  if (!is.null(outpath)) cat("Saved:", outpath, "\n")
}

############################################################
## Format brain imaging exposure regions
############################################################

cat("\n[4] Formatting brain imaging exposure regions...\n")

idp_loop <- get_coloc_exposure_table(mr_result, "IDP")
print_coloc_exposure_summary(mr_result, idp_loop, "IDP")

setup_parallel_with_progress(parallel_jobs)
for (i in seq_len(nrow(idp_loop))) {
  filename0 <- idp_loop$trait1[i]
  trait1_rawcol <- idp_loop$trait1_rawcol[i]

  cat("###### IDP ", i, "/", nrow(idp_loop), " ", filename0, " ######\n")
  outpath <- format_idp_one(
    filename0 = filename0,
    trait1_rawcol = trait1_rawcol,
    data_list = data_list,
    idp_info = idp_info,
    idp_snplist = idp_snplist,
    input_dir = idp_input_dir,
    output_dir = coloc_dir,
    window_size = window_size
  )
  if (!is.null(outpath)) cat("Saved:", outpath, "\n")
}

############################################################
## Format protein pQTL exposure regions
############################################################

cat("\n[5] Formatting protein pQTL exposure regions...\n")

pqtl_loop <- get_coloc_exposure_table(mr_result, "pqtl")
print_coloc_exposure_summary(mr_result, pqtl_loop, "pqtl")

for (i in seq_len(nrow(pqtl_loop))) {
  protein <- pqtl_loop$trait1[i]

  cat("###### pQTL ", i, "/", nrow(pqtl_loop), " ", protein, " ######\n")
  outpath <- format_protein_one(
    protein = protein,
    protein_map = protein_map,
    pgwas_dir = pgwas_dir,
    output_dir = coloc_dir,
    window_size = window_size,
    work_dir = here()
  )
  if (!is.null(outpath)) cat("Saved:", outpath, "\n")
}

############################################################
## Run coloc by streaming one FinnGen outcome at a time
############################################################

cat("\n[6] Running coloc analyses...\n")

run_coloc_by_outcome_streamed(
  mr_file = mr_file,
  target_classes = c("metab", "IDP", "panukb", "pqtl"),
  coloc_dir = coloc_dir,
  result_dir = result_dir,
  result_prefix = "coloc_all",
  pval_threshold = Inf,
  parallel_jobs = coloc_parallel_jobs,
  window_size = window_size,
  finn_sum_dir = "Finngen_GWAS_summary/R12",
  data_list_file = data_list_file,
  finn_snplist_file = here("processed_data/finngen_r12/finngen_r12_snplist_hg37_af.tsv"),
  diag_info_file = here("rawdata/Finndisease_MR_0335.xlsx"),
  finngen_manifest_file = "Finngen_GWAS_summary/finngen_R12_manifest.tsv",
  skip_existing = TRUE,
  preload_exposure_regions = TRUE,
  preload_scope = "outcome"
)

cat("\n========== Coloc Workflow Complete ==========\n")
cat("Exposure region dir:", coloc_dir, "\n")
cat("Result dir:", result_dir, "\n")
cat("Total time:", format(Sys.time() - start_time), "\n")