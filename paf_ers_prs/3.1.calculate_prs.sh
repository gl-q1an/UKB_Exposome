#!/bin/bash

############################################################
## Calculate disease-specific PRS at candidate P-value thresholds
############################################################

rscript_dir="tools/PRSice2"
summary_dir="rawdata/gwas_summaries_hg37"
ukb_dir="processed_data/plink_genotypes"
out_dir_main="processed_data/PRS_results"
file_list="rawdata/gwas_summary_files.txt"
threads=10

############################################################
## Run PRSice for each disease
############################################################

while IFS= read -r gwas_summary || [[ -n "$gwas_summary" ]]; do
  name=$(basename "$gwas_summary" _finn_r12_ukbsnp_hg37.tsv.gz)
  out_dir="${out_dir_main}/prs_${name}"
  mkdir -p "$out_dir"
  echo "Calculating PRS: ${name}"

  Rscript "${rscript_dir}/PRSice.R" \
    --prsice "${rscript_dir}/PRSice_linux" \
    --base "${summary_dir}/${gwas_summary}" \
    --snp uniq_id37 \
    --chr CHR \
    --bp BP_hg37 \
    --a1 A1 \
    --a2 A2 \
    --stat BETA \
    --beta \
    --pvalue P \
    --clump-kb 250kb \
    --clump-p 0.1 \
    --clump-r2 0.100000 \
    --fastscore \
    --bar-levels 5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1 \
    --no-regress \
    --target "${ukb_dir}/ukb_imp_chr#_v3" \
    --thread "$threads" \
    --binary-target T \
    --out "${out_dir}/prs_${name}_detail"
done < "$file_list"
