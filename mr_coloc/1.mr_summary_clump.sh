#!/bin/bash

set -euo pipefail

############################################################
## LD clumping for MR exposure summary statistics
############################################################

N_JOBS=20

BASE_DIR="${WORK_DIR}/processed_data/medelian_analysis"
PANUKB_IN="${BASE_DIR}/exp_sum/panukb"
PANUKB_OUT="${BASE_DIR}/exp_sum_clump/panukb_1kg"

mkdir -p "${PANUKB_OUT}"

combine_clump() {
  local clumped_file="$1"
  local input_file="$2"
  local output_file="$3"

  unset R_HOME
  Rscript -e "source('${FUNC_FILE}'); combine_clump_with_summary('${clumped_file}', '${input_file}', '${output_file}')"
}

export -f combine_clump
export R_ENV R_ENV_NAME FUNC_FILE

############################################################
## Clumping
############################################################

find "${PANUKB_IN}" -maxdepth 1 -name "*.EURsig_std_1kg.tsv" | sort | parallel --bar -j "${N_JOBS}" '
  name=$(basename {} .EURsig_std_1kg.tsv)
  echo "===== Begin Pan-UKB ${name} ====="
  plink2 \
    --bfile "'"${REF_BFILE}"'" \
    --clump {} \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-kb 10000 \
    --clump-snp-field SNP \
    --clump-field P \
    --threads 1 \
    --out "'"${PANUKB_OUT}"'/${name}_clump_1kg"
  combine_clump "'"${PANUKB_OUT}"'/${name}_clump_1kg.clumps" "{}" "'"${PANUKB_OUT}"'/${name}.EURsig_std_clumped_1kg.tsv"
  echo "===== End Pan-UKB ${name} ====="
'
