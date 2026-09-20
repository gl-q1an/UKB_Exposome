suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(purrr)
  library(stringr)
  library(here)
  library(qs)
  library(lubridate)
})

source(here("data_form/step1_func.R"))

################################################################################
# Common continuous variables
################################################################################

######################################################
## Tidy variables according to the data-processing codebook.
## The codebook mainly covers common-data variables.
## Process regular variables first.
## Handle variables that aggregate multiple fields later.
#####################################################

ukb_val_info <- read_xlsx(
  here("rawdata/ukb_exposome.xlsx"),
  col_types = c(
    "field_id" = "text"
  )
) |>
  filter(exp_cate != "Proteome") |> # Exclude proteomics.
  filter((exp_subcate!='Brain MRI' | is.na(exp_subcate))) |>  # Exclude brain MRI.
  filter(Level2 != "Derived accelerometry" | is.na(Level2)) |># Exclude accelerometry, which has already been processed.
  filter(exp_cate != "Metabolome") |> # Exclude metabolomics.
  filter(!grepl("\\|", field_id)) |>  # Special cases handled manually.
  mutate(field_id = as.numeric(field_id))

ukb_val_code <- read_xlsx(here("rawdata/ukb_codebook_0119.xlsx"), col_types = c("field_id" = "text")) |>
  filter(!grepl("\\|", field_id)) |>
  mutate(field_id = as.numeric(field_id))

ukb_val_code_info <- ukb_val_code |>
  drop_na(title) 
table(duplicated(ukb_val_code_info$unique_id))

ukb_sum <- read.csv(here("rawdata/ukb_data_summary.csv"))

df_ukb_common <- fread(here("processed_data/ukb_common.tsv"))

# Merge extracted columns with the metadata table to recover original UKB column names.
# Use unique_id as the new tidy variable name.

field_id_table <- data.frame(
  ukb_colname = names(df_ukb_common)[-1],
  field_id = str_extract(names(df_ukb_common)[-1], "(?<=p)\\d+"),
  stringsAsFactors = FALSE
) |>
  group_by(field_id) |>
  summarise(
    ukb_colname = paste(str_sort(ukb_colname, numeric = TRUE), collapse = "|"),
    .groups = "drop"
) |>
  mutate(field_id = as.numeric(field_id))

ukb_val_info <- ukb_val_info |>
  left_join(field_id_table, by='field_id')

table(is.na(ukb_val_info$ukb_colname))

# Create detailed type groups to make manual checks easier.

dir.create(here("tmp_files/ukb_data_extract"),
           recursive = TRUE, showWarnings = FALSE)

#######################################################################
# 1. Raw data are continuous or integer, coding is 0, and fields are not arrays.
#######################################################################

con_0code_0arr_vec <- ukb_sum |>
  filter((value_type %in% c(11, 31)) & (encoding_id == 0) & (arrayed == 0)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

ukb_0code_0arr_info <- ukb_val_info |>
  filter(field_id %in% con_0code_0arr_vec)
#--------------------------------------------------------------------------
# 1.1 Raw and exported variables are both continuous.
# Manual inspection showed that most variables need no change except log-transforming blood and urine measurements.
#--------------------------------------------------------------------------
con_0code_0arr_cc <- ukb_0code_0arr_info |>
  filter(coding_note=="Continuous")

result_0code_0arr_cc <- data.frame(eid = df_ukb_common$eid)
  
for (i in seq_len(nrow(con_0code_0arr_cc))) {
  cat(i,"/",nrow(con_0code_0arr_cc), "begin at",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")
  data_type <- con_0code_0arr_cc$TopLevel[i]
  old_col <- con_0code_0arr_cc$ukb_colname[i]
  new_col <- con_0code_0arr_cc$exp_code[i]
  
  if (old_col %in% colnames(df_ukb_common)) {
    if (data_type == "Biological samples") {
      result_0code_0arr_cc[[new_col]] <- lognmin(df_ukb_common[[old_col]])
    } else {
      result_0code_0arr_cc[[new_col]] <- df_ukb_common[[old_col]]
    }
  } else {
    warning(paste("Column", old_col, "does not exist in the data frame"))
    result_0code_0arr_cc[[new_col]] <- NA
  }
}

result_0code_0arr_cc <- result_0code_0arr_cc |>
  arrange(eid)

saveRDS(result_0code_0arr_cc, file=here("tmp_files/ukb_data_extract/con_0code_0arr_cc.RDS"))
rm(con_0code_0arr_cc,result_0code_0arr_cc,i,data_type,old_col,new_col)

#--------------------------------------------------------------------------
# 1.2 Raw variables are continuous but exported variables are binary.
#--------------------------------------------------------------------------
con_0code_0arr_cb <- ukb_0code_0arr_info |>
  filter(coding_note=="Binary")

code_0code_0arr_cb <- ukb_val_code |>
  filter(unique_id %in% con_0code_0arr_cb$unique_id)

## Split binary variables by sign automatically and handle the remaining variables manually.
tmp_vec <- code_0code_0arr_cb |>
  filter(coding_group %in% c(">0","=0")) |>
  pull(unique_id) |>
  unique()

## Continuous variables converted to binary classes using zero as the cutoff.
con_0code_0arr_cb1 <- con_0code_0arr_cb |>
  filter(unique_id %in% tmp_vec)

code_0code_0arr_cb1 <- code_0code_0arr_cb |>
  filter(unique_id %in% tmp_vec)

result_0code_0arr_cb1 <- data.frame(eid = df_ukb_common$eid)

for (i in seq_len(nrow(con_0code_0arr_cb1))) {
  cat("***",i,"/",nrow(con_0code_0arr_cb1), "begin at",format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  old_col <- con_0code_0arr_cb1$ukb_colname[i]
  new_col <- con_0code_0arr_cb1$exp_code[i]
  
  if (old_col %in% colnames(df_ukb_common)) {
    vec <- df_ukb_common[[old_col]]
    new_vec <- ifelse(
      is.na(vec),
      NA_character_,
      ifelse(vec > 0, "g1_gt0",
            ifelse(vec == 0, "g0_eq0", NA_character_))
    )

  print(table(new_vec))

  result_0code_0arr_cb1[[new_col]] <- factor(
    new_vec,
    levels = c("g0_eq0", "g1_gt0")
  )
  print(summary(result_0code_0arr_cb1[[new_col]]))
  cat("===========================\n")
  } else {
    warning(paste("Column", old_col, "does not exist in the data frame"))
  }
}

result_0code_0arr_cb1 <- result_0code_0arr_cb1 |>
  arrange(eid)

saveRDS(result_0code_0arr_cb1, file=here("tmp_files/ukb_data_extract/con_0code_0arr_cb1.RDS"))
rm(con_0code_0arr_cb1,code_0code_0arr_cb1,result_0code_0arr_cb1,i,old_col,new_col,vec,new_vec)

## Continuous variables exported as binary variables that require manual handling.
con_0code_0arr_cb2 <- con_0code_0arr_cb |>
  filter(!(unique_id %in% tmp_vec))

code_0code_0arr_cb2 <- code_0code_0arr_cb |>
  filter(!(unique_id %in% tmp_vec))

result_0code_0arr_cb2 <- data.frame(eid = df_ukb_common$eid)

### 24009
old_col <- con_0code_0arr_cb2$ukb_colname[1]
new_col <- con_0code_0arr_cb2$exp_code[1]
cat("---",old_col,"-",new_col,'---\n')
vec <- df_ukb_common[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse(vec > 500, "g1_gt500",
        ifelse(vec <= 500, "g0_le500", NA_character_))
)
print(table(new_vec))
result_0code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_le500", "g1_gt500")
)
print(summary(result_0code_0arr_cb2[[new_col]]))

### 20320 20321
for (i in 2:3) {
  old_col <- con_0code_0arr_cb2$ukb_colname[i]
  new_col <- con_0code_0arr_cb2$exp_code[i]
  cat("---",old_col,"-",new_col,'---\n')

  vec <- df_ukb_common[[old_col]]
  new_vec <- ifelse(
    is.na(vec),
    NA_character_,
    ifelse(vec > 0, "g1_gt0",
          ifelse(vec == 0, "g0_eq0", NA_character_))
  )

print(table(new_vec))

result_0code_0arr_cb2[[new_col]] <- factor(
  new_vec,
  levels = c("g0_eq0", "g1_gt0")
)
print(summary(result_0code_0arr_cb2[[new_col]]))
cat("------------------\n")
}

### 22335
old_col <- con_0code_0arr_cb2$ukb_colname[4]
new_col <- con_0code_0arr_cb2$exp_code[4]
cat("---",old_col,"-",new_col,'---\n')
vec <- df_ukb_common[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse((vec < 0 | vec > 75), "g1_lt0_OR_gt75",
        ifelse((vec >= 0 & vec <= 75), "g0_ge0_AND_le75", NA_character_))
)
print(table(new_vec))
result_0code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_ge0_AND_le75", "g1_lt0_OR_gt75")
)
print(summary(result_0code_0arr_cb2[[new_col]]))

### 22336
old_col <- con_0code_0arr_cb2$ukb_colname[5]
new_col <- con_0code_0arr_cb2$exp_code[5]
cat("---",old_col,"-",new_col,'---\n')
vec <- df_ukb_common[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse((vec < -30 | vec > 90), "g1_ltn30_OR_gt90",
        ifelse((vec >= -30 & vec <= 90), "g0_gen30_AND_le90", NA_character_))
)
print(table(new_vec))
result_0code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_gen30_AND_le90", "g1_ltn30_OR_gt90")
)
print(summary(result_0code_0arr_cb2[[new_col]]))

### 22337
old_col <- con_0code_0arr_cb2$ukb_colname[6]
new_col <- con_0code_0arr_cb2$exp_code[6]
cat("---",old_col,"-",new_col,'---\n')
vec <- df_ukb_common[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse((vec < 15 | vec > 75), "g1_lt15_OR_gt75",
        ifelse((vec >= 15 & vec <= 75), "g0_ge15_AND_le75", NA_character_))
)
print(table(new_vec))
result_0code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_ge15_AND_le75", "g1_lt15_OR_gt75")
)
print(summary(result_0code_0arr_cb2[[new_col]]))

result_0code_0arr_cb2 <- result_0code_0arr_cb2 |>
  arrange(eid)

saveRDS(result_0code_0arr_cb2, file=here("tmp_files/ukb_data_extract/con_0code_0arr_cb2.RDS"))
rm(con_0code_0arr_cb2,code_0code_0arr_cb2,result_0code_0arr_cb2,old_col,new_col,vec,new_vec)
rm(con_0code_0arr_cb,code_0code_0arr_cb,tmp_vec)

#--------------------------------------------------------------------------
# 1.3 Raw variables are continuous but exported variables are ordinal categorical variables.
#--------------------------------------------------------------------------
con_0code_0arr_co <- ukb_0code_0arr_info |>
  filter(coding_note=="Ordinal categorical")

code_0code_0arr_co <- ukb_val_code |>
  filter(unique_id %in% con_0code_0arr_co$unique_id)

## As above, inspect variables whose labels start with t1 or numbers.
tmp_vec <- code_0code_0arr_co |>
  filter(grepl("^t1", coding_group)) |>
  pull(unique_id) |>
  unique()

## Continuous variables converted to three categories by calculating cut points.
con_0code_0arr_co1 <- con_0code_0arr_co |>
  filter(unique_id %in% tmp_vec)

code_0code_0arr_co1 <- code_0code_0arr_co |>
  filter(unique_id %in% tmp_vec)

result_0code_0arr_co1 <- data.frame(eid = df_ukb_common$eid)

# These variables are mostly categorized by tertiles.
sink(here("tmp_files/output_log.txt"), append = TRUE, split = TRUE)
for (i in seq_len(nrow(con_0code_0arr_co1))) {
  cat("===",i,"/",nrow(con_0code_0arr_co1), "begin at",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"===\n")
  old_col <- con_0code_0arr_co1$ukb_colname[i]
  new_col <- con_0code_0arr_co1$exp_code[i]
  cat("||--",old_col,"-",new_col,'---\n')

  vec <- df_ukb_common[[old_col]]

  quantiles <- quantile(
    vec, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
  )

  min_val <- quantiles[1] 
  q33 <- quantiles[2]
  q67 <- quantiles[3]
  max_val <- quantiles[4]

  factor_levels <- c(
    paste0("t1_", as.integer(min_val), "_", as.integer(q33)),
    paste0("t2_", as.integer(q33), "_", as.integer(q67)),
    paste0("t3_", as.integer(q67), "_", as.integer(max_val))
  )

  new_vec <- dplyr::case_when(
    is.na(vec)              ~ NA_character_,
    vec <= q33              ~ factor_levels[1],
    vec > q33 & vec <= q67  ~ factor_levels[2],
    vec > q67               ~ factor_levels[3],
    TRUE                    ~ NA_character_
  )

  print(table(new_vec))

  prop_vec <- prop.table(table(new_vec))
  if (any(prop_vec > 0.40) || any(prop_vec < 0.25)) {
    cat("**Note Percent**:", paste0(names(prop_vec), "=", round(prop_vec * 100, 1), "%", collapse = ", "), "\n")
  }

  new_vec <- factor(
    new_vec,
    levels  = factor_levels,
    ordered = TRUE
  )

  result_0code_0arr_co1[[new_col]] <- new_vec

  print(summary(result_0code_0arr_co1[[new_col]]))

 cat("======================================================\n")
}
sink()

result_0code_0arr_co1 <- result_0code_0arr_co1 |>
  arrange(eid)

saveRDS(result_0code_0arr_co1, file=here("tmp_files/ukb_data_extract/con_0code_0arr_co1.RDS"))
rm(con_0code_0arr_co1,code_0code_0arr_co1,result_0code_0arr_co1,i,old_col,new_col,vec,new_vec)

## Continuous variables exported as ordinal categories that require manual handling.
con_0code_0arr_co2 <- con_0code_0arr_co |>
  filter(!(unique_id %in% tmp_vec))

code_0code_0arr_co2 <- code_0code_0arr_co |>
  filter(!(unique_id %in% tmp_vec))

result_0code_0arr_co2 <- data.frame(eid = df_ukb_common$eid)

### 26075
old_col <- con_0code_0arr_co2$ukb_colname[1]
new_col <- con_0code_0arr_co2$exp_code[1]
cat("---",old_col,"-",new_col,'---\n')
vec <- df_ukb_common[[old_col]]

new_vec <- dplyr::case_when(
  is.na(vec)   ~ NA_integer_,
  vec == 0     ~ 1,
  vec == 22    ~ 2,
  vec == 44    ~ 3,
  vec == 88    ~ 3,
  TRUE         ~ NA_integer_
)

print(table(new_vec))
result_0code_0arr_co2[[new_col]] <- factor(
    new_vec,
    levels = c(1, 2, 3),
    ordered = TRUE
)
print(summary(result_0code_0arr_co2[[new_col]]))

### 26105
old_col <- con_0code_0arr_co2$ukb_colname[2]
new_col <- con_0code_0arr_co2$exp_code[2]
cat("---",old_col,"-",new_col,'---\n')
vec <- df_ukb_common[[old_col]]

new_vec <- dplyr::case_when(
  is.na(vec)   ~ NA_integer_,
  vec == 0     ~ 1,
  vec == 40    ~ 2,
  vec == 80    ~ 3,
  vec == 160   ~ 3,
  TRUE         ~ NA_integer_
)

print(table(new_vec))
result_0code_0arr_co2[[new_col]] <- factor(
    new_vec,
    levels = c(1, 2, 3),
    ordered = TRUE
)
print(summary(result_0code_0arr_co2[[new_col]]))

result_0code_0arr_co2 <- result_0code_0arr_co2 |>
  arrange(eid)

saveRDS(result_0code_0arr_co2, file=here("tmp_files/ukb_data_extract/con_0code_0arr_co2.RDS"))
rm(con_0code_0arr_co2,code_0code_0arr_co2,result_0code_0arr_co2,old_col,new_col,vec,new_vec)

rm(con_0code_0arr_co,code_0code_0arr_co,tmp_vec,ukb_0code_0arr_info)

#######################################################################
# 2. Raw data are continuous or integer, coding is not 0, and fields are not arrays.
#######################################################################

con_1code_0arr_vec <- ukb_sum |>
  filter((value_type %in% c(11, 31)) & (encoding_id != 0) & (arrayed == 0)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

ukb_1code_0arr_info <- ukb_val_info |>
  filter(field_id %in% con_1code_0arr_vec)

#--------------------------------------------------------------------------
# 2.1 Raw variables are continuous with coded components, and exported variables are continuous.
#--------------------------------------------------------------------------
# The continuous subset is empty.
con_1code_0arr_cc <- ukb_1code_0arr_info |>
  filter(coding_note == "Continuous")

rm(con_1code_0arr_cc)

#--------------------------------------------------------------------------
# 2.2 Raw variables are continuous with coded components but exported variables are binary.
#--------------------------------------------------------------------------
# Binary variables.
### Note: coded variables are decoded during export and are therefore character variables.

con_1code_0arr_cb <- ukb_1code_0arr_info |>
  filter(coding_note == "Binary")

code_1code_0arr_cb <- ukb_val_code |>
  filter(unique_id %in% con_1code_0arr_cb$unique_id)

table(code_1code_0arr_cb$coding_label)

# Extract these columns first.
# Convert 'Do not know' and 'Prefer not to answer' to NA.
# Convert 'Less than once a year' to 0.5.
# Convert 'No pain or discomfort' to 0 and 'Pain with unspecified rating' to 1.

data_1code_0arr_cb <- df_ukb_common |>
  select(eid, all_of(unique(con_1code_0arr_cb$ukb_colname)))

data_1code_0arr_cb <- data_1code_0arr_cb %>%
  mutate(
    across(
      -1,
      ~ {
        x <- as.character(.)
        x[x == ""] <- NA_character_
        x[x %in% c("Do not know", "Prefer not to answer")] <- NA_character_
        x[x == "Less than once a year"] <- "0.5"
        x[x == "No pain or discomfort"] <- "0"
        x[x == "Pain with unspecified rating"] <- "1"
        as.numeric(x)
      }
    )
  )

# Automatically process variables separated by sign and manually handle the rest.
tmp_vec <- code_1code_0arr_cb |>
  filter(coding_group %in% c(">0","=0")) |>
  pull(unique_id) |>
  unique()

## Continuous variables converted to binary classes using zero as the cutoff.
con_1code_0arr_cb1 <- con_1code_0arr_cb |>
  filter(unique_id %in% tmp_vec)

code_1code_0arr_cb1 <- code_1code_0arr_cb |>
  filter(unique_id %in% tmp_vec)

result_1code_0arr_cb1 <- data.frame(eid = data_1code_0arr_cb$eid)

for (i in seq_len(nrow(con_1code_0arr_cb1))) {
  cat("***",i,"/",nrow(con_1code_0arr_cb1), "begin at",as.character(Sys.time()))
  old_col <- con_1code_0arr_cb1$ukb_colname[i]
  new_col <- con_1code_0arr_cb1$exp_code[i]
  
  if (old_col %in% colnames(data_1code_0arr_cb)) {
    vec <- data_1code_0arr_cb[[old_col]]
    new_vec <- ifelse(
      is.na(vec),
      NA_character_,
      ifelse(vec > 0, "g1_gt0",
            ifelse(vec == 0, "g0_eq0", NA_character_))
    )

  print(table(new_vec))

  result_1code_0arr_cb1[[new_col]] <- factor(
    new_vec,
    levels = c("g0_eq0", "g1_gt0")
  )
  print(summary(result_1code_0arr_cb1[[new_col]]))
  cat("===========================\n")
  } else {
    warning(paste("Column", old_col, "does not exist in the data frame"))
  }
}

result_1code_0arr_cb1 <- result_1code_0arr_cb1 |>
  arrange(eid)

saveRDS(result_1code_0arr_cb1, file=here("tmp_files/ukb_data_extract/con_1code_0arr_cb1.RDS"))
rm(con_1code_0arr_cb1,code_1code_0arr_cb1,result_1code_0arr_cb1,i,old_col,new_col,vec,new_vec)

## Continuous variables exported as binary variables that require manual handling.
con_1code_0arr_cb2 <- con_1code_0arr_cb |>
  filter(!(unique_id %in% tmp_vec))

code_1code_0arr_cb2 <- code_1code_0arr_cb |>
  filter(!(unique_id %in% tmp_vec))

result_1code_0arr_cb2 <- data.frame(eid = data_1code_0arr_cb$eid)

### 21044
old_col <- con_1code_0arr_cb2$ukb_colname[1]
new_col <- con_1code_0arr_cb2$exp_code[1]
cat("---",old_col,"-",new_col,'---\n')
vec <- data_1code_0arr_cb[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse(vec > 1, "g1_gt1",
        ifelse(vec <= 1, "g0_le1", NA_character_))
)
print(table(new_vec))
result_1code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_le1", "g1_gt1")
)
print(summary(result_1code_0arr_cb2[[new_col]]))

### 1160_1
old_col <- con_1code_0arr_cb2$ukb_colname[2]
new_col <- con_1code_0arr_cb2$exp_code[2]
cat("---",old_col,"-",new_col,'---\n')
vec <- data_1code_0arr_cb[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse(vec < 7, "g1_lt7",
         ifelse((vec >= 7 & vec <= 8), "g0_ge7_AND_le8", NA_character_))
)
print(table(new_vec))
result_1code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_ge7_AND_le8", "g1_lt7")
)
print(summary(result_1code_0arr_cb2[[new_col]]))

### 1160_2
old_col <- con_1code_0arr_cb2$ukb_colname[3]
new_col <- con_1code_0arr_cb2$exp_code[3]
cat("---",old_col,"-",new_col,'---\n')
vec <- data_1code_0arr_cb[[old_col]]
new_vec <- ifelse(
  is.na(vec),
  NA_character_,
  ifelse(vec > 8, "g1_gt8",
         ifelse((vec >= 7 & vec <= 8), "g0_ge7_AND_le8", NA_character_))
)
print(table(new_vec))
result_1code_0arr_cb2[[new_col]] <- factor(
    new_vec,
    levels = c("g0_ge7_AND_le8", "g1_gt8")
)
print(summary(result_1code_0arr_cb2[[new_col]]))

result_1code_0arr_cb2 <- result_1code_0arr_cb2 |>
  arrange(eid)

saveRDS(result_1code_0arr_cb2, file=here("tmp_files/ukb_data_extract/con_1code_0arr_cb2.RDS"))
rm(con_1code_0arr_cb2,code_1code_0arr_cb2,result_1code_0arr_cb2,old_col,new_col,vec,new_vec)
rm(con_1code_0arr_cb,code_1code_0arr_cb,tmp_vec)
rm(data_1code_0arr_cb)

#--------------------------------------------------------------------------
# 2.3 Raw variables are continuous but exported variables are ordinal categorical variables.
#--------------------------------------------------------------------------
con_1code_0arr_co <- ukb_1code_0arr_info |>
  filter(coding_note == "Ordinal categorical")

code_1code_0arr_co <- ukb_val_code |>
  filter(unique_id %in% con_1code_0arr_co$unique_id)

table(code_1code_0arr_co$coding_label)

# Extract these columns first.
# Convert 'Do not know' and 'Prefer not to answer' to NA.
# Convert 'Less than once a year' to 0.5.
# Convert 'No pain or discomfort' to 0 and 'Pain with unspecified rating' to 1.

data_1code_0arr_co <- df_ukb_common |>
  select(eid, all_of(unique(con_1code_0arr_co$ukb_colname)))

# Inspect the number of mappings because these cannot all be mapped manually.

label2num <- code_1code_0arr_co |>
  select(coding_label,coding_num) |>
  filter(coding_label!='-') |>
  mutate(label_num = paste(coding_label,coding_num,sep="::")) |>
  distinct(label_num, .keep_all = TRUE) |>
  mutate(coding_num = str_remove(coding_num, "^assign\\s+"))

label_map <- setNames(
  as.character(label2num$coding_num),
  label2num$coding_label
)

# ======== Diagnostic 1: check whether one label maps to multiple numeric values. ========
label_conflict <- label2num %>%
  group_by(coding_label) %>%
  summarise(
    n_unique_num = n_distinct(coding_num),
    unique_nums  = paste(unique(coding_num), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_unique_num > 1)
if (nrow(label_conflict) > 0) {
  cat("**Label Conflict Found**: The following labels map to multiple numeric values:\n")
  print(label_conflict)
} else {
  cat("**Label Conflict Check Passed**: No label mapping conflicts found\n")
}

data_1code_0arr_co <- data_1code_0arr_co %>%
  mutate(
    across(
      -1,
      ~ {
        x <- as.character(.x)
        x[x %in% names(label_map)] <- label_map[x[x %in% names(label_map)]]
        x[x %in% c("", "NA")] <- NA_character_
        x
      }
    )
  )

# ======== Diagnostic 2: check whether non-numeric unmapped values were silently converted to NA. ========
cols_to_check <- setdiff(names(data_1code_0arr_co), "eid")
unmapped_report <- list()
for (col_i in cols_to_check) {
  raw_vals <- unique(as.character(data_1code_0arr_co[[col_i]]))
  unmapped <- setdiff(raw_vals, c(names(label_map), "", "NA"))
  unmapped <- unmapped[!is.na(unmapped)]
  if (length(unmapped) > 0) {
    # Check whether all unmapped values are numeric strings.
    is_num <- grepl("^-?\\d+\\.?\\d*$", unmapped)
    non_num <- unmapped[!is_num]
    if (length(non_num) > 0) {
      unmapped_report[[col_i]] <- non_num
    }
  }
}
if (length(unmapped_report) > 0) {
  cat("**Unmapped Values Found**: The following columns contain non-numeric unmapped values that will be converted to NA by as.numeric():\n")
  for (nm in names(unmapped_report)) {
    cat("  ", nm, ": ", paste(unmapped_report[[nm]], collapse = ", "), "\n")
  }
} else {
  cat("**Unmapped Check Passed**: All values are either in label_map or have been handled correctly\n")
}

data_1code_0arr_co <- data_1code_0arr_co %>%
  mutate(across(-1, as.numeric))

## As above, inspect variables whose labels start with t1 or numbers.
tmp_vec <- code_1code_0arr_co |>
  filter(grepl("^t1", coding_group)) |>
  pull(unique_id) |>
  unique()

## Continuous variables converted to three categories by calculating cut points.
con_1code_0arr_co1 <- con_1code_0arr_co |>
  filter(unique_id %in% tmp_vec)

code_1code_0arr_co1 <- code_1code_0arr_co |>
  filter(unique_id %in% tmp_vec)

result_1code_0arr_co1 <- data.frame(eid = data_1code_0arr_co$eid)

# These variables are mostly categorized by tertiles.
sink(here("tmp_files/result_1code_0arr_co1_log.txt"), split = TRUE)
for (i in seq_len(nrow(con_1code_0arr_co1))) {
  cat("===",i,"/",nrow(con_1code_0arr_co1), "begin at",as.character(Sys.time()),"===\n")
  old_col <- con_1code_0arr_co1$ukb_colname[i]
  new_col <- con_1code_0arr_co1$exp_code[i]
  cat("||--",old_col,"-",new_col,'---\n')

  vec <- data_1code_0arr_co[[old_col]]

  quantiles <- quantile(
    vec, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
  )

  min_val <- quantiles[1]
  q33 <- quantiles[2]
  q67 <- quantiles[3]
  max_val <- quantiles[4]

  factor_levels <- c(
    paste0("t1_", as.integer(min_val), "_", as.integer(q33)),
    paste0("t2_", as.integer(q33), "_", as.integer(q67)),
    paste0("t3_", as.integer(q67), "_", as.integer(max_val))
  )

  new_vec <- dplyr::case_when(
    is.na(vec)              ~ NA_character_,
    vec <= q33              ~ factor_levels[1],
    vec > q33 & vec <= q67  ~ factor_levels[2],
    vec > q67               ~ factor_levels[3],
    TRUE                    ~ NA_character_
  )

  print(table(new_vec))
  prop_vec <- prop.table(table(new_vec))
  if (any(prop_vec > 0.40) || any(prop_vec < 0.25)) {
    cat("**Note Percent**:", paste0(names(prop_vec), "=", round(prop_vec * 100, 1), "%", collapse = ", "), "\n")
  }

  new_vec <- factor(
    new_vec,
    levels  = factor_levels,
    ordered = TRUE
  )

  result_1code_0arr_co1[[new_col]] <- new_vec

  print(summary(result_1code_0arr_co1[[new_col]]))
 
 cat("======================================================\n")
}
sink()

# Manually adjust three-category variables with problematic groupings.
# p2734_i0 - tidy_p2734
old_col <- 'p2734_i0'
new_col <- 'tidy_p2734'
cat("---",old_col,"-",new_col,'---\n')
vec <- data_1code_0arr_co[[old_col]]

new_vec <- dplyr::case_when(
  is.na(vec)           ~ NA_character_,
  vec == 0             ~ 't1_0_0',
  (vec > 0 & vec <= 2) ~ 't2_0_2',
  vec > 2              ~ 't3_2_22',
  TRUE                 ~ NA_character_
)

print(table(new_vec))
result_1code_0arr_co1[[new_col]] <- factor(
    new_vec,
    levels = c('t1_0_0', 't2_0_2', 't3_2_22'),
    ordered = TRUE
)
print(summary(result_1code_0arr_co1[[new_col]]))

# p2744_i0 - tidy_p2744
old_col <- 'p2744_i0'
new_col <- 'tidy_p2744'
cat("---",old_col,"-",new_col,'---\n')
vec <- data_1code_0arr_co[[old_col]]

new_vec <- dplyr::case_when(
  is.na(vec)           ~ NA_character_,
  vec <= 5             ~ 't1_1_5',
  (vec > 5 & vec <= 7) ~ 't2_5_7',
  vec > 7              ~ 't3_7_16',
  TRUE                 ~ NA_character_
)

print(table(new_vec))
result_1code_0arr_co1[[new_col]] <- factor(
    new_vec,
    levels = c('t1_1_5', 't2_5_7', 't3_7_16'),
    ordered = TRUE
)
print(summary(result_1code_0arr_co1[[new_col]]))

# -- p864_i0 - tidy_p864 ---
old_col <- 'p864_i0'
new_col <- 'tidy_p864'
cat("---",old_col,"-",new_col,'---\n')
vec <- data_1code_0arr_co[[old_col]]

new_vec <- dplyr::case_when(
  is.na(vec)           ~ NA_character_,
  vec <= 2             ~ 't1_0_2',
  (vec > 2 & vec <= 5) ~ 't2_2_5',
  vec > 5              ~ 't3_5_7',
  TRUE                 ~ NA_character_
)

print(table(new_vec))
result_1code_0arr_co1[[new_col]] <- factor(
    new_vec,
    levels = c('t1_0_2', 't2_2_5', 't3_5_7'),
    ordered = TRUE
)
print(summary(result_1code_0arr_co1[[new_col]]))

result_1code_0arr_co1 <- result_1code_0arr_co1 |>
  arrange(eid)

saveRDS(result_1code_0arr_co1, file=here("tmp_files/ukb_data_extract/con_1code_0arr_co1.RDS"))
rm(con_1code_0arr_co1,code_1code_0arr_co1,result_1code_0arr_co1,i,old_col,new_col,vec,new_vec)

## Continuous variables exported as ordinal categories that require manual handling.
con_1code_0arr_co2 <- con_1code_0arr_co |>
  filter(!(unique_id %in% tmp_vec))

code_1code_0arr_co2 <- code_1code_0arr_co |>
  filter(!(unique_id %in% tmp_vec))

result_1code_0arr_co2 <- data.frame(eid = df_ukb_common$eid)

sink(here("tmp_files/result_1code_0arr_co2_log.txt"), split = TRUE)
for (i in seq_len(nrow(con_1code_0arr_co2))) {
  cat("===",i,"/",nrow(con_1code_0arr_co2), "begin at",as.character(Sys.time()),"===\n")
  old_col <- con_1code_0arr_co2$ukb_colname[i]
  new_col <- con_1code_0arr_co2$exp_code[i]
  cat("||--",old_col,"-",new_col,'---\n')

  vec <- data_1code_0arr_co[[old_col]]

  max_val <- max(vec,na.rm = TRUE)

  new_vec <- dplyr::case_when(
    is.na(vec)           ~ NA_character_,
    vec == 0             ~ 'q1_0_0',
    (vec > 0 & vec <= 3) ~ 'q2_0_3',
    (vec > 3 & vec <= 6) ~ 'q3_3_6',
    vec > 6              ~ paste0("q4_6_", as.integer(max_val)),
    TRUE                 ~ NA_character_
  )

  print(table(new_vec))

  new_vec <- factor(
    new_vec,
    levels  = c('q1_0_0','q2_0_3','q3_3_6',paste0("q4_6_", as.integer(max_val))),
    ordered = TRUE
  )

  result_1code_0arr_co2[[new_col]] <- new_vec

  print(summary(result_1code_0arr_co2[[new_col]]))
 
 cat("======================================================\n")
}
sink()

saveRDS(result_1code_0arr_co2, file=here("tmp_files/ukb_data_extract/con_1code_0arr_co2.RDS"))
rm(con_1code_0arr_co2,code_1code_0arr_co2,result_1code_0arr_co2,i,old_col,new_col,vec,new_vec,min_val,max_val,q33,q67,quantiles,factor_levels)
rm(con_1code_0arr_co,code_1code_0arr_co,tmp_vec,ukb_1code_0arr_info,con_1code_0arr_vec,data_1code_0arr_co,label2num)

#######################################################################
# 3. Raw data are continuous or integer, coding is 0, and fields are arrays.
#######################################################################

con_0code_1arr_vec <- ukb_sum |>
  filter((value_type %in% c(11, 31)) & (encoding_id == 0) & (arrayed != 0)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

ukb_0code_1arr_info <- ukb_val_info |>
  filter(field_id %in% con_0code_1arr_vec)

con_0code_1arr <- ukb_0code_1arr_info
code_0code_1arr <- ukb_val_code |>
  filter(unique_id %in% con_0code_1arr$unique_id)

# Calculate means first because these variables are processed after averaging.

tmp_0code_1arr <- data.frame(eid = df_ukb_common$eid)
for (i in seq_len(nrow(con_0code_1arr))) {
  cat("===",i,"/",nrow(con_0code_1arr), "begin at",as.character(Sys.time()),"===\n")
  old_col_list_str <- con_0code_1arr$ukb_colname[i]
  old_col_list <- strsplit(old_col_list_str, "\\|")[[1]]
  new_col <- paste0('tmp_',con_0code_1arr$exp_code[i])
  cat("||--",old_col_list_str,"-",new_col,'---\n')

  df_tmp <- df_ukb_common[, ..old_col_list]

  new_vec <- rowMeans(df_tmp, na.rm = TRUE)

  tmp_0code_1arr[[new_col]] <- new_vec
 cat("======================================================\n")
}

result_0code_1arr <- data.frame(eid = df_ukb_common$eid)
# Automatically classify the two ordinal variables.

ordinal_exp_code_vec <- c('tidy_p399','tidy_p3064')
for (i in 1:2) {
  new_col <- ordinal_exp_code_vec[i]
  old_col <- paste0('tmp_',new_col)
  cat("||--",old_col,"-",new_col,'---\n')

  vec <- tmp_0code_1arr[[old_col]]

  quantiles <- quantile(
    vec, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
  )

  min_val <- quantiles[1] 
  q33 <- quantiles[2]
  q67 <- quantiles[3]
  max_val <- quantiles[4]

  factor_levels <- c(
    paste0("t1_", as.integer(min_val), "_", as.integer(q33)),
    paste0("t2_", as.integer(q33), "_", as.integer(q67)),
    paste0("t3_", as.integer(q67), "_", as.integer(max_val))
  )

  new_vec <- dplyr::case_when(
    is.na(vec)              ~ NA_character_,
    vec <= q33              ~ factor_levels[1],
    vec > q33 & vec <= q67  ~ factor_levels[2],
    vec > q67               ~ factor_levels[3],
    TRUE                    ~ NA_character_
  )

  print(table(new_vec))

  new_vec <- factor(
    new_vec,
    levels  = factor_levels,
    ordered = TRUE
  )

  result_0code_1arr[[new_col]] <- new_vec
  print(summary(result_0code_1arr[[new_col]]))
  cat("======================================================\n")
}

# Manually classify the remaining five binary variables in two groups.
# c('tidy_p4079','tidy_p102','tidy_p4080')
binary_ukb_col_vec <- c('tidy_p4079','tidy_p102','tidy_p4080')
check_point_vec <- c(80,100,120)

for (i in 1:3) {
  new_col <- binary_ukb_col_vec[i]
  old_col <- paste0('tmp_',new_col)
  cat("***",i,"/ 3",new_col, "begin at",as.character(Sys.time()))

  check_point <- check_point_vec[i]
  
  vec <- tmp_0code_1arr[[old_col]]
  new_vec <- ifelse(
    is.na(vec),
    NA_character_,
    ifelse(vec > check_point, paste0("g1_gt",check_point),
          ifelse(vec <= check_point, paste0("g0_le",check_point), NA_character_))
  )
  print(table(new_vec))
  result_0code_1arr[[new_col]] <- factor(
    new_vec,
    levels = paste0(c("g0_le","g1_gt"),check_point)
  )
  print(summary(result_0code_1arr[[new_col]]))
  cat("===========================\n")
}

# c('tidy_p5086','tidy_p5087')
binary_ukb_col_vec <- c('tidy_p5086','tidy_p5087')

for (i in 1:2) {
  new_col <- binary_ukb_col_vec[i]
  old_col <- paste0('tmp_',new_col)
  cat("***",i,"/ 2",new_col, "begin at",as.character(Sys.time()))
  
  vec <- tmp_0code_1arr[[old_col]]
  new_vec <- ifelse(
    is.na(vec),
    NA_character_,
    ifelse(vec > 0.5, "g1_gt05",
          ifelse((vec >= 0 & vec <= 0.5), "g0_0_05", NA_character_))
  )
  print(table(new_vec))
  result_0code_1arr[[new_col]] <- factor(
    new_vec,
    levels = c("g0_0_05","g1_gt05")
  )
  print(summary(result_0code_1arr[[new_col]]))
  cat("===========================\n")
}

result_0code_1arr <- result_0code_1arr |>
  arrange(eid)

saveRDS(result_0code_1arr, file=here("tmp_files/ukb_data_extract/con_0code_1arr.RDS"))
rm(code_0code_1arr,con_0code_1arr,df_tmp,result_0code_1arr,tmp_0code_1arr,ukb_0code_1arr_info,
   binary_ukb_col_vec,check_point,check_point_vec,con_0code_1arr_vec,factor_levels,i,max_val,min_val,
   new_col,new_vec,old_col,old_col_list,old_col_list_str,ordinal_exp_code_vec,ordinal_ukb_col_vec,q33,q67,quantiles,vec)

#######################################################################
# 4. Raw data are continuous or integer, coding is not 0, and fields are arrays.
#######################################################################

con_1code_1arr_vec <- ukb_sum |>
  filter((value_type %in% c(11, 31)) & (encoding_id != 0) & (arrayed != 0)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

ukb_1code_1arr_info <- ukb_val_info |>
  filter(field_id %in% con_1code_1arr_vec)

con_1code_1arr <- ukb_1code_1arr_info
code_1code_1arr <- ukb_val_code |>
  filter(unique_id %in% con_1code_1arr$unique_id)

# No variables fall into this group; continuous-variable processing is complete.

#######################################################################
# 5. Merge data and check whether expected column names are present.
#######################################################################

file_list <- c(
  'tmp_files/ukb_data_extract/con_0code_0arr_cc.RDS',
  'tmp_files/ukb_data_extract/con_0code_0arr_cb1.RDS',
  'tmp_files/ukb_data_extract/con_0code_0arr_cb2.RDS',
  'tmp_files/ukb_data_extract/con_0code_0arr_co1.RDS',
  'tmp_files/ukb_data_extract/con_0code_0arr_co2.RDS',
  'tmp_files/ukb_data_extract/con_1code_0arr_cb1.RDS',
  'tmp_files/ukb_data_extract/con_1code_0arr_cb2.RDS',
  'tmp_files/ukb_data_extract/con_1code_0arr_co1.RDS',
  'tmp_files/ukb_data_extract/con_1code_0arr_co2.RDS',
  'tmp_files/ukb_data_extract/con_0code_1arr.RDS'
)

file_list <- here(file_list)

data_list <- lapply(file_list, readRDS)

merged_data <- data_list %>%
  reduce(inner_join, by = "eid")

ukb_val_info_conti <- ukb_val_info |>
  filter(value_type %in% c(11,31))

table(ukb_val_info$value_type,useNA = "always")

not_in_merged <- setdiff(ukb_val_info_conti$exp_code, names(merged_data))
not_in_tidy <- setdiff(names(merged_data), ukb_val_info_conti$exp_code)

cat("Present in exp_code but missing from merged_data:", length(not_in_merged), "\n")
print(not_in_merged)

cat("\nPresent in merged_data but missing from exp_code:", length(not_in_tidy), "\n")
print(not_in_tidy)

# Check completed.
# Data extraction for this step is complete.

saveRDS(merged_data, file=here("processed_data/ukb_common_continus.RDS"))

################################################################################
# Special continuous variables
################################################################################

# Read input data.
ukb_val_info <- read_xlsx(
  here("rawdata/ukb_exposome.xlsx"),
  col_types = c(
    "field_id" = "text"
  )
) |>
  filter(grepl("\\|", field_id))

ukb_val_code <- read_xlsx(here("rawdata/ukb_codebook_0119.xlsx"), col_types = c("field_id" = "text")) |>
  filter(grepl("\\|", field_id))

ukb_val_code_info <- ukb_val_code |>
  drop_na(title) 
table(duplicated(ukb_val_code_info$unique_id))

ukb_sum <- read.csv(here("rawdata/ukb_data_summary.csv"))

df_ukb_common <- fread(here('processed_data/ukb_common.tsv'))

# Extract all field_ids that may be involved in this step.

field_id_table <- data.frame(
  ukb_colname = names(df_ukb_common)[-1],
  field_id = str_extract(names(df_ukb_common)[-1], "(?<=p)\\d+"),
  stringsAsFactors = FALSE
) |>
  group_by(field_id) |>
  summarise(
    ukb_colname = paste(str_sort(ukb_colname, numeric = TRUE), collapse = "|"),
    .groups = "drop"
) |>
  mutate(field_id = as.numeric(field_id))

field_id_vec <- ukb_val_info$field_id |>
  str_split("\\|") |>
  unlist() |>
  unique() |>
  as.integer()

field_id_table <- field_id_table |>
  filter(field_id %in% field_id_vec)

df_ukb_common <- df_ukb_common |>
  select(eid,all_of(field_id_table$ukb_colname))

head(df_ukb_common)

# Iterate over each row in ukb_val_info.
for (i in seq_along(ukb_val_info)) {
  
  cat('===',i,' BEGIN ===\n')

  # Get the field_id vector separated by |.
  field_id_vec <- str_split(ukb_val_info$field_id[i], "\\|")[[1]]
  
  # Generate the new column name.
  new_colname <- paste0('tidy_p', ukb_val_info$unique_id[i])
  
  # Initialize the temporary-column list.
  temp_cols <- c()
  
  # Iterate over each field_id.
  for (id in field_id_vec) {
    
    # Get the corresponding UKB column name.
    old_colname <- field_id_table %>%
      filter(field_id == id) %>%
      pull(ukb_colname)
    
    # Continue when a matching column name is found.
    if (length(old_colname) > 0 && !is.na(old_colname)) {
      
      # Create a temporary column name.
      tmp_colname <- paste0("tmp_", old_colname)
      
      # Calculate tertiles and create groups while ignoring NA.
      df_ukb_common[[tmp_colname]] <- cut(
        df_ukb_common[[old_colname]],
        breaks = quantile(df_ukb_common[[old_colname]], 
                         probs = c(0, 1/3, 2/3, 1), 
                         na.rm = TRUE),
        labels = c("t1", "t2", "t3"),
        include.lowest = TRUE
      )
      
      # Store the temporary column name.
      temp_cols <- c(temp_cols, tmp_colname)
    }
  }
  
  # Merge all temporary columns.
  if (length(temp_cols) > 0) {
    
    # Check whether any participant has values in multiple columns.
    # Select columns using data.table syntax.
    temp_data <- df_ukb_common[, ..temp_cols]
    non_na_counts <- rowSums(!is.na(temp_data))
    conflict_rows <- which(non_na_counts > 1)
    
    if (length(conflict_rows) > 0) {
      warning(paste0("Column '", new_colname, "' has ", length(conflict_rows), 
                     " participants with values in multiple field_ids. Using the first non-NA value."))
    }
    
    # Merge by taking the first non-NA value.
    df_ukb_common[[new_colname]] <- do.call(coalesce, as.list(temp_data))
    
    # Convert the new column to a factor.
    df_ukb_common[[new_colname]] <- factor(
      df_ukb_common[[new_colname]], 
      levels = c("t1", "t2", "t3"),ordered = TRUE
    )
    
    # Remove temporary columns.
    df_ukb_common[, (temp_cols) := NULL]
  }
  cat('===',i,' END ===\n')
}

result_df <- df_ukb_common |>
  select(eid,starts_with("tidy_p")) |>
  arrange(eid)

head(result_df)

saveRDS(result_df, file=here("processed_data/ukb_common_continus_special.RDS"))

################################################################################
# Common categorical variables
################################################################################

# Read input data.
ukb_val_info <- read_xlsx(
  here("rawdata/ukb_exposome.xlsx"),
  col_types = c(
    "field_id" = "text"
  )
) |>
  filter(exp_cate != "Proteome") |> # Exclude proteomics.
  filter((exp_subcate!='Brain MRI' | is.na(exp_subcate))) |>  # Exclude brain MRI.
  filter(Level2 != "Derived accelerometry" | is.na(Level2)) |># Exclude accelerometry, which has already been processed.
  filter(exp_cate != "Metabolome") |> # Exclude metabolomics.
  filter(!grepl("\\|", field_id)) |>  # Special cases handled manually.
  mutate(field_id = as.numeric(field_id))

ukb_val_code <- read_xlsx(here("rawdata/ukb_codebook_0119.xlsx"), col_types = c("field_id" = "text")) |>
  filter(!grepl("\\|", field_id)) |>
  mutate(field_id = as.numeric(field_id))

ukb_val_code_info <- ukb_val_code |>
  drop_na(title) 
table(duplicated(ukb_val_code_info$unique_id))

ukb_sum <- read.csv(here("rawdata/ukb_data_summary.csv"))

df_ukb_common <- fread(here("processed_data/ukb_common.tsv"))

field_id_table <- data.frame(
  ukb_colname = names(df_ukb_common)[-1],
  field_id = str_extract(names(df_ukb_common)[-1], "(?<=p)\\d+"),
  stringsAsFactors = FALSE
) |>
  group_by(field_id) |>
  summarise(
    ukb_colname = paste(str_sort(ukb_colname, numeric = TRUE), collapse = "|"),
    .groups = "drop"
) |>
  mutate(field_id = as.numeric(field_id)) # This ordering is important for extracting work-related variables later.

ukb_val_info <- ukb_val_info |>
  left_join(field_id_table, by='field_id')

# Create detailed type groups to make manual checks easier.

dir.create(here("tmp_files/ukb_data_extract"),
           recursive = TRUE, showWarnings = FALSE)

#######################################################################
# 1. Raw data are categorical and fields are not arrays.
#######################################################################

cat_0arr_vec <- ukb_sum |>
  filter((value_type %in% c(21)) & (arrayed == 0)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

ukb_cat_0arr_info <- ukb_val_info |>
  filter(field_id %in% cat_0arr_vec)

cat_code_0arr <- ukb_val_code |>
  filter(unique_id %in% ukb_cat_0arr_info$unique_id)

# Processing first maps coding_label values onto coding_num values.
# For Binary coding_note values, check the mapping and then convert to factor.
# For Ordinal categorical coding_note values, check that there are more than two classes and then convert to ordered factors.

result_cat_0arr <- data.frame(eid = df_ukb_common$eid)
sink(here("tmp_files/cat_0array_log.txt"), split = TRUE)
for (i in seq_len(nrow(ukb_cat_0arr_info))) {
  i_id <- ukb_cat_0arr_info$unique_id[i]
  cat("===",i,"/",nrow(ukb_cat_0arr_info), "begin at",as.character(Sys.time()),"===\n")

  code_df_tmp <- cat_code_0arr |>
    filter(unique_id==i_id)

  data_class <- code_df_tmp$coding_note[1]
  number_class <- length(na.omit(unique(code_df_tmp$coding_num[code_df_tmp$coding_num != "NA"])))

  old_col <- ukb_cat_0arr_info$ukb_colname[i]
  new_col <- ukb_cat_0arr_info$exp_code[i]
  cat("||--",old_col,"-",new_col,'--\n')
  cat("||--",data_class,"-",number_class,'--\n')

  label_map <- setNames(
    as.character(code_df_tmp$coding_num),
    code_df_tmp$coding_label
  )

  vec <- df_ukb_common[[old_col]]

  vec2 <- as.character(vec)
  # This step preserves True and False labels so they do not become TRUE and FALSE.

  if (is.logical(vec)) {
    vec2[vec2 == "TRUE"] <- "True"
    vec2[vec2 == "FALSE"] <- "False"
  }

  idx <- vec2 %in% names(label_map)
  vec2[idx] <- label_map[vec2[idx]]
  vec2[vec2 %in% c("", "NA")] <- NA_character_

  print(table(vec2))

  vec2 <- as.numeric(vec2)

  result_cat_0arr[['tmp_col']] <- vec2

  # Process binary variables.
  if (data_class == 'Binary') {
    # Check whether number_class equals 2.
    if (number_class != 2) {
      cat("WARNING:", new_col, "- Binary variable should have 2 classes, but has", number_class, "\n")
    }
    unique_vals <- unique(vec2[!is.na(vec2)])
    if (!all(unique_vals %in% c(0, 1))) {
      cat("WARNING:", new_col, "- Binary variable contains values other than 0, 1, and NA\n")
    }
    result_cat_0arr[[new_col]] <- factor(result_cat_0arr[['tmp_col']], levels = c(0, 1))
  }
  
  # Process ordinal categorical variables.
  if (data_class == 'Ordinal categorical') {
    if (number_class <= 2) {
      cat("WARNING:", new_col, "- Ordinal categorical variable should have >2 classes, but has", number_class, "\n")
    }
    unique_levels <- sort(unique(vec2[!is.na(vec2)]))
    result_cat_0arr[[new_col]] <- factor(result_cat_0arr[['tmp_col']],
                                          levels = unique_levels,
                                          ordered = TRUE)
  }
  
  print(summary(result_cat_0arr[[new_col]]))
  cat('=======================================\n')

  # Remove temporary columns.
  result_cat_0arr[['tmp_col']] <- NULL
}
sink()

result_cat_0arr <- result_cat_0arr |>
  arrange(eid)

saveRDS(result_cat_0arr, file=here("tmp_files/ukb_data_extract/result_cat_0arr.RDS"))
rm(cat_0arr_vec,ukb_cat_0arr_info,cat_code_0arr,result_cat_0arr,i,i_id,code_df_tmp,data_class,number_class,
   old_col,new_col,label_map,vec,vec2,idx,unique_vals,non_na_vals,unique_levels)

#######################################################################
# 2. Raw data are categorical and fields are arrays.
#######################################################################

cat_1arr_vec <- ukb_sum |>
  filter((value_type %in% c(21)) & (arrayed != 0)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

ukb_cat_1arr_info <- ukb_val_info |>
  filter(field_id %in% cat_1arr_vec)

cat_code_1arr <- ukb_val_code |>
  filter(unique_id %in% ukb_cat_1arr_info$unique_id)

# These 27 variables are work-related, so reuse the previous work-processing logic.
# The codebook does not explicitly mark how work variables should be encoded.
# Download coding 497 from the UKB website.
# Also read work-related columns from the original file.

col_extract_vec <- unique(unlist(strsplit(ukb_cat_1arr_info$ukb_colname, "|", fixed = TRUE)))

ukb_work_data <- df_ukb_common |>
  select(eid, all_of(col_extract_vec)) # The columns have already been ordered above.

code_df <- fread(here("rawdata/coding497.tsv"))
correspond_code <- get_top_level_mapping(code_df)

job_mapping <- correspond_code %>%
  select(work_content = 1, work_category = 2, unique_sig=3) %>%
  distinct() |>
  arrange(work_category)

work_categories <- sort(unique(job_mapping$work_category))

################
# Later checks identified two additional columns that need to be extracted.
# p22599 Number of jobs held
# p22661 Number of gap periods
###############

ukb_work_extra_data <- fread(here("processed_data/ukb_work_extra.tsv"))
setcolorder(ukb_work_extra_data, c("eid",sort_ukb_cols(names(ukb_work_extra_data)[-1])))

ukb_work_data <- ukb_work_extra_data |>
  inner_join(ukb_work_data, by='eid') |>
  arrange(eid)

names(ukb_work_data)

##### Data tidying is complete here; the remaining code handles vectorization and parallelization. #####

result_category <- process_work_data_vectorized(
  data = ukb_work_data,
  cutoff_year_col = "p53_i0",
  mode = "work_category",
  job_mapping = job_mapping,
  work_categories = work_categories
)

# Note that the codes are changed here.
# Recode to 1, 2, 3, 4, and 5 according to alphabetical order.
# Align explicitly by alphabetical order of work_categories to avoid unstable dcast column order.
result_category <- result_category %>%
  select(eid, all_of(work_categories))
names(result_category)[-1] <- paste0('tidy_p22601_', seq_along(work_categories)) 

# result_category <- readRDS(here("tmp_files/ukb_data_extract/work_category.RDS"))
# These 15 columns require data-type conversion.
work_categories[9]
# The ninth variable is categorized into office-work tertiles; the others are binary.
cols_to_factor <- paste0("tidy_p22601_", setdiff(1:15, 9))
result_category <- result_category %>%
  mutate(across(
    all_of(cols_to_factor),
    ~ factor(ifelse(.x > 0, "g1_gt0", "g0_0"), levels = c("g0_0","g1_gt0"))
  ))

# Process tidy_p22601_9 separately.
old_col <- 'tidy_p22601_9'
new_col <- 'tidy_p22601_9'
cat("||--",old_col,"-",new_col,'---\n')
vec <- result_category[[old_col]]
quantiles <- quantile(
  vec, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
)
min_val <- quantiles[1]
q33 <- quantiles[2]
q67 <- quantiles[3]
max_val <- quantiles[4]
factor_levels <- c(
  paste0("t1_", as.integer(min_val), "_", as.integer(q33)),
  paste0("t2_", as.integer(q33), "_", as.integer(q67)),
  paste0("t3_", as.integer(q67), "_", as.integer(max_val))
)
new_vec <- dplyr::case_when(
  is.na(vec)              ~ NA_character_,
  vec <= q33              ~ factor_levels[1],
  vec > q33 & vec <= q67  ~ factor_levels[2],
  vec > q67               ~ factor_levels[3],
  TRUE                    ~ NA_character_
)
print(table(new_vec))
new_vec <- factor(
  new_vec,
  levels  = factor_levels,
  ordered = TRUE
)
result_category[[new_col]] <- new_vec

summary(result_category)
saveRDS(result_category, file=here("tmp_files/ukb_data_extract/work_category.RDS"))

# result_gap <- process_work_data_vectorized(
#   data = ukb_work_data,
#   cutoff_year_col = "p53_i0",
#   mode = "gap"
# )
# This field was not requested for extraction and has no clear variable name.
# The function can still extract it if needed.

# The third type requires weights during calculation, so map strings to weight values first.
work_exp_list1 <- c(22606:22615)

setDT(ukb_work_data)
for (num in work_exp_list1) {
  # Match column names such as p22606, p22606_a0, and p22606_a1, but not p226060.
  cat('Processing...',num,'\n')
  pattern <- paste0("^p", num, "(?![0-9])")
  matching_cols <- grep(pattern, names(ukb_work_data), perl = TRUE, value = TRUE)
  
  if (length(matching_cols) > 0) {
    new_cols <- paste0("n", matching_cols)
    
    for (i in seq_along(matching_cols)) {
      old_col <- matching_cols[i]
      new_col <- new_cols[i]
      
      # Use data.table syntax for speed.
      ukb_work_data[, (new_col) := fcase(
        get(old_col) %in% c("Do not know", "Rarely/never"), 0,
        get(old_col) == "Sometimes", 1,
        get(old_col) == "Often", 2,
        default = NA_real_
      )]
    }
  }
}

names(ukb_work_data)

result_exposure1 <- process_work_data_vectorized(
  data = ukb_work_data,
  cutoff_year_col = "p53_i0",
  mode = "work_exposure",
  job_mapping = job_mapping,
  work_categories = work_categories,
  exp_list = work_exp_list1
)

names(result_exposure1)
# result_exposure1 <- readRDS(here("tmp_files/ukb_data_extract/work_exposure1.RDS"))
# Variables 6 to 9 require three categories, while variables 10 to 15 require binary categories.

for (i in 22606:22609) {
  old_col <- paste0('tidy_p',i)
  new_col <- paste0('tidy_p',i)
  cat("||--",old_col,"-",new_col,'---\n')
  vec <- result_exposure1[[old_col]]
  quantiles <- quantile(
    vec, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE
  )
  min_val <- quantiles[1]
  q33 <- quantiles[2]
  q67 <- quantiles[3]
  max_val <- quantiles[4]
  factor_levels <- c(
    paste0("t1_", as.integer(min_val), "_", as.integer(q33)),
    paste0("t2_", as.integer(q33), "_", as.integer(q67)),
    paste0("t3_", as.integer(q67), "_", as.integer(max_val))
  )
  new_vec <- dplyr::case_when(
    is.na(vec)              ~ NA_character_,
    vec <= q33              ~ factor_levels[1],
    vec > q33 & vec <= q67  ~ factor_levels[2],
    vec > q67               ~ factor_levels[3],
    TRUE                    ~ NA_character_
  )
  print(table(new_vec))
  new_vec <- factor(
    new_vec,
    levels  = factor_levels,
    ordered = TRUE
  )
  result_exposure1[[new_col]] <- new_vec
}

cols_to_factor1 <- paste0("tidy_p", 22610:22615)
result_exposure1 <- result_exposure1 %>%
  mutate(across(
    all_of(cols_to_factor1),
    ~ factor(ifelse(.x > 0, "g1_gt0", "g0_0"), levels = c("g0_0","g1_gt0"))
  ))

summary(result_exposure1)
saveRDS(result_exposure1, file=here("tmp_files/ukb_data_extract/work_exposure1.RDS"))

#----------------------------
# 22616 and 22620.
#---------------------------
work_exp_list2 <- c(22616,22620)
setDT(ukb_work_data)
for (num in work_exp_list2) {
  # Match column names such as p22606, p22606_a0, and p22606_a1, but not p226060.
  cat('Processing...',num,'\n')
  pattern <- paste0("^p", num, "(?![0-9])")
  matching_cols <- grep(pattern, names(ukb_work_data), perl = TRUE, value = TRUE)
  
  if (length(matching_cols) > 0) {
    new_cols <- paste0("n", matching_cols)
    
    for (i in seq_along(matching_cols)) {
      old_col <- matching_cols[i]
      new_col <- new_cols[i]
      
      # Use data.table syntax for speed.
      ukb_work_data[, (new_col) := fcase(
        get(old_col) == "No", 0,
        get(old_col) == "Yes", 1,
        default = NA_real_
      )]
    }
  }
}

names(ukb_work_data)

result_exposure2 <- process_work_data_vectorized(
  data = ukb_work_data,
  cutoff_year_col = "p53_i0",
  mode = "work_exposure",
  job_mapping = job_mapping,
  work_categories = work_categories,
  exp_list = work_exp_list2
)

names(result_exposure2) # No renaming is needed here.

# result_exposure2 <- readRDS(here("tmp_files/ukb_data_extract/work_exposure2.RDS"))
cols_to_factor2 <- paste0("tidy_p", c(22616,22620))
result_exposure2 <- result_exposure2 %>%
  mutate(across(
    all_of(cols_to_factor2),
    ~ factor(ifelse(.x > 0, "g1_gt0", "g0_0"), levels = c("g0_0","g1_gt0"))
  ))

summary(result_exposure2)

saveRDS(result_exposure2, file=here("tmp_files/ukb_data_extract/work_exposure2.RDS"))

################################################
# Merge all results above into one dataset.
# As above, check whether each expected column is present.
###############################################

file_list <- c(
  'tmp_files/ukb_data_extract/result_cat_0arr.RDS',
  'tmp_files/ukb_data_extract/work_category.RDS',
  'tmp_files/ukb_data_extract/work_exposure1.RDS',
  'tmp_files/ukb_data_extract/work_exposure2.RDS'
)

file_list <- here(file_list)

data_list <- lapply(file_list, readRDS)

merged_data <- data_list %>%
  reduce(inner_join, by = "eid")

ukb_val_info_cat <- ukb_val_info |>
  filter(value_type %in% c(21))

table(ukb_val_info$value_type,useNA = "always")

not_in_merged <- setdiff(ukb_val_info_cat$exp_code, names(merged_data))
not_in_tidy <- setdiff(names(merged_data), ukb_val_info_cat$exp_code)

cat("Present in exp_code but missing from merged_data:", length(not_in_merged), "\n")
print(not_in_merged)

cat("\nPresent in merged_data but missing from exp_code:", length(not_in_tidy), "\n")
print(not_in_tidy) # Removed from the original file? tidy_p22601_8.

# Check completed.
# Data extraction for this step is complete.

merged_data <- merged_data |>
  arrange(eid)

saveRDS(merged_data, file=here("processed_data/ukb_common_cat.RDS"))

################################################################################
# Common multi-category variables
################################################################################

# Read input data.
ukb_val_info <- read_xlsx(
  here("rawdata/ukb_exposome.xlsx"),
  col_types = c(
    "field_id" = "text"
  )
) |>
  filter(exp_cate != "Proteome") |> # Exclude proteomics.
  filter((exp_subcate!='Brain MRI' | is.na(exp_subcate))) |>  # Exclude brain MRI.
  filter(Level2 != "Derived accelerometry" | is.na(Level2)) |># Exclude accelerometry, which has already been processed.
  filter(exp_cate != "Metabolome") |> # Exclude metabolomics.
  filter(!grepl("\\|", field_id)) |>  # Special cases handled manually.
  mutate(field_id = as.numeric(field_id))

ukb_val_code <- read_xlsx(here("rawdata/ukb_codebook_0119.xlsx"), col_types = c("field_id" = "text")) |>
  filter(!grepl("\\|", field_id)) |>
  mutate(field_id = as.numeric(field_id))

ukb_val_code_info <- ukb_val_code |>
  drop_na(title) 
table(duplicated(ukb_val_code_info$unique_id))

ukb_sum <- read.csv(here("rawdata/ukb_data_summary.csv"))

df_ukb_common <- fread(here("processed_data/ukb_common.tsv"))

# Merge extracted columns with the metadata table to recover original UKB column names.
# Use unique_id as the new tidy variable name.

field_id_table <- data.frame(
  ukb_colname = names(df_ukb_common)[-1],
  field_id = str_extract(names(df_ukb_common)[-1], "(?<=p)\\d+"),
  stringsAsFactors = FALSE
) |>
  group_by(field_id) |>
  summarise(
    ukb_colname = paste(str_sort(ukb_colname, numeric = TRUE), collapse = "|"),
    .groups = "drop"
) |>
  mutate(field_id = as.numeric(field_id)) 

ukb_val_info <- ukb_val_info |>
  left_join(field_id_table, by='field_id')

dir.create(here("tmp_files/ukb_data_extract"),
           recursive = TRUE, showWarnings = FALSE)

#######################################################################
# Multi-category data processing.
#######################################################################

multi_cat_vec <- ukb_sum |>
  filter(value_type %in% c(22)) |>
  filter(field_id %in% ukb_val_info$field_id) |>
  pull(field_id)

multi_cat_info <- ukb_val_info |>
  filter(field_id %in% multi_cat_vec)

multi_cat_code <- ukb_val_code |>
  filter(unique_id %in% multi_cat_info$unique_id)

# The data structure in this section is relatively special.
# Only values marked as the string 'NA' are treated as missing; unmarked values are coded as 0.
# table(multi_cat_info$coding_note) shows that all tidied variables are binary.
# Therefore, only the case and control definitions need to be specified.
# Except for several variables at the beginning, the remaining variables can be processed automatically.

multi_cat_code$coding_num[is.na(multi_cat_code$coding_num) & is.na(multi_cat_code$title)] <- "0"
multi_cat_code$coding_label[is.na(multi_cat_code$coding_label)] <- "-"

# Some multi-category rows require special handling.
col_extract_vec <- unique(unlist(strsplit(multi_cat_info$ukb_colname, "|", fixed = TRUE)))

ukb_multicat_data <- df_ukb_common |>
  select(eid, all_of(col_extract_vec)) 

ukb_multicat_extra_data <- fread(here("processed_data/ukb_multicat_extra.tsv"))
ukb_multicat_data <- ukb_multicat_extra_data |>
  inner_join(ukb_multicat_data, by='eid') |>
  arrange(eid)

#---------------------------
# Variables related to 6147.
#---------------------------

p6147_control <- ukb_multicat_data %>%
  filter(p2207_i0 == "No") %>%
  pull(eid)

# Initialize the result data frame.
result_p6147 <- data.frame(eid = ukb_multicat_data$eid)

for (i in 1:4) {
  uniq <- paste0('6147_', i)
  cat(paste0(i, "/4 Processing ", uniq, " begin at ", Sys.time(), "\n"))
  
  code_tmp <- multi_cat_code |>
    filter(unique_id == uniq)
  
  mapping <- setNames(code_tmp$coding_num, code_tmp$coding_label)
  
  new_col_name <- paste0('tidy_p6147_', i)
  
  df_tmp <- ukb_multicat_data[, .(eid, p6147_i0)]
  df_tmp[, is_control := eid %in% p6147_control]
  
  df_with_data <- df_tmp[!is.na(p6147_i0) & p6147_i0 != ""]
  
  if (nrow(df_with_data) > 0) {
    df_split <- df_with_data[, .(
      coded_value = unlist(strsplit(p6147_i0, "\\|"))
    ), by = .(eid, is_control)]
    
    df_split[, coded_value := trimws(coded_value)]
    df_split[, coded_value := ifelse(coded_value %in% names(mapping), 
                                      mapping[coded_value], 
                                      coded_value)]
    
    df_agg <- df_split[, .(
      has_one = any(coded_value == "1", na.rm = TRUE),
      all_na = all(coded_value %in% c("NA", "EMPTY") | is.na(coded_value))
    ), by = .(eid, is_control)]
    
    # Classification logic.
    df_agg[, status := fcase(
      has_one, "1",
      is_control & all_na, "0",
      all_na, NA_character_,
      default = NA_character_
    )]
  } else {
    df_agg <- data.table(eid = integer(), status = character())
  }
  
  df_no_data <- df_tmp[is.na(p6147_i0) | p6147_i0 == ""]
  df_no_data[, status := ifelse(is_control, "0", NA_character_)]
  
  conflict_count <- 0
  if (nrow(df_agg) > 0) {
    conflict_count <- df_agg[has_one == TRUE & is_control == TRUE, .N]
  }
  
  if (conflict_count > 0) {
    warning(paste0("Warning: ", uniq, "--", conflict_count, "has 1 with CONTRL\n"))
  }
  
  df_final <- rbindlist(list(
    df_agg[, .(eid, status)],
    df_no_data[, .(eid, status)]
  ))
  
  result_p6147 <- result_p6147 %>%
    left_join(as.data.frame(df_final), by = "eid") %>%
    rename(!!new_col_name := status)
  
  result_p6147[[new_col_name]] <- factor(result_p6147[[new_col_name]], levels = c("0", "1"))
  
  cat(paste0(i, "/4 ", uniq, " completed at ", Sys.time(), "\n"))
}

summary(result_p6147) # This is consistent with the previous result.

# Save the result.
saveRDS(result_p6147, here("tmp_files/ukb_data_extract/result_p6147.rds"))

rm(p6147_control,result_p6147,i,uniq,mapping,new_col_name,df_tmp,df_with_data,
  df_split,df_agg,df_no_data,conflict_count,df_final)

#---------------------------
# Other data.
# The remaining variables are simpler to integrate.
# Directly assign cases when any value is 1, controls when all values are 0, and NA otherwise.
#---------------------------

loop_vec <- multi_cat_info |>
    filter(field_id!=6147) |>
    pull(unique_id) |>
    unique()

result_df <- data.frame(eid = ukb_multicat_data$eid)

for (i in seq_along(loop_vec)) {

  uniq <- loop_vec[i]
  cat(i, "/", length(loop_vec), " ", uniq,
      " begin at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      

  ## ---- 1. Get mapping information ----
  code_tmp <- multi_cat_code |>
    filter(unique_id == uniq)

  mapping <- setNames(code_tmp$coding_num, code_tmp$coding_label)

  field_id <- code_tmp$field_id[1]
  new_col_name <- paste0('tidy_p',uniq)

  pattern <- paste0("^p", field_id, "(?!\\d)")
  matched_cols <- grep(pattern, names(ukb_multicat_data),
                        perl = TRUE, value = TRUE)

  ## ---- 2. Extract eid and all related columns ----
  dt_tmp <- as.data.table(
    ukb_multicat_data[, c("eid", matched_cols), with = FALSE]
  )

  ## ---- 3. Convert wide to long and split | ----
  dt_long <- melt(
    dt_tmp,
    id.vars = "eid",
    measure.vars = matched_cols,
    value.name = "raw_value",
    na.rm = TRUE
  )

  dt_long <- dt_long[
    !is.na(raw_value) & raw_value != ""
  ][
    , .(coded_value = unlist(strsplit(as.character(raw_value), "\\|"))),
    by = eid
  ]

  dt_long[, coded_value := trimws(coded_value)]

  ## ---- 4. Map labels ----
  dt_long[
    coded_value %in% names(mapping),
    coded_value := mapping[coded_value]
  ]

  ## ---- 5. Aggregate at the eid level ----
  dt_agg <- dt_long[, .(
    has_one   = any(coded_value == "1", na.rm = TRUE),
    all_empty = all(coded_value == "0"),
    all_na    = all(coded_value == "NA" | is.na(coded_value)),
    has_na    = any(coded_value == "NA"),
    has_empty = any(coded_value == "0")
  ), by = eid]

  ## ---- 6. Warn about conflicts ----
  conflict_count <- dt_agg[has_na & has_empty, .N]
  if (conflict_count > 0) {
    cat(
      paste0("**Warning: ", uniq,
             " has NA coexisting with 0 (", conflict_count, ")\n")
    )
  }

  ## ---- 7. Classification logic ----
  dt_agg[, status := fcase(
    has_one, "1",
    all_empty, "0",
    all_na, NA_character_,
    default = NA_character_
  )]

  ## ---- 8. EIDs without any data ----
  eid_with_data <- unique(dt_agg$eid)

  dt_no_data <- data.table(
    eid = setdiff(dt_tmp$eid, eid_with_data),
    status = NA_character_
  )

  ## ---- 9. Merge results ----
  dt_final <- rbindlist(list(
    dt_agg[, .(eid, status)],
    dt_no_data
  ), use.names = TRUE)

  result_df <- result_df %>%
    left_join(as.data.frame(dt_final), by = "eid") %>%
    rename(!!new_col_name := status)

  result_df[[new_col_name]] <-
    factor(result_df[[new_col_name]], levels = c("0", "1"))

  cat( "DONE!", uniq,
      " completed at ", as.character(Sys.time()), "\n")
}

summary(result_df)

saveRDS(result_df, here("tmp_files/ukb_data_extract/result_multicat.rds"))

################################################
# Merge all results above into one dataset.
# As above, check whether each expected column is present.
###############################################

file_list <- c(
  'tmp_files/ukb_data_extract/result_p6147.rds',
  'tmp_files/ukb_data_extract/result_multicat.rds'
)

file_list <- here(file_list)

data_list <- lapply(file_list, readRDS)

merged_data <- data_list %>%
  reduce(inner_join, by = "eid")

ukb_val_info_multicat <- ukb_val_info |>
  filter(value_type %in% c(22))

table(ukb_val_info$value_type,useNA = "always")

not_in_merged <- setdiff(ukb_val_info_multicat$exp_code, names(merged_data))
not_in_tidy <- setdiff(names(merged_data), ukb_val_info_multicat$exp_code)

cat("Present in exp_code but missing from merged_data:", length(not_in_merged), "\n")
print(not_in_merged)

cat("\nPresent in merged_data but missing from exp_code:", length(not_in_tidy), "\n")
print(not_in_tidy)

# Check completed.
# Data extraction for this step is complete.

merged_data <- merged_data |>
  arrange(eid)

saveRDS(merged_data, file=here("processed_data/ukb_common_multicat.RDS"))

########################################################
## This is the final common-data step.
## Therefore, check all merged common-data results again.
########################################################

file_list <- c(
  'processed_data/ukb_common_continus.RDS',
  'processed_data/ukb_common_continus_special.RDS',
  'processed_data/ukb_common_cat.RDS',
  'processed_data/ukb_common_multicat.RDS'
)

file_list <- here(file_list)

data_list <- lapply(file_list, readRDS)

merged_data <- data_list %>%
  reduce(left_join, by = "eid")

ukb_val_info <- read_xlsx(
  here("rawdata/ukb_exposome.xlsx"),
  col_types = c(
    "field_id" = "text"
  )
) |>
  filter(exp_cate != "Proteome") |> # Exclude proteomics.
  filter((exp_subcate!='Brain MRI' | is.na(exp_subcate))) |> # Exclude brain MRI.
  filter(Level2 != "Derived accelerometry" | is.na(Level2)) |># Exclude accelerometry, which has already been processed.
  filter(exp_cate != "Metabolome") |> # Exclude metabolomics.
  mutate(field_id = as.numeric(field_id))

table(ukb_val_info$value_type,useNA = "always")

not_in_merged <- setdiff(ukb_val_info$exp_code, names(merged_data))
not_in_tidy <- setdiff(names(merged_data), ukb_val_info$exp_code)

cat("Present in exp_code but missing from merged_data:", length(not_in_merged), "\n")
print(not_in_merged)

cat("\nPresent in merged_data but missing from exp_code:", length(not_in_tidy), "\n")
print(not_in_tidy)

# merged_data <- merged_data |>
#   select(-any_of(setdiff(not_in_tidy, "eid")))

merged_data <- merged_data |>
  arrange(eid)

# Check completed.
# Data extraction for this step is complete.

saveRDS(merged_data, file=here("processed_data/ukb_common_all.RDS"))

# test <- readRDS(here("processed_data/ukb_common_all.RDS"))
# fwrite(test,here("processed_data/ukb_data_extract_common119.tsv"),sep = '\t')

################################################################################
# Other exposure data
################################################################################

ukb_val_info <- read_xlsx(
  here("rawdata/ukb_exposome.xlsx"),
  col_types = c("field_id" = "text"))

##################################################
## Proteomics data.
## These data were already processed during extraction.
## Proteomics variables are analyzed separately in most downstream analyses.
## Save them into the cache for later use.
#################################################

protein_data <- fread(here("processed_data/ukb_protein.tsv"))
protein_covar <- fread(here("processed_data/ukb_protein_covar.tsv"))

saveRDS(protein_data, file=here("processed_data/ukb_protein.RDS"))
saveRDS(protein_covar, file=here("processed_data/ukb_protein_covar.RDS"))

##################################################
## Metabolomics data.
## Because proteomics has special sampling characteristics, although the time point is related to baseline,
## only about 50,000 participants are available, so extract this subset separately.
## Use this subset for separate downstream analyses.
#################################################

metab_data <- fread(here('processed_data/ukb_metab.tsv'))

old_cols <- names(metab_data)[-1]

metab_data[, paste0("tidy_", old_cols) := lapply(.SD, lognmin), .SDcols = old_cols]

metab_tidy_data <- metab_data |>
  select(eid,starts_with('tidy_'))

saveRDS(metab_tidy_data, file=here("processed_data/ukb_metab.RDS"))

##################################################
## Imaging data do not require special processing.
## Save them for convenient downstream merging.
#################################################

ukb_val_info_image <- ukb_val_info |>
  filter(Time_variableid == 'p53_i2')

image_data <- fread(here('processed_data/ukb_mri.tsv'))

saveRDS(image_data, file=here("processed_data/ukb_mri.RDS"))

##################################################
## Accelerometry data.
#################################################

acc_data <- fread(here('processed_data/ukb_acc.tsv'))
saveRDS(acc_data, file=here("processed_data/ukb_acc.RDS"))

################################################
## Merge and check.
################################################

file_list <- c(
  'processed_data/ukb_common_all.RDS',
  'processed_data/ukb_metab.RDS',
  'processed_data/ukb_mri.RDS',
  'processed_data/ukb_acc.RDS'
)

file_list <- here(file_list)

data_list <- lapply(file_list, readRDS)

merged_data <- data_list |>
  reduce(left_join, by = "eid")

not_in_merged <- setdiff(ukb_val_info$exp_code, c(names(merged_data),names(protein_data)))
not_in_tidy <- setdiff(c(names(merged_data),names(protein_data)), ukb_val_info$exp_code)

cat("Present in exp_code but missing from merged_data:", length(not_in_merged), "\n")
print(not_in_merged)

cat("\nPresent in merged_data but missing from exp_code:", length(not_in_tidy), "\n")
print(not_in_tidy)

merged_data <- merged_data |>
  arrange(eid)

saveRDS(merged_data,here("processed_data/ukb_all_common_data.RDS"))

qsave(merged_data,here("processed_data/ukb_all_common_data.qs"))
qsave(protein_data,here("processed_data/ukb_all_protein_data.qs"))
qsave(protein_covar,here("processed_data/ukb_all_protein_covar.qs"))

######################################################
## Death events.
#####################################################

death_info <- fread(here('processed_data/ukb_deathtime.tsv'))

# Keep the first available value.

death_info <- death_info |>
  select(eid,p40000_i0,p40001_i0)

saveRDS(death_info,here("processed_data/ukb_death_info.RDS"))

################################################################################
# Covariate data
################################################################################

ukb_covar_data <- fread(here('processed_data/ukb_covar.tsv'))

##############################################
## Tidy basic covariates.
##############################################

ukb_covar_data <- ukb_covar_data |>
  mutate(
    tidy_p33 = make_date(
      year = p34,
      month = month(parse_date_time(p52, orders = "b")),
      day = 15
    )
  )

class(ukb_covar_data$tidy_p33)
summary(ukb_covar_data$tidy_p33)

# Sex is binary; binary covariates are treated as regular factors with 0 as the reference level.
ukb_covar_data <- ukb_covar_data |>
  mutate(tidy_p31 = factor(p31,levels = c("Female", "Male")))

class(ukb_covar_data$tidy_p31)
summary(ukb_covar_data$tidy_p31)

# tidy_p54_i0 and tidy_p54_i2 represent baseline center and follow-up imaging center.
ukb_covar_data <- ukb_covar_data |>
  mutate(
    tidy_p54_i0 = factor(na_if(p54_i0, "")),
    tidy_p54_i2 = factor(na_if(p54_i2, ""))
  )
class(ukb_covar_data$tidy_p54_i0)
class(ukb_covar_data$tidy_p54_i2)
summary(ukb_covar_data$tidy_p54_i0)
summary(ukb_covar_data$tidy_p54_i2)

# tidy_p21000_i0 stores ethnicity data.
ukb_covar_data <- ukb_covar_data |>
  mutate(tidy_p21000_i0 = case_when(
    p21000_i0 %in% c("White", "British", "Irish", "Any other white background") ~ "White",
    TRUE ~ "NotWhite"),
    tidy_p21000_i0 = factor(tidy_p21000_i0,levels = c("NotWhite","White"))
  )
class(ukb_covar_data$tidy_p21000_i0)
summary(ukb_covar_data$tidy_p21000_i0)

# Caucasian group.
ukb_covar_data <- ukb_covar_data |>
  mutate(tidy_p22006 = case_when(
    p22006 %in% c("Caucasian") ~ 1,
    TRUE ~ 0),
    tidy_p22006 = factor(tidy_p22006,levels = c(1,0))
  )

# Other groups.
ukb_covar_data <- ukb_covar_data %>%
  mutate(
    tidy_p21001_i0 = p21001_i0,  # BMI
    tidy_p74_i0 = p74_i0,        # Fasting_time  
    tidy_p20282_i0 = factor(p20282_i0,levels = 1:25),  # NMR Batch
    tidy_p25000_i2 = p25000_i2  # TIV
  )

ukb_covar_result <- ukb_covar_data %>%
  select(eid, starts_with("tidy"))

summary(ukb_covar_result)

#########################################################
# Handle missing values.
#########################################################

# Impute BMI and fasting time using sex-specific medians.
check <- ukb_covar_result[is.na(ukb_covar_result$tidy_p21001_i0),]

ukb_covar_result <- ukb_covar_result %>%
  group_by(tidy_p31) %>%
  mutate(
    tidy_p21001_i0 = if_else(
      is.na(tidy_p21001_i0),
      median(tidy_p21001_i0, na.rm = TRUE),
      tidy_p21001_i0
    )
  ) %>%
  ungroup()

summary(ukb_covar_result)

check1 <- ukb_covar_result |>
  filter(eid %in% check$eid)

ukb_covar_result <- ukb_covar_result %>%
  mutate(
    tidy_p74_i0 = if_else(
      is.na(tidy_p74_i0),
      median(tidy_p74_i0, na.rm = TRUE),
      tidy_p74_i0
    )
  )

summary(ukb_covar_result)

names(ukb_covar_result)

ukb_covar_result <- ukb_covar_result |>
  rename(
    BirthDate = tidy_p33,
    Sex = tidy_p31,
    AccCeni0 = tidy_p54_i0,
    AccCeni2 = tidy_p54_i2,
    Ethnicity = tidy_p21000_i0,
    Ethnicity_Ca = tidy_p22006,
    BMIi0 = tidy_p21001_i0,
    FastingTime = tidy_p74_i0,
    MetabBatch = tidy_p20282_i0,
    TIV = tidy_p25000_i2
  )

class(ukb_covar_result$BirthDate)
class(ukb_covar_result$Sex)
class(ukb_covar_result$AccCeni0)
class(ukb_covar_result$AccCeni2)
class(ukb_covar_result$Ethnicity)
class(ukb_covar_result$BMIi0)
class(ukb_covar_result$FastingTime)
class(ukb_covar_result$MetabBatch)
class(ukb_covar_result$TIV)

summary(ukb_covar_result)

saveRDS(ukb_covar_result,here("processed_data/ukb_covar_common.RDS"))

##############################################
## Calculate age at each time point.
##############################################
rm(list=ls())
df_covar <- readRDS(here("processed_data/ukb_covar_common.RDS"))
df_time <- fread(here("processed_data/ukb_time.tsv"))

time_cols <- setdiff(names(df_time), "eid")

df_covar <- df_covar %>%
  mutate(BirthDate = as.Date(BirthDate))

df_time <- df_time %>%
  mutate(across(-eid, as.Date))

# Merge and calculate age in completed years.
df_age <- df_time %>%
  left_join(df_covar %>% select(eid, BirthDate), by = "eid") %>%
  mutate(
    across(
      all_of(time_cols),
      ~ {
        y_diff <- year(.) - year(BirthDate)
        not_had_birthday <- format(., "%m%d") < format(BirthDate, "%m%d")
        y_diff - as.integer(not_had_birthday)
      },
      .names = "Age_{.col}"
    )
  ) %>%
  select(eid, starts_with("Age_"))

df_age <- df_age %>%
  arrange(eid)

saveRDS(df_age,here("processed_data/Age_at_each_time.RDS"))

##############################################
## Calculate sample age.
##############################################

ukb_covar2_data <- fread(here('processed_data/ukb_covar.tsv'))
df_time <- fread(here("processed_data/ukb_time.tsv")) |>
  select(eid,p53_i0)
ukb_covar2_data <- ukb_covar2_data |>
  inner_join(df_time,by='eid')
# 1. Convert date columns to Date format if needed.
# p30162_i0 belongs to basophil measurements, but some dates precede follow-up, so exclude this time variable.
date_cols <- c("p53_i0", "p30621_i0", "p23658_i0", "p30512_i0")

for(col in date_cols) {
  if(!inherits(ukb_covar2_data[[col]], "Date")) {
    ukb_covar2_data[[col]] <- as.Date(ukb_covar2_data[[col]])
  }
}

# 2. Impute missing values using the median of each column.
for(col in date_cols) {
  # Calculate the column median excluding NA.
  median_date <- median(ukb_covar2_data[[col]], na.rm = TRUE)
  # Impute missing values.
  ukb_covar2_data[[col]][is.na(ukb_covar2_data[[col]])] <- median_date
}

# 3. Create new day-difference columns.
ukb_covar2_data <- ukb_covar2_data %>%
  mutate(
    SamAge_BloodBioch = round(as.numeric(difftime(p30621_i0, p53_i0, units = "days")) / 365.25, 1),
    SamAge_BloodMetab = round(as.numeric(difftime(p23658_i0, p53_i0, units = "days")) / 365.25, 1),
    SamAge_UrineBioch = round(as.numeric(difftime(p30512_i0, p53_i0, units = "days")) / 365.25, 1)
  )

ukb_covar2_data <- ukb_covar2_data %>%
  select(eid, starts_with("SamAge_"))

summary(ukb_covar2_data)

saveRDS(ukb_covar2_data,here("processed_data/covar_sample_age.RDS"))

# Check p20282 against metabolomics data to identify missing values.
# Revisit later to decide whether imputation is needed.

# Merge covariates together.

file_list <- c(
  'processed_data/ukb_covar_common.RDS',
  'processed_data/Age_at_each_time.RDS',
  'processed_data/covar_sample_age.RDS'
)

file_list <- here(file_list)

data_list <- lapply(file_list, readRDS)

merged_data <- data_list %>%
  reduce(inner_join, by = "eid")

summary(merged_data)

saveRDS(merged_data,here("processed_data/ukb_all_covar_common.RDS"))
qsave(merged_data,here("processed_data/ukb_all_covar_common.qs"))
