library(data.table)
library(dplyr)

expo_df <- fread("path/to/exposure_data.csv")
target_file_lst <- c("disease1", "disease2", "disease3")
expo_f_lst <- c("expo1","expo2","expo3")

for (disease in target_file_lst) {
  matched_data <- fread(paste0("path/to/matched_results/",disease,".csv"))

  bins <- seq(-16, 0, by = 1)
  lower_bounds <- abs(bins[2:length(bins)])
  upper_bounds <- abs(bins[1:(length(bins)-1)])
  bins_level <- paste0(upper_bounds, "_", lower_bounds)

  matched_data <- matched_data[matched_data$yrs <= 16 & matched_data$yrs >= 0, ]
  matched_data$yrs <- -matched_data$yrs
  matched_data$bin_old <- cut(matched_data$yrs, bins)
  matched_data$bin_new <- gsub("\\(|\\)|\\[|\\]", "", matched_data$bin_old) |>
    strsplit(",") |>
    sapply(function(x) paste0(abs(as.numeric(x[1])), "_", abs(as.numeric(x[2]))))

  disease_result <- data.frame(expo = expo_f_lst)
  for (bin in bins_level) {disease_result[[bin]] <- NA}
  for (bin in bins_level) {
    data <- matched_data[matched_data$bin_new==bin,]
    case_ids <- data[data$y == 1,]$id
    control_ids <- data[data$y == 0,]$id
    
    case_data <- expo_df[expo_df$id %in% case_ids,][,-1]
    case_means <- colMeans(case_data, na.rm = TRUE)

    control_data <- expo_df[expo_df$id %in% control_ids,][,-1]
    control_means <- colMeans(control_data, na.rm = TRUE)
    control_sd <- apply(control_data, 2, sd, na.rm = TRUE)

    case_control_zscore <- (case_means-control_means)/control_sd
    disease_result[[bin]] <- case_control_zscore
  }

  write.csv(disease_result,paste0("path/to/z_score_results/",disease,".csv"),row.names=FALSE)
}
