library(survival)
library(data.table)
library(broom)

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

    all_model_results <- list()

    for (bin in bins_level) {
        data <- matched_data[matched_data$bin_new==bin,]
        for (expo in expo_f_lst) {
            temp_data <- merge(
                data[,c("id","y","subclass")],
                expo_df[,c("id",..expo)],
                by="id"
            )
            setnames(temp_data, old = expo, new = "expo_value")

            model <- clogit(y ~ expo_value + strata(subclass), data = temp_data)
            res <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
            res <- res[, c("estimate","std.error","statistic", "conf.low", "conf.high", "p.value")]
            res$bin <- bin
            res$exposure <- expo
            all_model_results[[paste(bin, expo, sep = "_")]] <- res
        }
    }

    final_res_table <- do.call(rbind, all_model_results)
    write.csv(final_res_table,paste0("path/to/or_results/",disease,".csv"),row.names=FALSE)
}
