library(MatchIt)
library(dplyr)

cov_df <- read.csv("path/to/cov_data.csv")
cov_df_selected <- cov_df %>% mutate(across(c(Sex, Race), as.factor),across(c(Age, BMI), ~ as.numeric(scale(.x))))
target_file_lst <- c("disease_chapter1","disease_chapter2","disease_chapter3")
ratio <- 1
my_formula <- "Age+BMI"

for (tgt_file in target_file_lst) {
    tmp_tgt_df <- read.csv(paste0("path/to/disease_data/",tgt_file,".csv"))
    df <- merge(tmp_tgt_df,cov_df_selected,by="id")

    set.seed(2024)
    
    match.it <-  matchit(
      formula = as.formula(paste("y ~", my_formula)),
      data = df,
      method = "nearest",
      ratio = ratio,distance = "mahalanobis",link = "logit",
      exact = ~ Sex+Race
    )
    df.match <- match.data(match.it)

    result <- df.match %>%
      group_by(subclass) %>%
      mutate(
        yrs = yrs[y == 1]
      ) %>%
      ungroup()
    
    write.csv(result, paste0("path/to/matched_results/",tgt_file,".csv"), row.names = FALSE)
}
