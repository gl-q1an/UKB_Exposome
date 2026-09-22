# ==============================================================================
# Pipeline: Dose-Response Analysis using Restricted Cubic Splines (RCS)
==============================================================================

library(rms)
library(data.table)

# 1. Configuration and Initialization ------------------------------------------
valid_exposure_vars <- setdiff(names(exposure), "eid")

my_covariates <- c(
  'Sex',
  'AccCeni0',
  'Age_p53_i0'
)

time_var  <- "diag_time_p53_i0"
event_var <- "status"

df2 <- data_merged

# 2. Data Preprocessing --------------------------------------------------------
# Truncate continuous exposure variables at 1st and 99th percentiles to handle outliers
for(v in valid_exposure_vars){
  q1  <- quantile(df2[[v]], 0.01, na.rm = TRUE)
  q99 <- quantile(df2[[v]], 0.99, na.rm = TRUE)
  df2[[v]] <- pmin(pmax(df2[[v]], q1), q99)
}

# Ensure categorical covariates are properly encoded as factors
df2$AccCeni0  <- factor(df2$AccCeni0)
df2$Sex       <- factor(df2$Sex)
df2$Ethnicity <- factor(df2$Ethnicity)

# Initialize data distribution settings for the rms package
dd <- datadist(df2)
options(datadist = "dd")

# 3. Iterative RCS Modeling ----------------------------------------------------
results_list <- list()
error_list   <- list()

for(exposure_var in valid_exposure_vars){
  
  cat("\n====================================\n")
  cat("Processing variable:", exposure_var, "\n")
  cat("====================================\n")
  
  tryCatch({
    
    # Calculate quantiles for knots (10th, 50th, 90th percentiles)
    knots <- quantile(
      df2[[exposure_var]],
      probs = c(0.10, 0.50, 0.90),
      na.rm = TRUE
    )
    
    # Construct model formula dynamically based on unique knots count
    if(length(unique(knots)) < 3){
      cat("Warning: Duplicated knots detected. Defaulting to df = 3\n")
      fml <- as.formula(sprintf(
        "Surv(%s, %s) ~ rcs(%s, df = 3) + %s",
        time_var, event_var, exposure_var,
        paste(my_covariates, collapse = " + ")
      ))
    } else {
      fml <- as.formula(sprintf(
        "Surv(%s, %s) ~ rcs(%s, c(%f,%f,%f)) + %s",
        time_var, event_var, exposure_var,
        knots[1], knots[2], knots[3],
        paste(my_covariates, collapse = " + ")
      ))
    }
    
    # Fit the Cox proportional hazards model
    fit_rcs <- cph(
      fml,
      data = df2,
      x = TRUE, y = TRUE,
      method = "breslow"
    )
    
    # Set the 5th percentile as the reference value (HR = 1)
    ref_value <- quantile(df2[[exposure_var]], probs = 0.05, na.rm = TRUE)
    dd$limits[[exposure_var]][2] <- ref_value
    fit_rcs <- update(fit_rcs) # Update model matrix with new reference limit
    
    # Extract ANOVA p-values (Overall and Non-linear effects)
    an        <- anova(fit_rcs)
    p_overall <- an[1, "P"]
    
    nonlin_idx <- grep("Nonlinear", rownames(an))
    p_nonlin   <- if(length(nonlin_idx) > 0) an[nonlin_idx, "P"] else NA
    
    cat("Overall P   =", round(p_overall, 4), "\n")
    cat("Non-linear P =", round(p_nonlin, 4), "\n")
    
    # Predict HR estimates across the continuous spectrum of exposure
    HR_plot <- Predict(
      fit_rcs,
      name = exposure_var,
      fun = exp,
      ref.zero = TRUE
    )
    
    pred_df <- as.data.frame(HR_plot)
    
    # Append metadata to the predicted values dataframe
    pred_df$exposure_variable <- exposure_var
    pred_df$p_overall         <- p_overall
    pred_df$p_nonlinear        <- p_nonlin
    pred_df$ref_value         <- ref_value
    
    if(length(knots) == 3){
      pred_df$knot_10th <- knots[1]
      pred_df$knot_50th <- knots[2]
      pred_df$knot_90th <- knots[3]
    }
    
    # Standardize output column headers
    colnames(pred_df)[1] <- "exposure_value"
    colnames(pred_df)[6] <- "HR"
    colnames(pred_df)[7] <- "lower_CI"
    colnames(pred_df)[8] <- "upper_CI"
    
    results_list[[exposure_var]] <- pred_df
    
  }, error = function(e){
    # Error log capture
    cat("❌ ERROR in:", exposure_var, "\n")
    cat("Message:", conditionMessage(e), "\n")
    
    error_list[[exposure_var]] <<- list(
      variable = exposure_var,
      error = conditionMessage(e)
    )
  })
}

# 4. Merge and Save Results ----------------------------------------------------
final_results <- rbindlist(results_list, fill = TRUE)
error_results <- rbindlist(error_list, fill = TRUE)

# Show data preview
print(head(final_results, 20))
print(dim(final_results))

# Export outputs using relative pathing for repository portability
output_dir <- "output"
if(!dir.exists(output_dir)) dir.create(output_dir)

fwrite(final_results, file.path(output_dir, "dose_response_results.csv"))
if(nrow(error_results) > 0){
  fwrite(error_results, file.path(output_dir, "dose_response_errors.log"))
}
