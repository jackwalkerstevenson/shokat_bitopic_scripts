# plot data from primary CFU experiments
# Jack Stevenson started 2025-09-08
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(MASS) # for negative binomial model glm.nb
library(multcomp) # for multiple comparison test glht
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_CFU.R"
source(params_path)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_CFU.R", output_dir)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# write timestamped input file to output
input_path <- str_glue("{input_dir}{input_filename}")
doseplotr::file_copy_to_dir(input_path, output_dir)
# import, preprocess and report data-----------------------------------------------
all_data <- readxl::read_excel(input_path) |>
  janitor::clean_names() |> 
  tidyr::pivot_longer(
    cols = matches("(cfu_g|bfu_e)_\\d+"), # get columns starting with bfu or cfu
    names_to = "measurement", # store column names in "measurement"
    names_pattern = "(.+)_\\d+", # just get "bfu_e" or "cfu_g" without number
    values_to = "count"
  ) |> 
  dplyr::mutate(measurement = as.factor(measurement)) |> 
  doseplotr::filter_validate_reorder("treatment", treatments) |> 
  dplyr::mutate(
    count_normalized_percent = count / mean(count[treatment == "DMSO"]) * 100,
    .by = measurement)
# report processed data
write_csv(all_data,
          fs::path(output_dir,
                   str_glue("CFU_data_{get_timestamp()}.csv")))
# exploratory stats, deprecated-----------------------------------------------
# data_10uM <- all_data |> 
#   dplyr::filter(dose_nm == 10000 | treatment == "DMSO")
# data_10uM_bfu <- data_10uM |> 
#   dplyr::filter(measurement == "bfu_e"|treatment == "DMSO")
# data_10uM_cfu <- data_10uM |> 
#   dplyr::filter(measurement == "cfu_g"|treatment == "DMSO")
# # plot everything in 10 µM to start
# jitter_plot_10uM <- data_10uM |> 
#   ggplot(aes(x = count, y = measurement, color = treatment)) +
#   # geom_histogram()
#   geom_jitter() +
#   scale_color_manual(values = color_map_treatments) +
#   theme_prism()
# # calculate basic summary stats
# summary_stats_10uM <- data_10uM |> 
#   group_by(measurement, treatment) |> 
#   summarize(
#     n = n(),
#     mean = mean(count),
#     sd = sd(count),
#     variance = var(count),
#     cv_percent = (sd/mean) * 100,
#     variance_to_mean = variance / mean,
#     .groups = "drop"
#    )
# goal
# forget about poisson model, already tried it, seems overdispersed
# for each combination of measurement and dose_nm (treating as whole separate datasets):
# fit negative binomial model
# e.g. MASS::glm.nb(count ~ treatment, data = appropriate_data_subset)
# run Dunnett test, get log FC estimates and P values for each treatment vs DMSO
# e.g. multcomp::glht(appropriate_data_subset, linfct = mcp(treatment = "Dunnett"), alternative = "two.sided")
# make a dataframe with columns "measurement", dose_nm" "treatment", "dunnett_log_fc_estimate", "dunnett_adjusted_pval"
# bind to other dataframes from other measurement/dose combos

# poisson model testing
# poisson_model_10uM_bfu <- glm(count ~ treatment, 
#                      data = data_10uM_bfu, 
#                      family = poisson(link = "log"))
# residual_deviance <- poisson_model_10uM_bfu$deviance
# df_residual <- poisson_model$df.residual
# overdispersion_param <- residual_deviance / df_residual

# negative binomial model testing
# nb_model_10uM_bfu <- MASS::glm.nb(count ~ treatment, data = data_10uM_bfu)
# # dunnett_test_10uM_bfu <- multcomp::glht(poisson_model_10uM_bfu, 
# dunnett_test_10uM_bfu <- multcomp::glht(nb_model_10uM_bfu, 
#                               linfct = mcp(treatment = "Dunnett"),
#                               alternative = "two.sided")
# dunnett_ci_10uM_bfu <- confint(dunnett_test_10uM_bfu, level = 0.95)
# plot(poisson_model_10uM_bfu$fitted.values, residuals(poisson_model, type = "deviance"),
#      xlab = "Fitted Values", ylab = "Deviance Residuals",
#      main = "Residuals vs Fitted")
# abline(h = 0, lty = 2)

# function to fit negative binomial, run dunnett test and return info by treatment-----------------------------------------------
analyze_data_subset <- function(data, measurement_name, dose_level){
  tryCatch({
    test_data <- data |> 
      dplyr::filter(measurement == measurement_name | treatment == "DMSO",
                    dose_nm == dose_level | treatment == "DMSO")
    # fit negative binomial model
    nb_model <- MASS::glm.nb(count ~ treatment, data = test_data)
    # run Dunnett multiple comparison test
    dunnett_test <- multcomp::glht(nb_model, 
                                   linfct = mcp(treatment = "Dunnett"),
                                   alternative = "two.sided")
    dunnett_summary <- summary(dunnett_test) # get test info
    coef_estimates <- coef(dunnett_summary) # coefs are named by treatment " - DMSO"
    treatment_names <- names(coef_estimates) |> # extract treatment names
      str_extract(".*(?= - DMSO)"[1])
    coef_estimates <- unname(coef_estimates) # names no longer needed
    p_values <- dunnett_summary$test$pvalues # extract p-values
    results_df <- data.frame(
      measurement = as.factor(measurement_name),
      dose_nm = dose_level,
      treatment = as.factor(treatment_names),
      dunnett_log_fc_estimate = coef_estimates,
      dunnett_adjusted_pval = p_values
    )
    annotated_results_df <- results_df |> 
      dplyr::mutate(
        significance = dplyr::case_when(
          dunnett_adjusted_pval < 0.001 ~ "***",
          dunnett_adjusted_pval < 0.01 ~ "**",
          dunnett_adjusted_pval < 0.05 ~ "*",
          .default = NA
        )
      )
    return(annotated_results_df)
    }, error = function(e){
      warning(str_glue("Analysis failed for measurement {measurement_name}, dose {dose_level}: {e$message}"))
      return(data.frame(
        measurement = character(0),
        dose_nm = numeric(0),
        treatment = character(0),
        dunnett_log_fc_estimate = numeric(0),
        dunnett_adjusted_pval = numeric(0)
        ))
      }
  )
}
# function to fit and test models separately on each measurement/dose subset-----------------------------------------------
fit_nb_dunnett_all_data <- function(data){
  # get all combinations of measurement and dose
  measurement_dose_combinations <- all_data |>
    dplyr::select(measurement, dose_nm) |>
    dplyr::distinct() |> 
    dplyr::filter(!is.na(dose_nm)) # remove the DMSO-only subset
  
  all_results <- list()
  for(i in 1:nrow(measurement_dose_combinations)){
      current_measurement <- measurement_dose_combinations$measurement[i]
      current_dose <- measurement_dose_combinations$dose_nm[i]
      
      cat(paste("Processing:", current_measurement, "at dose", current_dose, "nm\n"))
      
      # analyze this subset
      subset_results <- analyze_data_subset(data, current_measurement, current_dose)
      
      # Store results
      if(nrow(subset_results) > 0) {
        all_results[[i]] <- subset_results
      }
    }
    
    # Combine all results
    if(length(all_results) > 0) {
      final_results <- bind_rows(all_results)
      return(final_results)
    } else {
      warning("No successful analyses completed")
      return(data.frame(
        measurement = character(0),
        dose_nm = numeric(0),
        treatment = character(0),
        dunnett_log_fc_estimate = numeric(0),
        dunnett_adjusted_pval = numeric(0)
      ))
    }
  }
# fit models and report stats-----------------------------------------------
nb_dunnett_stats <- fit_nb_dunnett_all_data(all_data)
write.csv(nb_dunnett_stats, fs::path(output_dir, str_glue("CFU_negative_binomial_dunnett_stats_{doseplotr::get_timestamp()}.csv")))
# set up plot constants-----------------------------------------------
bar_plot_width = 7
bar_plot_height = 5
summary_data_with_stats <- all_data |> 
  dplyr::group_by(measurement, dose_nm, treatment) |> 
  dplyr::summarize(
    n = n(),
    mean_count_normalized_percent = mean(count_normalized_percent),
    sem_count_normalized_percent = sd(count_normalized_percent) / sqrt(n)
  ) |> 
  dplyr::left_join(nb_dunnett_stats,
            by = join_by(measurement, dose_nm, treatment))
# plot 1 µM only-----------------------------------------------
summary_data_with_stats |> 
  dplyr::filter(treatment == "DMSO" | dose_nm == 1000) |> 
  ggplot(aes(x = treatment, y = mean_count_normalized_percent, fill = measurement, label = significance)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_count_normalized_percent - sem_count_normalized_percent,
                    ymax = mean_count_normalized_percent + sem_count_normalized_percent),
                width = 0.5,
                position = position_dodge(width = .9)) +
  # stat_summary(fun = "mean",
  #              geom = "col",
  #              position = "dodge") +
  # stat_summary(fun.data = "mean_se",
  #              geom = "errorbar",
  #              position = position_dodge(width = 1),
  #              width = 0.5) + 
  # geom_text(position = position_dodge(width = 1),
  #           vjust = 0.5,
  #           na.rm = TRUE) +
  scale_fill_manual(values = color_map_cfu,
                    labels = display_names_cfu) +
  theme_prism() +
  labs(x = "treatment (1 µM)",
       y = "colony-forming units (% of control)",
       title = "Inhibition of primary PBMCs",
       fill = "colony type") +
  theme(plot.background = element_blank(),
        legend.title = element_text(),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust = 0))
ggsave(str_glue(
  "{output_dir}/CFU_1uM_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = bar_plot_width, height = bar_plot_height)
# plot 10 µM only-----------------------------------------------
summary_data_with_stats |> 
  dplyr::filter(treatment == "DMSO" | dose_nm == 10000) |> 
  ggplot(aes(x = treatment, y = mean_count_normalized_percent, fill = measurement, label = significance)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_count_normalized_percent - sem_count_normalized_percent,
                    ymax = mean_count_normalized_percent + sem_count_normalized_percent),
                width = 0.5,
                position = position_dodge(width = .9)) +
  # stat_summary(fun = "mean",
  #              geom = "col",
  #              position = "dodge") +
  # stat_summary(fun.data = "mean_se",
  #              geom = "errorbar",
  #              position = position_dodge(width = 1),
  #              width = 0.5) + 
  # geom_text(position = position_dodge(width = 1),
  #           vjust = 0.5,
  #           na.rm = TRUE) +
  scale_fill_manual(values = color_map_cfu,
                    labels = display_names_cfu) +
  theme_prism() +
  labs(x = "treatment (10 µM)",
       y = "colony-forming units (% of control)",
       title = "Inhibition of primary PBMCs",
       fill = "colony type") +
  theme(plot.background = element_blank(),
        legend.title = element_text(),
        axis.text.x = element_text(angle = -45, vjust = 0.6, hjust = 0))
ggsave(str_glue(
  "{output_dir}/CFU_10uM_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = bar_plot_width, height = bar_plot_height)