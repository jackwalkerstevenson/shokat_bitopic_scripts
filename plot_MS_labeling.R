# plot mass spec labeling data
# Jack Stevenson started 2025-05-27
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MS_labeling.R"
source(params_path)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MS_labeling.R", output_dir)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
raw_data <- readxl::read_excel(fs::path(input_dir, input_filename))
all_data <- raw_data |> 
  tidyr::pivot_longer(
    # find deduplication suffixes introduced by read_excel, e.g. "...1"
    cols = matches(".*\\.\\.\\.[0-9]+"),
    names_to = "treatment", # column names will become the value of "treatment"
    names_pattern = "(.*)\\.\\.\\.[0-9]+", # capture name before dedupe suffix
    values_to = "percent_labeling"
  ) |> 
  janitor::clean_names()
# report processed data
write_csv(all_data,
          fs::path(output_dir,
                   str_glue("MS_labeling_all_data_{get_timestamp()}.csv")))
# fit models-----------------------------------------------
#chatgpt
# Nonlinear model function
fit_model <- function(df) {
  nls(
    percent_labeling ~ plateau * (1 - exp(-k * elapsed_time_from_first_injection_min)),
    data = df,
    start = list(plateau = 100, k = 0.001)
  )
}

# # test data for one model
# data_62 <- all_data |> 
#   filter(treatment == "KL-ELN5-62")
# test <- fit_model(data_62)

# test data for one model
data_59 <- all_data |> 
  filter(treatment == "KL-ELN5-59")
test <- fit_model(data_59)

# Fit models for each treatment-----------------------------------------------
# Function to safely fit exponential model (Claude)
fit_exponential_safe <- function(data) {
  max_time <- max(data$elapsed_time_from_first_injection_min)
  k_estimate <- 0.01  # approx for majority labeling in 100 min
  
  tryCatch({
    model <- nls(percent_labeling ~ 100 * (1 - exp(-k * elapsed_time_from_first_injection_min)), 
                 data = data,
                 start = list(k = k_estimate))
    
    # Generate predictions
    time_seq <- seq(0, max_time, length.out = 100)
    predictions <- predict(model, newdata = data.frame(elapsed_time_from_first_injection_min = time_seq))
    
    return(data.frame(
      elapsed_time_from_first_injection_min = time_seq, 
      predicted = predictions,
      treatment = unique(data$treatment)
    ))
    
  }, error = function(e) {
    # Return empty data frame if fit fails
    return(data.frame(
      elapsed_time_from_first_injection_min = numeric(0), 
      predicted = numeric(0),
      treatment = character(0)
    ))
  })
}
# plot labeling by compound-----------------------------------------------
all_data |> 
  ggplot(aes(x = elapsed_time_from_first_injection_min,
             y = percent_labeling,
             color = treatment)) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar") +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 2
  ) +
  labs(x = "time (min)",
       y = "percent of protein labeled") +
  theme_prism() +
  theme(plot.background = element_blank()) # transparent background

ggsave(str_glue(
  "{output_dir}/MS_labeling_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 7, height = 5)
