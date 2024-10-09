# plot nanobret dose-response data
# Jack Stevenson started 2024-10-03
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
library(openxlsx) # for excel import and merged cell handling
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_nanobret.R"
source(params_path)
input_path <- str_glue("{input_directory}{input_filename}")
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_nanobret.R", output_directory)
# write timestamped input file to output
doseplotr::file_copy_to_dir(input_path, output_directory)
# set up input and output directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
# function to remove duplicate suffixes from values that used to be colnames
remove_dup_suffix <- function(text){
  # if string contains an underscore and then a numeral, remove it
  stringr::str_replace(text, "_\\d+$", "")
}
raw_data <- input_path |> 
  excel_sheets() |> 
  set_names() |> 
  purrr::map(function(sheet_name){
    # split merged colname cells into duplicates, not treating them as colnames
    temp_df <- read.xlsx(xlsxFile = input_path,
                         sheet = sheet_name,
                         colNames = FALSE,
                         fillMergedCells = TRUE)
    # treat them as colnames again. now in unhappy duplicate colname state
    colnames(temp_df) <-  make.names(temp_df[1,], unique = FALSE)
    temp_df <- temp_df[-1, ] |> # remove first row, which was used as colnames
      janitor::clean_names() |> # make names nice and add duplicate suffixes
      tidyr::pivot_longer( # pivot all target observations
        cols = 2:last_col(), # all columns after the first one
        names_to = "target",
        values_to = "response"
      ) |> 
      # make treatment name into a column
      tidyr::pivot_longer(
        cols = 1,
        names_to = "treatment",
        values_to = "dose_nM"
      )
    return(temp_df)
  }) |> 
  bind_rows() |>  # combine data from all sheets into one dataframe
  # remove leftover duplicate suffixes from former colnames
  dplyr::mutate(
    target = stringr::str_replace(target, "_\\d+$", ""),
    dose_nM = as.numeric(dose_nM),
    response = as.numeric(response),
    # todo: replace this with a filter_validate_reorder
    treatment = fct_inorder(treatment),
    target = fct_inorder(target)
  ) |> 
  doseplotr::make_log_dose()

# report raw data
write_csv(raw_data,
          fs::path(output_directory,
                   str_glue("nanobret_raw_data_{get_timestamp()}.csv")))
# plot just one target as a test-------------------------------------------
data <- raw_data |>
  filter_trt_tgt(trt = "x_asciminib_bodipy_n_m", tgt = "dmso")
x_min <- floor(min(data$log_dose))
x_max <- ceiling(max(data$log_dose))
x_limits <- c(x_min, x_max)
data_summary <- data |>
  dplyr::group_by(.data$treatment, .data$log_dose) |> # group by treatment
  doseplotr::summarize_response()
# this doesn't work rn
model_predictions <- data |>
  dplyr::group_by(.data$treatment, .data$log_dose) |>
  summarize_models(response_col = response_col, rigid = TRUE) |>
  get_predictions(response_col = "mean_response") # name for plotting
test_plot_data |>
  ggplot(aes(x = log_dose, y = response)) +
  geom_point() +
  theme_prism() +
  theme(legend.title = element_text()) +
  labs(x = "[asciminib-BODIPY tracer] (logmolar)",
       y = "BRET ratio (mBU)",
       title = "DMSO")
