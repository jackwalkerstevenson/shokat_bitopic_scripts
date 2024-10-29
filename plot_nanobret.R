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
  # if string ends in an underscore and then a numeral, remove it
  stringr::str_replace(text, "_\\d+$", "")
}
raw_data <- input_path |> 
  readxl::excel_sheets() |> 
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
# function to plot specified targets of a treatment-----------------------------------------------
plot_nanobret_trt_tgts <- function(data, trt, tgts, plot_title){
  trt_display_name = display_names_treatments[trt]
  trt_data <- data |>
    doseplotr::filter_trt_tgt(trt=trt, tgt = tgts)
  trt_tgts <- unique(trt_data$target) # names of targets actually present
  x_min <- floor(min(trt_data$log_dose))
  x_max <- ceiling(max(trt_data$log_dose))
  x_limits <- c(x_min, x_max)
  trt_data_summary <- trt_data |> 
    group_by(target, log_dose) |> 
    doseplotr::summarize_response()
  model_predictions <- trt_data |>
    dplyr::group_by(target, log_dose) |> 
    doseplotr::summarize_models(response_col="response", bounded=FALSE) |>
    doseplotr::get_predictions(response_col="response")
  plot_output <- trt_data_summary |>
    ggplot(aes(x = log_dose, y = mean_response, color = target, shape = target)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymax = .data$mean_response + .data$sem,
                      ymin = .data$mean_response - .data$sem,
                      width = .data$w)) +
    geom_line(data = model_predictions, aes(y = response), linewidth = .75, alpha = 0.8) +
    scale_x_continuous(breaks = scales::breaks_width(1),
                       minor_breaks = minor_breaks_log(x_limits[1],
                                                       x_limits[2])) +
    scale_shape_manual(values = shape_map_targets,
                       name = "competitor",
                       labels = display_names_targets) +
    scale_color_manual(values = color_map_targets,
                       name = "competitor",
                       labels = display_names_targets) +
    guides(
      # ggprism guides to end at last tick
      x = ggprism::guide_prism_offset_minor(),
      y = ggprism::guide_prism_offset(),
    ) +
    theme_prism() +
    theme(plot.background = element_blank(), # transparent
          legend.title = element_text()) +
    labs(x = str_glue("log10[{trt_display_name}] (M)"),
         y = "BRET ratio (mBU)",
         title = str_glue("{plot_title}"))
  return(plot_output)
}
# function to plot each target of a treatment along with control----------------
plot_nanobret_trt <- function(data, trt){
  trt_display_name = display_names_treatments[trt]
  trt_data <- data |>
    doseplotr::filter_trt_tgt(trt=trt)
  trt_tgts <- unique(trt_data$target) # names of all targets
  noncontrol_targets <- trt_tgts[trt_tgts != control_target]
  print(noncontrol_targets)
  for(tgt in noncontrol_targets){
    plot_nanobret_trt_tgts(data, trt=trt,
                           tgts = c(control_target, tgt),
                           plot_title = display_names_targets[tgt]) |> 
      doseplotr::save_plot(
        filename=str_glue("{output_directory}/plot_nanobret_trt_{janitor::make_clean_names(trt)}_tgt_{janitor::make_clean_names(tgt)}_{get_timestamp()}.{plot_type}")
        )
  }
}
# loop to plot all targets of each treatment along with control-----------------
treatments <- unique(raw_data$treatment)
for(trt in treatments){
  plot_nanobret_trt(raw_data, trt) 
}