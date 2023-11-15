# plot_CTG spinoff 2023-09-12 for plotting multiple targets by linetype
# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(readxl) # for excel file handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis transforms
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(doseplotr) # you bet
# import global parameters and clear environment--------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_plot_CTG.R")
# import and preprocess data----------------------------------------------------
# create input and output directories, since git doesn't track empty directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
if(exists("input_filename")){ # if using single file import, read and preprocess
  plate_data <- readr::read_csv(input_filename) |>
    preprocess_plate_data()
} else {plate_data <- import_plates(input_directory)} # else import platemaps
# filter and validate imported data---------------------------------------------
if(exists("treatments")){
  plate_data <- filter_validate_reorder(plate_data, "treatment", treatments)
}
if(exists("targets")){
  plate_data <- filter_validate_reorder(plate_data, "target", targets)
}
plot_data <- plate_data |> filter(log_dose != -Inf) # data without untreated
# generate data-dependent global plot parameters--------------------------------
if (!exists("treatments")){ # if treatments not specified, use all treatments
  treatments <- as.vector(unique(plot_data$treatment))}
if (!exists("targets")){ # if targets not specified, use all targets
  targets <- as.vector(unique(plot_data$target))}
# find x-axis min/max values for consistent zoom window between all plots
if(override_x_lim){
  x_limits <- manual_x_lim
} else{
  x_min <- floor(min(plot_data$log_dose))
  x_max <- ceiling(max(plot_data$log_dose))
  x_limits <- c(x_min, x_max)
}
# x_limits <- c(-12,-4) # manual x limit backup
# fit models and report model parameters----------------------------------------
model_summary <- summarize_models(plot_data,
                                  response_col = "response_norm",
                                  rigid = rigid) |> # use global rigid parameter
  dplyr::select(-model) |>  # remove actual model from report
  mutate(across(where(is.numeric), \(x){signif(x, digits = 4)}))
write_csv(model_summary,
          str_glue("output/CTG_model_summary_{get_timestamp()}.csv"))
# helper to plot multiple targets at once---------------------------------------
plot_multitarget <- function(df, tgts,
                        response_col = "response",
                        color_map = NULL,
                        rigid = FALSE,
                        x_limits = NULL,
                        legend_title = ggplot2::waiver(),
                        legend_labels = ggplot2::waiver(),
                        ...){
  # get data for specified targets
  data <- df |> filter_trt_tgt(tgt = tgts)
  # calculate x limits if not specified
  if(is.null(x_limits)){
    x_min <- floor(min(data$log_dose))
    x_max <- ceiling(max(data$log_dose))
    x_limits <- c(x_min, x_max)}
  # summarize actual data to plot points and error bars
  data_summary <- data |>
    dplyr::group_by(.data$target, .data$treatment, .data$log_dose) |> # group by treatment
    summarize_response(response_col = response_col)
  # fit models to data to plot model prediction curves
  model_predictions <- data |>
    dplyr::group_by(.data$target, .data$treatment, .data$log_dose) |>
    summarize_models(response_col = response_col, rigid = rigid) |>
    get_predictions(response_col = "mean_response") # name for plotting
  # set up color parameters based on number of treatments
  num_treatments <- length(unique(data_summary$treatment))
  vr <- viridis_range(num_treatments)
  vr_begin <- vr[[1]]
  vr_end <- vr[[2]]
  vr_option <- vr[[3]]
  # plot points and error bars from the summarized data
  p <- {ggplot2::ggplot(data_summary,
                        ggplot2::aes(x = .data$log_dose,
                                     y = .data$mean_response,
                                     shape = .data$treatment,
                                     color = .data$treatment)) +
      ggplot2::geom_point(size = 3) +
      {if(is.null(color_map)){ # optionally manually specify colors
        viridis::scale_color_viridis(discrete = TRUE, option = vr_option,
                                     begin = vr_begin, end = vr_end,
                                     name = legend_title,
                                     labels = legend_labels)
      } else{
        ggplot2::scale_color_manual(values = color_map,
                                    name = legend_title,
                                    labels = legend_labels)
      }}} |>
    base_dose_response(x_limits = x_limits,
                       legend_title = legend_title,
                       legend_labels = legend_labels,
                       ...)
  # plot model predictions for models that were fit successfully
  p <- p +
    ggplot2::geom_line(data = model_predictions, aes(linetype = target), linewidth = .75, alpha = 0.8)
  return(p)
}
# plot data for multiple targets------------------------------------------
treatment_legend_title <- if(override_treatment_title){
  treatment_title} else "treatment"
# get all treatments for specified targets to set legend length
tgt_treatments <- as.vector(unique((plot_data |>
                                      filter_trt_tgt(tgt = targets))$treatment))
if(manually_relabel_treatments){
  tgt_treatments <- Vectorize(get_display_name, vectorize.args = "name")(
    tgt_treatments, display_names_treatments, TRUE)}
tgt <- targets
plot_multitarget(plot_data, tgt,
              rigid = rigid, # global param: rigid low-dose asymptote
              grid = grid, # global param: background grid on plot
              no_legend = no_legend, # global param: whether to omit legend
              x_limits = get_if(x_limits, global_x_lim),
              response_col = "response_norm", # CTG uses response_norm
              ylab = "luminescence (% of untreated)", # CTG = luminescence assay
              legend_title = treatment_legend_title,
              # if relabeling treatments, get display names for legend
              legend_labels = get_if(display_names_treatments,
                                     manually_relabel_treatments,
                                     otherwise = ggplot2::waiver()),
              #if relabeling, get display name for title
              plot_title = "Cell viability",
              # if manually setting colors of treatments, get color map
              color_map = get_if(color_map_treatments,
                                 manually_recolor_treatments),
              # if manually setting shapes of treatments, get shape map
              shape_map = get_if(shape_map_treatments,
                                 manually_reshape_treatments)
  )|>
    save_plot(
      str_glue("output/CTG_multitarget_{get_timestamp()}.{plot_type}"),
      no_legend = no_legend,
      width = 8,
      legend_len = longest(c(tgt_treatments, treatment_legend_title)))
