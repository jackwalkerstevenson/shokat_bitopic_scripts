#' ---
#'title: "plot_selectscreen"
#'author: "Jack Stevenson"
#' ---
#'script for analyzing Thermo SelectScreen data
#'copied updates from plot_CTG.R 2023-08-09
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# import global parameters---------------------------------------------------------
params_path <- "parameters/parameters_plot_selectscreen.R"
source(params_path)
# import and tidy data---------------------------------------------------------
# create input and output directories, since git doesn't track empty directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
plot_data <- import_selectscreen(input_filename)
doseplotr::file_copy_to_dir(input_filename, output_directory) # write timestamped input file
# filter and validate imported data---------------------------------------------
if(exists("treatments")){
  plot_data <- plot_data |> filter_validate_reorder("treatment", treatments)
}
if(exists("targets")){
  plot_data <- plot_data |> filter_validate_reorder("target", targets)
}
# report raw data and parameters-----------------------------------------------
write_csv(plot_data,
          fs::path(output_directory,
                   str_glue("selectscreen_raw_data_{get_timestamp()}.csv")))
doseplotr::file_copy_to_dir(params_path, output_directory)
doseplotr::file_copy_to_dir(scales_path, output_directory)
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
# fit models to output EC values------------------------------------------------
model_summary <- summarize_models(plot_data,
                                  response_col = "response",
                                  rigid = rigid) |> # use global rigid parameter
  select(-model) |>  # remove actual model from report
  mutate(across(where(is.numeric), \(x){signif(x, digits = 4)}))
write_csv(model_summary,
          str_glue("output/selectscreen_model_summary_{get_timestamp()}.csv"))
# plot data for each treatment separately---------------------------------------
for (trt in treatments){
  # get all targets for this treatment to set legend length
  trt_targets <- as.vector(unique((plot_data |>
                                     filter_trt_tgt(trt = trt))$target))
  if(manually_relabel_targets){
    trt_targets <- Vectorize(get_display_name, vectorize.args = "name")(
      trt_targets, display_names_targets, TRUE)}
  plot_treatment(plot_data, trt, rigid = rigid, grid = grid,
                 no_legend = no_legend,
                 color_map = get_if(color_map_targets,
                                    manually_recolor_targets),
                 shape_map = get_if(shape_map_targets,
                                    manually_reshape_targets),
                 x_limits = get_if(x_limits, global_x_lim),
                 response_col = "response",
                 ylab = "kinase activity (% of untreated)",
                 legend_title = "cell line",
                 legend_labels = get_if(display_names_targets,
                                        manually_relabel_targets,
                                        otherwise = ggplot2::waiver()),
                 plot_title = doseplotr::get_display_name(trt,
                                                          display_names_treatments,
                                                          manually_relabel_treatments)
  ) |> 
    save_plot(
      str_glue("output/selectscreen_treatment_{trt}_{get_timestamp()}.{plot_type}"),
      legend_len = longest(trt_targets))
}  
# plot data for each target separately------------------------------------------
for (tgt in targets){ 
  # get all treatments for this target to set legend length
  tgt_treatments <- as.vector(unique((plot_data |>
                                        filter_trt_tgt(tgt = tgt))$treatment))
  if(manually_relabel_treatments){
    tgt_treatments <- Vectorize(get_display_name, vectorize.args = "name")(
      tgt_treatments, display_names_treatments, TRUE)}
  plot_target(plot_data, tgt,
              rigid = rigid, # global rigid low-dose asymptote parameter
              grid = grid, # global grid plotting parameter
              no_legend = no_legend,
              x_limits = get_if(x_limits, global_x_lim),
              response_col = "response", # selectscreen uses response
              ylab = "kinase activity (% of untreated)",
              legend_title = "treatment",
              legend_labels = get_if(display_names_treatments,
                                     manually_relabel_treatments,
                                     otherwise = ggplot2::waiver()),
              plot_title = get_display_name(tgt,
                                            display_names_targets,
                                            manually_relabel_targets),
              color_map = get_if(color_map_treatments,
                                 manually_recolor_treatments),
              shape_map = get_if(shape_map_treatments,
                                 manually_reshape_treatments)
  )|>
    save_plot(
      str_glue("output/selectscreen_target_{tgt}_{get_timestamp()}.{plot_type}"),
      legend_len = longest(tgt_treatments))
}
