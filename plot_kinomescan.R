#' ---
#'title: "plateplotr kinomescan"
#'author: "Jack Stevenson"
#'date: "2023-07-20"
#' ---
#'edited version of plateplotr for dealing with kinomescan data
#'copied from plot_selectscreen.R 2023-07-20

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(patchwork) # for plot organization
library(assertthat) # for QC assertions
library(doseplotr) # you bet
options(dplyr.summarise.inform = FALSE)
# import global parameters and clear environment--------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_plot_kinomescan.R")
# import and tidy data---------------------------------
all_data <- import_kinomescan(str_glue("input/{input_filename}")) |>
  filter(log_dose != -Inf) |> # drop 0 values
  filter_trt_tgt(treatments, targets) |> 
  mutate(treatment = fct_relevel(treatment, treatments),
         target = fct_relevel(target, targets),)
# filter and validate imported data---------------------------------------------
if(exists("treatments")){
  all_data <- filter_validate_reorder(all_data, "treatment", treatments)
}
if(exists("targets")){
  all_data <- filter_validate_reorder(all_data, "target", targets)
}
plot_data <- all_data |> filter(log_dose != -Inf) # data without untreated
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
          str_glue("output/kinomescan_model_summary_{get_timestamp()}.csv"))
# plot data for each treatment separately---------------------------------------
target_legend_title <- if(override_target_title){
  target_title} else "cell line"
for (trt in treatments){
  # get all targets for this treatment to set legend length
  trt_targets <- as.vector(unique((plot_data |>
                                     filter_trt_tgt(trt = trt))$target))
  if(manually_relabel_targets){
    trt_targets <- Vectorize(get_display_name, vectorize.args = "name")(
      trt_targets, display_names_targets, TRUE)}
  plot_treatment(plot_data, trt,
                 rigid = rigid, # global param: rigid low-dose asymptote
                 grid = grid, # global param: background grid on plot
                 no_legend = no_legend, # global param: whether to omit legend
                 x_limits = get_if(x_limits, global_x_lim),
                 response_col = "response_norm",
                 ylab = "signal (% of untreated)",
                 legend_title = target_legend_title,
                 # if relabeling targets, get display names for legend
                 legend_labels = get_if(display_names_targets,
                                        manually_relabel_targets,
                                        otherwise = ggplot2::waiver()),
                 #if relabeling, get display name for title
                 plot_title = get_display_name(trt,
                                               display_names_treatments,
                                               manually_relabel_treatments),
                 # if manually setting colors of targets, get color map
                 color_map = get_if(color_map_targets,
                                    manually_recolor_targets),
                 # if manually setting shapes of targets, get shape map
                 shape_map = get_if(shape_map_targets,
                                    manually_reshape_targets)
  ) |> 
    save_plot(
      str_glue("output/kinomescan_treatment_{trt}_{get_timestamp()}.{plot_type}"),
      legend_len = if(no_legend) 0 else{
        longest(c(trt_targets, target_legend_title))})
}  
# plot data for each target separately------------------------------------------
treatment_legend_title <- if(override_treatment_title){
  treatment_title} else "treatment"
for (tgt in targets){ 
  # get all treatments for this target to set legend length
  tgt_treatments <- as.vector(unique((plot_data |>
                                        filter_trt_tgt(tgt = tgt))$treatment))
  if(manually_relabel_treatments){
    tgt_treatments <- Vectorize(get_display_name, vectorize.args = "name")(
      tgt_treatments, display_names_treatments, TRUE)}
  plot_target(plot_data, tgt,
              rigid = rigid, # global param: rigid low-dose asymptote
              grid = grid, # global param: background grid on plot
              no_legend = no_legend, # global param: whether to omit legend
              x_limits = get_if(x_limits, global_x_lim),
              response_col = "response_norm",
              ylab = "signal (% of untreated)",
              legend_title = treatment_legend_title,
              # if relabeling treatments, get display names for legend
              legend_labels = get_if(display_names_treatments,
                                     manually_relabel_treatments,
                                     otherwise = ggplot2::waiver()),
              #if relabeling, get display name for title
              plot_title = get_display_name(tgt,
                                            display_names_targets,
                                            manually_relabel_targets),
              # if manually setting colors of treatments, get color map
              color_map = get_if(color_map_treatments,
                                 manually_recolor_treatments),
              # if manually setting shapes of treatments, get shape map
              shape_map = get_if(shape_map_treatments,
                                 manually_reshape_treatments)
  )|>
    save_plot(
      str_glue("output/kinomescan_target_{tgt}_{get_timestamp()}.{plot_type}"),
      no_legend = no_legend,
      legend_len = longest(c(tgt_treatments, treatment_legend_title)))
}
