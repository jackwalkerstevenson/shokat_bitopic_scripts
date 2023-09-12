#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2023"
#' ---
#'plateplotr generates dose-response plots from plate reader data.
#' 
#'Input: platemap spreadsheets, Excel or CSV, formatted as per the plater package
#'
#'An import template is provided to make it easy to make platemaps.
#' 
#'Variables expected in input platemaps:
#' 
#'- 'treatment': name of treatment used
#'
#'note: if a dilution series includes a vehicle (zero-concentration) well, it should be labeled as the same treatment
#'  
#'- 'conc_nM' or 'conc_uM': concentration of treatment used in µM or nM
#'- 'target': target of treatment, e.g. cell line or purified protein
#'- 'read_norm': normalized plate reader data (calculated in the platemap)
#'- 'replicate': replicate of condition (included for possible future QC)
#'
#'How to use plateplotr:
#'
#'1. Fill in "treatments.R" with the list of treatments you want to plot
#'2. Make a copy of the import platemap for each plate of data you want to import
#'3. Fill out each platemap with treatment(s), target(s) (e.g. cell line) and concentrations used
#'4. Copy the corresponding raw plate reader data into each platemap
#'5. Copy input platemaps into a directory called "input" in the same directory as this script
#'6. Run plot.R
#'7. Output plots will be created in an "output" folder in the working directory

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(readxl) # for excel file handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis transforms
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(patchwork) # for plot organization
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
  select(-model) |>  # remove actual model from report
  mutate(across(where(is.numeric), \(x){signif(x, digits = 4)}))
write_csv(model_summary,
          str_glue("output/CTG_model_summary_{get_timestamp()}.csv"))
# plot untreated data by target for QC------------------------------------------
if(plate_data |> filter(log_dose == -Inf) |> nrow() > 0){ # only if there are untreated wells
  # set color parameters for treatments
  vr <- viridis_range(length(treatments))
  vr_begin <- vr[[1]]
  vr_end <- vr[[2]]
  vr_option <- vr[[3]] # color scale
  untreated_data_by_tgt <-  plate_data |>
    filter(log_dose == -Inf) |> # untreated data only
    group_by(target) |>
    summarize_response(response_col = "response")
  untreated_data_by_tgt_and_trt <- plate_data |>
    filter(log_dose == -Inf) |> # untreated data only
    group_by(target, treatment) |>
    summarize_response(response_col = "response")
  p <- untreated_data_by_tgt_and_trt |>
    ggplot(aes(y = target, x = mean_response, color = treatment)) +
    geom_point(size = pt_size, alpha = 0.8) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(labels = \(x) format(x, scientific = TRUE)) +
    {if(!manually_recolor_treatments){ # optionally manually specify colors
      viridis::scale_color_viridis(discrete = TRUE, option = vr_option,
                                   begin = vr_begin, end = vr_end)
    } else{
      ggplot2::scale_color_manual(values = color_map_treatments)
    }} +
    geom_errorbar(aes(xmax = mean_response+sem, xmin = mean_response-sem,
                      width = w)) +
    theme_prism() +
    theme(plot.background = element_blank()) + # need for transparent background
    labs(x = "CTG luminescense",
         title = "Growth of untreated control wells") +
    geom_point(data = untreated_data_by_tgt,
             size = 6, alpha = 0.3, color = "black") +
    geom_text(data = untreated_data_by_tgt, # label average points with value
              aes(label = mean_response |>
                    signif(digits = 3) |>
                    format(scientific = TRUE)),
              color = "black",
              show.legend = FALSE, # no legend for the text labels
              position = position_nudge(y = .3)) +
    guides(color = guide_legend(override.aes = list(size = pt_size)))
  save_plot(p, str_glue("output/CTG_QC_untreated_{get_timestamp()}.{plot_type}"),
            width = 14, height = 8)
}
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
                 response_col = "response_norm", # CTG uses response_norm
                 ylab = "luminescence (% of untreated)", # CTG = luminescence assay
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
      str_glue("output/CTG_treatment_{trt}_{get_timestamp()}.{plot_type}"),
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
              response_col = "response_norm", # CTG uses response_norm
              ylab = "luminescence (% of untreated)", # CTG = luminescence assay
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
      str_glue("output/CTG_target_{tgt}_{get_timestamp()}.{plot_type}"),
      no_legend = no_legend,
      legend_len = longest(c(tgt_treatments, treatment_legend_title)))
}
