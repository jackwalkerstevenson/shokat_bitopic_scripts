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

# prepare global variables---------------------------------
# the order of the treatment list is the order they will be plotted
source("parameters/treatments.R") # import list of treatments to include in plots
source("parameters/targets.R") # import list of targets to include in plots
# create input and output directories, since git doesn't track empty directories
dir.create("input/", showWarnings = FALSE) # do nothing if directory already exists
dir.create("output/", showWarnings = FALSE)
input_directory <- "input/"
plot_type <- "pdf"
no_legend <- FALSE # global variable for removing all legends from plots
# import data----------------------------------------------------------------
input_filename <- "input/2023-06-21 Ivan raw data names edited.csv"
plate_data <- readr::read_csv(input_filename) |>
  preprocess_plate_data() |>
# plate_data <- import_plates(input_directory) |>
  filter(treatment %in% treatments) |> # take only specified treatments
  filter(target %in% targets) # take only specified targets
# assert that all treatments listed are actually present in imported data
imported_treatments <- unique(plate_data$treatment)
for(treatment in treatments){
  assert_that(treatment %in% imported_treatments,
  msg = glue::glue("treatment '{treatment}' from the list of treatments to plot was not found in imported data"))
}
plate_data <- plate_data |> 
  mutate(target = fct_relevel(target, targets)) |>
  mutate(treatment = fct_relevel(treatment, treatments)) # treatments in order of input list
plot_data <- plate_data |> filter(log_dose != -Inf)
# generate global parameters for all plots------------------------------------------
if (is.null(targets)){
  targets <- as.vector(unique(plot_data$target))}
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plot_data$log_dose))
x_max <- ceiling(max(plot_data$log_dose))
x_limits <- c(x_min, x_max)
global_x_lim <- TRUE # whether to use these global x limits
# x_limits <- c(-11,-5) # manual x limit backup
# create logistic minor breaks for all treatments
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# set default font size for plots
font_base_size <- 14 # 14 is theme_prism default
# set default point size for plots
pt_size = 3
# fit models and output model parameters----------------------------------------
model_summary <- summarize_models(plot_data, response_col = "response_norm")
write_csv(model_summary |> select(-model), # remove actual model from report
          str_glue("output/plate_model_summary_{get_timestamp()}.csv"))
# plot untreated data by target for QC------------------------------------------
# set color parameters for treatments
color_scale <- "viridis"
if(length(treatments) > 6){
  color_scale <- "turbo"
}
vr <- viridis_range(length(treatments))
vr_begin <- vr[[1]]
vr_end <- vr[[2]]
vr_option <- vr[[3]]
pt_size = 3
p <- plate_data |>
  filter(log_dose == -Inf) |>
  group_by(target, treatment) |>
  summarize_response(response_col = "response") |>
  ggplot(aes(y = target, x = mean_response, color = treatment)) +
  geom_point(size = pt_size, alpha = 0.8) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(labels = \(x) format(x, scientific = TRUE)) +
  scale_color_viridis(option = vr_option, discrete = TRUE,
                      begin = vr_begin, end = vr_end) +
  geom_errorbar(aes(xmax = mean_response+sem, xmin = mean_response-sem,
                    width = w)) +
  theme_prism() +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "CTG luminescense",
       title = "Growth of untreated control wells") +
  geom_point(data = plate_data |>
           filter(log_dose == -Inf) |>
           group_by(target) |>
           summarize_response(response_col = "response"),
           size = 6, alpha = 0.3, color = "black") +
  guides(color = guide_legend(override.aes = list(size = pt_size)))
save_plot(p, str_glue("output/plate_QC_untreated_{get_timestamp()}.{plot_type}"),
          width = 14, height = 8)
# plot data for each treatment separately----------------------------------------
for (trt in treatments){
  plot_treatment(plot_data, trt,
                 if(global_x_lim){x_limits = x_limits},
                 response_col = "response_norm") |>
    save_plot(
      str_glue("output/plate_treatment_{trt}_{get_timestamp()}.{plot_type}"),
      legend_len = longest(targets))
}  
# plot data for each target separately------------------------------------------
for (tgt in targets){
  tgt_treatments <- as.vector(unique((plot_data |> filter_trt_tgt(tgt = tgt))$treatment))
  # tgt_treatments_test <- unique(plot_data$treatment)
  plot_target(plot_data, tgt,
                 if(global_x_lim){x_limits = x_limits},
                 response_col = "response_norm") |>
    save_plot(
      str_glue("output/plate_target_{tgt}_{get_timestamp()}.{plot_type}"),
      # legend_len = longest(treatments))
      legend_len = longest(tgt_treatments))
}