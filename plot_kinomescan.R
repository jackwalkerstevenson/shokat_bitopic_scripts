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
# set global variables---------------------------------------------------------
input_filename <- "2023-06-27 kinomescan raw data edited for import.csv"
plot_type <- "pdf" # file type for saved output plots
dir.create("output/", showWarnings = FALSE) # silently create output directory
source("parameters/treatments.R")
source("parameters/targets.R")
# import and tidy data---------------------------------
all_data <- import_kinomescan(str_glue("input/{input_filename}")) |>
  filter(log_dose != -Inf) |> # drop 0 values
  filter_trt_tgt(treatments, targets) |> 
  mutate(treatment = fct_relevel(treatment, treatments),
         target = fct_relevel(target, targets),)
# fit models to output EC values------------------------------------------------
EC_summary <- summarize_models(all_data, activity_col = "response_norm")
write_csv(EC_summary, str_glue("output/kinomescan_model_summary_{get_timestamp()}.csv"))
# generate global parameters for all plots------------------------------------------
pt_size = 3 # size for all geom_point
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(all_data$log_dose))
x_max <- ceiling(max(all_data$log_dose))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all conc plots
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# helper function to add ggplot objects common to all plots--------------------
source("dose_response_global.R")
# helper function to plot one treatment----------------------------------------
plot_treatment <- function(data, trt,
                           viridis_begin = 1, viridis_end = 0, viridis_option = "viridis"){
  plate_summary <- data |> 
    filter(treatment == trt) |>  # get data from one treatment to work with
    group_by(target, log_dose) |>   # get set of replicates for each condition
    summarize_response(response_col = "response_norm")
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log_dose, y = mean_response, color = target)) +
      geom_point(aes(shape = target), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_response+sem, ymin = mean_response-sem, width = w)) +
      # use drm method from drc package to fit dose response curve
      #L.4 method is 4-param logistic curve. The more common LL.4() works on nonlog
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1)} |>
    dose_response_global(x_limits) +
    theme(aspect.ratio = 1) +
    scale_color_viridis(discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    labs(title = trt,
         y = "response")
}
# plot data for each treatment separately----------------------------------------
for (trt in treatments){
  trt_data <- filter_trt_tgt(all_data, trt)
  num_tgts <- length(unique(trt_data$target))
  vr <- viridis_range(num_tgts)
  vr_begin <- vr[[1]]
  vr_end <- vr[[2]]
  vr_option <- vr[[3]]
  p <- plot_treatment(trt_data, trt, viridis_option = vr_option,
                      viridis_begin = vr_begin, viridis_end = vr_end)
  # save plot with manually optimized aspect ratio
  save_plot(p, str_glue("output/kinomescan_treatment_{trt}_{get_timestamp()}.{plot_type}"),
            legend_len = longest(targets))
}  
# set color parameters for target plots--------------------------------------
vr <- viridis_range(length(treatments))
vr_begin <- vr[[1]]
vr_end <- vr[[2]]
vr_option <- vr[[3]]
# plot data for each target separately------------------------------------------
for (tgt in targets){
  plate_summary <- all_data |>
    filter(target == tgt) |>
    group_by(treatment, log_dose) |> # group into replicates for each condition
    summarize_response(response_col = "response_norm")
  # bracket ggplot so it can be piped to helper function
  p <- {ggplot(plate_summary, aes(x = log_dose, y = mean_response, color = treatment)) +
      geom_point(aes(shape = treatment), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_response+sem, ymin = mean_response-sem,
                        width = w),
                    alpha = alpha_val) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1, alpha = alpha_val)} |>
    dose_response_global(x_limits) +
    scale_color_viridis(discrete = TRUE, option = vr_option,
                        begin = vr_begin, end = vr_end) +
    labs(title = tgt,
         y = "kinase activity (%)")
  save_plot(p, str_glue("output/kinomescan_target_{tgt}_{get_timestamp()}.{plot_type}"),
            legend_len = longest(treatments))
}
# plot data for all targets at once-----------------------------------------
vr <- viridis_range(length(treatments))
vr_begin <- vr[[1]]
vr_end <- vr[[2]]
vr_option <- vr[[3]]
plate_summary <- all_data |>
  group_by(target, treatment, log_dose) |> # group into replicates for each condition
  summarize_response(response_col = "response_norm")

p <- {ggplot(plate_summary,aes(x = log_dose, y = mean_response, color = treatment)) +
    geom_point(aes(shape = treatment), size = pt_size) +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_response+sem, ymin = mean_response-sem, width = w), alpha = alpha_val) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = target), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, linewidth = 1, alpha = alpha_val)} |>
  dose_response_global(x_limits) +
  scale_color_viridis(option = vr_option, discrete = TRUE,
                      begin = vr_begin, end = vr_end) +
  labs(title = "All data",
       y = "kinase activity (%)")
save_plot(p, str_glue("output/kinomescan_all_data_{get_timestamp()}.{plot_type}"),
          legend_len = longest(treatments))