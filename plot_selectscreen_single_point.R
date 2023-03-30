#' ---
#'title: "plateplotr selectscreen single point"
#'author: "Jack Stevenson"
#'date: "2023-03-06"
#' ---
#'edited version of plateplotr for dealing with SelectScreen data
#'started 2023-03-06
# load required libraries-----------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(assertthat) # for QC assertions
options(dplyr.summarise.inform = FALSE)
# set global variables----------------------------------------------------------
source("parameters/treatments.R")
source("parameters/targets.R")
source("import_selectscreen.R")
source("viridis_range.R")
source("scatter_plot.R")
dir.create("output/", showWarnings = FALSE) # silently create output directory
plot_type <- "pdf" # file type for saved output plots
input_filename <- "ZLYTE_compiled_results_single_point.csv"
all_data <- import_selectscreen(input_filename, treatments)
# helper summary function-------------------------------------------------------
inhibition_summarize <- function(x){
  x %>% group_by(target) %>%
    summarize(across(c(treatment, pct_inhibition)), # keep variables not in current group
              # scale error bars by size of target group: workaround for dodge issue
              bar_size = .1 * n()) %>%
    group_by(target, treatment) %>%
    summarize(bar_size, # keep variable not in current group
      # standard error for error bars = standard deviation / square root of n
      sem = sd(pct_inhibition, na.rm = TRUE)/sqrt(n()),
      # mean value for plotting
      mean_pct_inhibition = mean(pct_inhibition),
  )}
# aesthetic parameters---------------------------------------------------------
font_base_size <- 14
text_factor <- font_base_size / 130 # assume font base size 14
# global axis labels
xlab = "percent inhibition"
ylab = "target kinase"
vr <- viridis_range(length(treatments))
vr_begin <- vr[1]
vr_end <- vr[2]
print(str_glue("viridis range calculated {vr_begin} to {vr_end}"))
# scatter plot-----------------------------------------------------------------
scatter_plot(all_data, file_name = "single_point_all_targets_scatter",
             title = "Single-point SelectScreen inhibition",
             caption = "Note: compound concentrations not equal between kinases",
             xlab = xlab, ylab = ylab,
             viridis_begin = vr_begin, viridis_end = vr_end,
             width = 10, height = 9, pt_size = 4,)
# multiple scatter plots---------------------------------
all_targets <- distinct(all_data["target"])$target
for(t in all_targets){
  text_width <- text_factor * str_length(t)
  print(str_glue("plotting single target {t}"))
  all_data %>% filter(target == t) %>%
    scatter_plot(file_name = str_glue("target_scatter_plot_{t}"),
                 title = "Single-point SelectScreen inhibition",
                 caption = "Note: compound concentrations not equal between kinases",
                 xlab = xlab, ylab = ylab,
                 viridis_begin = vr_begin, viridis_end = vr_end,
                 width = 7 + text_width, height = 3,
                 pt_size = 5, alpha = 0.7)
}