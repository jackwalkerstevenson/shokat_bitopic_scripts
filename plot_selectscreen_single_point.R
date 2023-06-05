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
library(doseplotr)
options(dplyr.summarise.inform = FALSE)
# set global variables----------------------------------------------------------
source("parameters/treatments.R")
source("scatter_plot.R")
dir.create("output/", showWarnings = FALSE) # silently create output directory
plot_type <- "pdf" # file type for saved output plots
input_filename <- "ZLYTE_single_point.csv"
all_data <- import_selectscreen(input_filename) |>
  filter_trt_tgt(trt = treatments)
data_conc_labeled <- all_data |> 
  mutate(target = str_glue("{target} ({Compound.Conc} nM cpd)"))
# helper summary function-------------------------------------------------------
inhibition_summarize <- function(x){
  x |> group_by(target) |>
    mutate(bar_size = .1 * n()) |>
    # keep variables not in current group
    group_by(target, treatment, Compound.Conc, bar_size) |>
    summarize(# standard error for error bars = standard deviation / square root of n
      sem = sd(pct_inhibition, na.rm = TRUE)/sqrt(n()),
      # mean value for plotting
      mean_pct_inhibition = mean(pct_inhibition),
    )}
# helper plot function----------------------------------------------------------
plot_raster <- function(plot){
  plot +
    theme_prism() +
    # reinstitute legend title dropped by theme_prism
    theme(legend.title = element_text()) +
    # remove space from side of table
    scale_x_discrete(expand = expansion(mult = c(0, 0.05))) +
    scale_y_discrete(limits = rev, # reverse y axis to put first on top
                     # remove space from bottom of table
                     expand = expansion(mult = c(0, 0))) +
    scale_fill_viridis(option = "magma",
                       begin = .25, end = 1,
                       limits = c(-5,100)) +
    geom_tile(aes(fill = mean_pct_inhibition)) +
    geom_text()
}
# raster plot for each treatment, all concs-------------------------------------
for(t in treatments){
  trt_data <- all_data |>
    filter_trt_tgt(trt = t) |> 
    inhibition_summarize()
    ggplot(trt_data, aes(x = factor(Compound.Conc),
               y = target,
               label = mean_pct_inhibition)) |> 
      plot_raster() +
      labs(x = "[compound] (nM)",
         y = "target kinase",
         title = str_glue("Single-point SelectScreen potency of {t}"),
         fill = "percent inhibition")
  ggsave(str_glue("output/single_pt_raster_all_targets_all_concs_{t}.pdf"),
         width = 10, height = 10)
  # raster plot for each conc of the treatment, less labeled
  trt_concs <- unique(trt_data$Compound.Conc) # list of concs for this treatment
  for(conc in trt_concs){
    trt_data |>
      filter(Compound.Conc == conc) |> 
      ggplot(aes(x = factor(Compound.Conc),
                 y = target,
                 label = mean_pct_inhibition)) |> 
      plot_raster() +
      labs(x = str_glue("{t}, {conc} nM"),
           y = "target kinase",
           title = str_glue("Single-point SelectScreen potency of {t}"),
           fill = "percent inhibition")
    ggsave(str_glue("output/single_pt_raster_all_targets_{conc}_nM_{t}.pdf"),
           width = 8, height = 10)
    }
}
# aesthetic parameters for scatter plots----------------------------------------
font_base_size <- 14
text_factor <- font_base_size / 130 # assume font base size 14
# global axis labels
xlab = "percent inhibition"
ylab = "target kinase"
vr <- viridis_range(length(treatments))
vr_begin <- vr[1]
vr_end <- vr[2]
# scatter plot for all targets for each treatment/conc--------------------------
for(t in treatments){
  trt_data <- all_data |> 
    filter(treatment == t)
  trt_concs <- unique(trt_data$Compound.Conc) # list of concs for this treatment
  for(conc in trt_concs){
    trt_data |> 
      filter(Compound.Conc == conc) |> 
      scatter_plot(file_name =
                     str_glue("single_pt_scatter_all_targets_{t}_{conc}_nM"),
                   title = "Single-point SelectScreen inhibition",
                   # caption = "Note: compound concentrations not equal between kinases",
                   xlab = xlab, ylab = ylab,
                   viridis_begin = vr_begin, viridis_end = vr_end,
                   width = 10, height = 9, pt_size = 4,)
  }
}
# scatter plots for all treatments and targets-----------------------------------
# assuming same conc per target for all treatments
data_conc_labeled |> 
  scatter_plot(file_name =
                 str_glue("single_pt_scatter_all_targets_all_treatments"),
               title = "Single-point SelectScreen inhibition",
               caption = "Note: compound concentrations not equal between kinases",
               xlab = xlab, ylab = ylab,
               viridis_begin = vr_begin, viridis_end = vr_end,
               width = 10, height = 9, pt_size = 4)
# # individual scatter plots by target---------------------------------
# targets <- unique(all_data$target)
# for(t in all_targets){
#   for(conc in concs){
#     text_width <- text_factor * str_length(t)
#     print(str_glue("plotting single target {t}"))
#     all_data |>
#       filter(target == t) |>
#       filter(Compound.Conc == conc) |> 
#       scatter_plot(file_name = str_glue("target_scatter_plot_{t}_{conc}_nM"),
#                    title = "Single-point SelectScreen inhibition",
#                    caption = "Note: compound concentrations not equal between kinases",
#                    xlab = xlab, ylab = ylab,
#                    viridis_begin = vr_begin, viridis_end = vr_end,
#                    width = 7 + text_width, height = 3,
#                    pt_size = 5, alpha = 0.7)
#   }
# 
# }