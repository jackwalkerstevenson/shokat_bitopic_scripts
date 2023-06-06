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
source("parameters/targets.R")
source("scatter_plot.R")
dir.create("output/", showWarnings = FALSE) # silently create output directory
plot_type <- "pdf" # file type for saved output plots
input_filename <- "ZLYTE_single_point.csv"
all_data <- import_selectscreen(input_filename) |>
  filter_trt_tgt(trt = treatments, tgt = target_list) |> 
  # sort treatments, targets in list order
  mutate(treatment = fct_relevel(treatment, treatments)) |> 
  mutate(target = fct_relevel(target, target_list))
data_conc_labeled <- all_data |> 
  mutate(target = str_glue("{target} ({Compound.Conc} nM cpd)"))
# aesthetic parameters for scatter plots----------------------------------------
font_base_size <- 14
text_factor <- font_base_size / 130 # assume font base size 14
# global axis labels
xlab = "percent inhibition"
ylab = "target kinase"
vr <- viridis_range(length(treatments))
vr_begin <- vr[1]
vr_end <- vr[2]
# helper summary function-------------------------------------------------------
inhibition_summarize <- function(x){
  x |> group_by(treatment, target) |>
    mutate(bar_size = .2 * n()) |>
    ungroup() |> 
    # keep variables not in current group
    group_by(target, treatment, Compound.Conc, bar_size) |>
    summarize(# standard error for error bars = standard deviation / square root of n
      sem = sd(pct_inhibition, na.rm = TRUE)/sqrt(n()),
      # mean value for plotting
      mean_pct_inhibition = mean(pct_inhibition),
    )}
# helper raster plot function----------------------------------------------------------
raster_helper <- function(plot){
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
                       limits = c(-5,103)) +
    geom_tile(aes(fill = mean_pct_inhibition)) +
    geom_text()
}
# raster plots for individual treatments----------------------------------------
# raster plot for each treatment, all concs
for(t in treatments){
  trt_data <- all_data |>
    filter_trt_tgt(trt = t) |> 
    inhibition_summarize()
    ggplot(trt_data, aes(x = factor(Compound.Conc),
               y = target,
               label = mean_pct_inhibition)) |> 
      raster_helper() +
      labs(x = "[compound] (nM)",
         y = "target kinase",
         title = str_glue("Single-point SelectScreen potency of {t}"),
         fill = "percent inhibition")
  ggsave(str_glue("output/single_pt_raster_all_targets_all_concs_{t}.pdf"),
         bg = "transparent", width = 10, height = 10)
  # raster plot for each conc of the treatment, less labeled
  trt_concs <- unique(trt_data$Compound.Conc) # list of concs for this treatment
  for(conc in trt_concs){
    trt_data |>
      filter(Compound.Conc == conc) |> 
      ggplot(aes(x = factor(Compound.Conc),
                 y = target,
                 label = mean_pct_inhibition)) |> 
      raster_helper() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = str_glue("{t}, {conc} nM"),
           y = "target kinase",
           title = str_glue("Single-point SelectScreen potency of {t}"),
           fill = "percent inhibition")
    ggsave(str_glue("output/single_pt_raster_all_targets_{conc}_nM_{t}.pdf"),
           bg = "transparent", width = 8, height = 10)
    }
}
# raster plot for multiple treatments together at one conc each----------------
all_data |> 
  filter(Compound.Conc %in% c(3.1, 14.8)) |> # manually select EC90 concs
  # mutate(treatment = str_glue("{treatment}, {Compound.Conc} nM")) |> 
  inhibition_summarize() |> 
  ggplot(aes(x = treatment, y = target, label = mean_pct_inhibition)) |> 
  raster_helper() +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.05)),
                   # WARNING MANUAL CONCENTRATION ANNOTATION
                   labels = c("PonatiLink-2-7-10" = "PonatiLink-2-7-10, 3.1 nM",
                              "ponatinib + asciminib" =
                                "ponatinib + asciminib, 14.8 nM")) +
  labs(y = "target kinase",
       title = str_glue("SelectScreen potency at Abl1 ~EC90"),
       fill = "percent inhibition")
ggsave("output/single_pt_raster_EC90_comparison.pdf",
         bg = "transparent", width = 10, height = 10)
# bar plot for multiple treatments at one conc each-----------------------------
geom_barwidth <- 0.75
all_data |> 
  filter(Compound.Conc %in% c(3.1, 14.8)) |> # manually select EC90 concs
  # mutate(treatment = str_glue("{treatment}, {Compound.Conc} nM")) |> 
  inhibition_summarize() |> 
  ggplot(aes(y = target,
             x = mean_pct_inhibition,
             fill = treatment,
             label = mean_pct_inhibition)) +
  scale_fill_viridis(discrete = TRUE, begin = vr_begin, end = vr_end,
                     # WARNING MANUAL CONCENTRATION LABELS
                     labels = c("ponatinib + asciminib" =
                                  "ponatinib + asciminib, 14.8 nM",
                                "PonatiLink-2-7-10" =
                                  "PonatiLink-2-7-10, 3.1 nM")) +
  theme_prism() +
  theme(plot.background = element_blank()) + # need for transparent background
  scale_y_discrete(limits = rev) +
  geom_bar(stat = "identity",
           width = geom_barwidth,
           position = position_dodge2(reverse = TRUE)) +
  geom_errorbar(aes(xmin = mean_pct_inhibition - sem,
                    xmax = mean_pct_inhibition + sem),
                width = geom_barwidth,
                position = position_dodge2(reverse = TRUE,
                                           width = geom_barwidth)) +
  # remove extra space around bars
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "target kinase",
       x = "percent inhibition",
       title = str_glue("SelectScreen potency at Abl1 ~EC90"),)
ggsave("output/single_pt_bar_EC90_comparison.pdf",
       bg = "transparent", width = 10, height = 14)
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