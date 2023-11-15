# compare values across different experiments
# Jack Stevenson started 2023-10-04
# clear environment and import global parameters--------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_compare_longitudinal.R")
# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(patchwork) # for plot organization
library(assertthat) # for QC assertions
library(doseplotr) # you bet
library(scales) # for plot scales
# import data and set global variables------------------------------------------
all_data <- readxl::read_excel(input_filename)
# filter and validate for data of interest--------------------------------------
plot_data <- all_data |>
  filter_validate_reorder("treatment", treatments) |>
  filter_validate_reorder("target", targets)
# plot IC50s across experiments-----------------------------------------------
plot_data |> 
  # group_by(experiment) |> 
  ggplot(aes(x = treatment, y = IC50_nM, color = experiment)) +
  scale_y_continuous(trans = c("log10", "reverse")) +
  geom_point() +
  geom_line()
  # theme_prism()
# factor test-----------------------------------------------
test_data <- plot_data |> mutate(experiment = as.factor(experiment))
test_data |> 
  # group_by(experiment) |>
  ggplot(aes(x = treatment, y = IC50_nM, color = experiment)) +
  # ggplot(aes(x = treatment, y = IC50_nM)) +
  # scale_y_continuous(trans = c("log10", "reverse")) +
  geom_point() +
  geom_line()
  # geom_line(group = experiment)
# theme_prism()
# manually calculate ratio between treatments-----------------------------------
reference_IC50s <- plot_data |>
  filter(treatment == reference_treatment) |> 
  group_by(experiment) |> 
  summarize(ref_IC50_nM = IC50_nM)
plot_data <- plot_data |> 
  left_join(reference_IC50s, by = "experiment") |> 
  filter(treatment != reference_treatment) |>
  mutate(fold_vs_ref_IC50 = IC50_nM/ref_IC50_nM)