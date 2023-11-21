# JWS visualizing growth of pooled screens. data from Josh Reynolds 2023-11-21
# load libraries----------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(scales) # for plot scales
library(doseplotr) # you bet
library(patchwork) # for compositing plots
# import global parameters and clear environment--------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_screen_growth.R")
# import and preprocess data----------------------------------------------------
# create input and output directories, since git doesn't track empty directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
original_data <- readr::read_csv(input_filename, skip_empty_rows = TRUE)
all_data <- original_data |> 
  filter(!is.na(log_count)) |> # filter for actual counts
  rename(treatment = drug, dose_nM = dose)
  # filter_validate_reorder("treatment", names(concs_to_plot)) |> 
  # filter_validate_reorder("id", manual_id)
# check what conditions exist
test_count <- original_data |> 
  count(drug, dose)
# aesthetic chunk for growth plots----------------------------------------------
growth_plot <- function(){
  list(
    geom_point(size = pt_size),
    geom_line(),
    scale_y_continuous(labels = label_comma(accuracy = 1, big.mark = ",")),
    theme_prism(),
    labs(x= "day",
         y = "calculated cell count")
  )
}
# chunk for manual naming-----------------------------------------------
manual_color_shape <- function(display_names){
  list(
    scale_color_manual(values = color_map_treatments,
                       labels = display_names),
    scale_shape_manual(values = shape_map_treatments,
                         labels = display_names)
  )}
# plot everything--------------------------------------------------------------
ggplot(all_data, aes(x = day, y = log_count, color = id)) +
  growth_plot() +
  scale_color_manual(values = manual_color_four_each) +
  labs(title = "Growth of screen pools")
# function to plot specific concs-----------------------------------------------
plot_concs <- function(data, concs, identifier, log_scale = FALSE){
  filtered_data <- data |> 
    filter(dose_nM == concs[treatment])
  p <- ggplot(filtered_data, aes(x = day, y = log_count,
                                 color = treatment, shape = treatment)) +
    growth_plot() +
    manual_color_shape(display_names_treatments) +
    labs(title = str_glue("Growth of screen pools at {identifier}"))
  if(log_scale){p <- p +
    scale_y_continuous(trans = "log10",
                       labels = label_comma(accuracy = 1, big.mark = ","))}
  scale <- if(log_scale){"log"} else{"linear"}
  save_plot(p,
            str_glue("output/screen_growth_{identifier}_{scale}_{get_timestamp()}.{plot_type}"),
            width = 8, height = 5)
  return(p)
}
# plot each IC set-----------------------------------------------
{plot_concs(all_data, concs_IC25, "IC25") +
plot_concs(all_data, concs_IC50, "IC50") +
plot_concs(all_data, concs_IC75, "IC75") +
plot_concs(all_data, concs_IC90, "IC90")} |> 
  save_plot(str_glue("output/screen_growth_linear_{get_timestamp()}.{plot_type}"),
            width = 16, height = 10)
{plot_concs(all_data, concs_IC25, "IC25", log_scale = TRUE) +
plot_concs(all_data, concs_IC50, "IC50", log_scale = TRUE) +
plot_concs(all_data, concs_IC75, "IC75", log_scale = TRUE) +
plot_concs(all_data, concs_IC90, "IC90", log_scale = TRUE)} |> 
  save_plot(str_glue("output/screen_growth_log_{get_timestamp()}.{plot_type}"),
            width = 16, height = 10)
