#' ---
#'title: "plot MTD"
#'author: "Jack Stevenson"
#'date: "2024-09-11"
#' ---
#'started 2024-09-11
#'this script plots animal body weight measurements from an MTD study
#'
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(pals) # for color palettes
library(scales) # for nice breaks
library(doseplotr) # you bet
# set up-----------------------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_MTD.R"
source(params_path)
dir.create("output/", showWarnings = FALSE) # silently create output directory
# import and process data-------------------------------------------------------
animal_key <- readxl::read_excel(key_path)
all_data <- readxl::read_excel(input_path) |> 
  # tidy up wide data
  tidyr::pivot_longer(cols = where(is.numeric),
               names_to = "day",
               values_to = "BW_change_proportion") |> 
  dplyr::mutate(day = as.numeric(day),
                BW_change_percent = 100 * BW_change_proportion) |> 
  # add info about dose by animal
  left_join(animal_key, by = "animal") |> 
  # change dose to factor
  dplyr::mutate(dose_mg_kg = as.factor(dose_mg_kg))
  # todo: assign rank order within groups by ending BW
  # in order to semirandomly assign line types to differentiate individuals
  # dplyr::group_by(dose_mg_kg, animal) |> 
  # dplyr::mutate(rank_order = row_number()
# line plot of BW by animal-----------------------------------------------------
all_data |> 
  ggplot(aes(x = day, y = BW_change_percent, color = dose_mg_kg)) +
  geom_line(aes(group = animal), alpha = 0.3) +
  geom_point(size = 2, alpha = 0.7) +
  theme_prism() +
  scale_x_continuous(breaks = scales::breaks_width(7)) +
  scale_y_continuous(limits = c(-25, 25), breaks = scales::breaks_width(5)) +
  scale_color_manual(values = pals::brewer.blues(n = 4)[2:4]) +
  theme(legend.title = element_text()) +
  # todo: only horizontal grid
  theme(panel.grid = element_line(color = "black",
                                  linetype = "dotted",
                                  linewidth = 0.2),
        plot.background = element_blank()) +
  labs(y = "body weight change (%)",
       color = "dose (mg/kg)",
       title = str_glue("Body weight change in {study_name}"))
ggsave(
  str_glue("{output_dir}/MTD_BW_{study_name}_{get_timestamp()}.{plot_type}"),
  width = 10, height = 6)
