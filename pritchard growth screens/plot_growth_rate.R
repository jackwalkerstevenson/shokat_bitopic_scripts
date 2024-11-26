#' ---
#'title: "plot growth rate" author: "Jack Stevenson" date: "2024-09-10"
#' ---
#'started 2024-09-10 this script plots results from Pritchard lab growth screen
#'data from 2024-09. its input is data extracted from the environment generated
#'by Josh's script as edited by Jack: 02_growth_tracking_plotting_JWS.R, which
#'runs after 01_growth_tracking_processing_JWS.R, which is unedited from Josh's
#'version
#'

# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(doseplotr) # you bet
library(ggbreak) # for broken plot axis
# set up-----------------------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_growth_rate.R"
source(params_path)
dir.create("output/", showWarnings = FALSE) # silently create output directory
# import and preprocess data-----------------------------------------------
# import all CSVs in input directory
CSV_paths <- get_paths_with_ext(input_dir, "csv")
assertthat::assert_that(length(CSV_paths) > 0,
                        msg = glue::glue("No CSV files found in directory ",
                                         "\"{dir}\""))
# empty dataframe
all_data <- tibble()
# for each file in CSV_paths, read_csv and add to the dataframe
for(file in CSV_paths){
  temp_data <- readr::read_csv(file, name_repair = "universal")
  all_data <- dplyr::bind_rows(all_data, temp_data)
}

all_data <- all_data |> 
  dplyr::rename(all_of(c(treatment = "drug",
                         dose_nM = "dose",
                         theoretical_count = "log_count"))) |> 
  dplyr::mutate(combo_asc_dose_nM = asc_dose_key_nM[library],
                treatment = treatments[treatment],
                relative_growth = theoretical_count / start_live_cell_count) |> 
  # create display names for treatments
  # if combo asc, change treatment name to + asciminib
  # if combo asc, display name dose is original dose + combo asc dose
  dplyr::mutate(treatment = ifelse(combo_asc_dose_nM == 0,
                                   treatment,
                                   str_glue("{treatment} + asciminib"))) |> 
  dplyr::mutate(
    display_name = ifelse(combo_asc_dose_nM == 0,
                          str_glue("{treatment} {dose_nM} nM"),
                          str_glue("{treatment} {dose_nM} nM + {combo_asc_dose_nM} nM")))
  filtered_data <- all_data |> 
    doseplotr::filter_validate_reorder("treatment", treatments) |>
    doseplotr::filter_validate_reorder("display_name", display_names) |> 
    dplyr::filter(theoretical_count < count_plot_limit)
# report input, raw data and parameters-----------------------------------------
write_csv(all_data,
          fs::path(output_dir,
                   str_glue("growth_rate_raw_data_{get_timestamp()}.csv")))
doseplotr::file_copy_to_dir(params_path, output_dir)
doseplotr::file_copy_to_dir(scales_path, output_dir)
doseplotr::file_copy_to_dir("pritchard growth screens/plot_growth_rate.R", output_dir)
# line plot of growth----------------------------------------------------------
filtered_data |>
  dplyr::group_by(display_name, day) |> 
  ggplot(aes(x = day,
             y = relative_growth,
             color = display_name,
             shape = display_name,
             linetype = display_name)) +
  stat_summary(fun.data = "mean_se",
               geom = "errorbar") +
  stat_summary(fun = "mean",
               geom = "line") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 3,
               alpha = 0.7) +
  scale_color_manual(values = color_map_display_names, name = "treatment") +
  scale_shape_manual(values = shape_map_display_names, name = "treatment") +
  scale_linetype_manual(values = linetype_map_display_names, name = "treatment") +
  theme_prism() +
  theme(plot.background = element_blank(),
        legend.title = element_text()) +
  labs(y = "relative cell growth (fold change)",
       title = "Potency against BCR::ABL1 mutant library")

ggsave(str_glue("{output_dir}/sat_mut_growth_rate_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 9, height = 5)
# # working on line plot of growth by treatment with broken y axis---------------------------
# test_treatment <- "ponatinib"
# filtered_data |>
#   dplyr::filter(treatment == test_treatment) |> 
#   dplyr::group_by(display_name, day) |>
#   ggplot(aes(x = day,
#              y = theoretical_count,
#              color = treatment,
#              shape = treatment,
#              linetype = display_name)) +
#   stat_summary(fun.data = "mean_se",
#                geom = "errorbar") +
#   stat_summary(fun = "mean",
#                geom = "line") +
#   stat_summary(fun = "mean",
#                geom = "point",
#                size = 3,
#                alpha = 0.7) +
#   scale_y_break(breaks = c(80, 300),
#                 scales = .5
#                 ) +
#   scale_color_manual(values = color_map_treatments) +
#   scale_shape_manual(values = shape_map_treatments) +
#   scale_linetype_manual(values = linetype_map_display_names) +
#   theme_prism() +
#   guides(color = guide_legend(order = 1),
#          shape = guide_legend(order = 1),
#          linetype = guide_legend(order = 2)) +
#   theme(plot.background = element_blank()) +
#   labs(y = "total theoretical live cell count",
#        title = "Growth of mutant library")
# 
# ggsave(str_glue("{output_dir}/sat_mut_growth_rate_trt_{doseplotr::get_timestamp()}.{plot_type}"),
#        width = 12, height = 7)
