# plot thermochemistry data from bitopic MD experiments
# Jack Stevenson started 2025-07-31
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MD_thermo.R"
source(params_path)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MD/plot_MD_thermo.R", output_dir)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# import input file
input_path <- str_glue("{input_dir}{input_filename}")
doseplotr::file_copy_to_dir(input_path, output_dir) # write timestamped file
# import key file
key_path <- str_glue("{input_dir}{key_filename}")
doseplotr::file_copy_to_dir(key_path, output_dir) # write timestamped file
# import, preprocess and report data-----------------------------------------------
raw_thermo_data <- readxl::read_excel(input_path)
key_data <- readxl::read_excel(key_path)
thermo_data <- raw_thermo_data |>
  dplyr::left_join(key_data, by = join_by(kenneth_id)) |> 
  doseplotr::filter_validate_reorder("compound_name_full", compounds) |> 
  dplyr::mutate(entropy_per_mass = entropy / molecular_weight,
                entropy_per_peg = entropy / linker_length_atoms)
# report processed data
write_csv(thermo_data,
          fs::path(output_dir,
                   str_glue("thermo_data_{get_timestamp()}.csv")))
# set up plot constants-----------------------------------------------
scatter_plot_width <- 7
scatter_plot_height <- 5
# plot entropy vs linker length-----------------------------------------------
thermo_data |> 
  ggplot(aes(x = linker_length_atoms, y = entropy, color = series)) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    width = 2,
    alpha = 0.5) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 2) +
  # stat_summary(
  # fun = "mean",
  # geom = "line",
  # alpha = 0.5) +
  scale_color_manual(values = color_map_series) +
  theme_prism() +
  theme(plot.background = element_blank()) + # transparent
  labs(x = "linker length (atoms)",
       y = "compound entropy (cal/mol-Kelvin)",
       title ="Entropy vs linker length")
ggsave(str_glue(
  "{output_dir}/entropy_vs_linker_atoms_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)