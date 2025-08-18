# plot thermochemistry data from bitopic MD experiments
# depends on output from plot_MD_RMSF
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
doseplotr::file_copy_to_dir("plot_MD/plot_MD_thermo.R", output_dir) # write timestamped code
doseplotr::file_copy_to_dir(params_path, output_dir) # write timestamped params
input_path <- str_glue("{input_dir}{input_filename}")
doseplotr::file_copy_to_dir(input_path, output_dir) # write timestamped input
RMSF_path <- str_glue("{input_dir}{RMSF_filename}")
doseplotr::file_copy_to_dir(RMSF_path, output_dir) # write timestamped RMSF file
IC50_path <- str_glue("{input_dir}{IC50_filename}")
doseplotr::file_copy_to_dir(IC50_path, output_dir) # write timestamped IC50 file
key_path <- str_glue("{input_dir}{key_filename}")
doseplotr::file_copy_to_dir(key_path, output_dir) # write timestamped key file
# import, preprocess and report data-----------------------------------------------
raw_thermo_data <- readxl::read_excel(input_path)
key_data <- readxl::read_excel(key_path)

# import RMSF/IC50 data from other experiments
RMSF_data <- readr::read_csv(RMSF_path) |> 
  dplyr::left_join(key_data, by = join_by("compound_name_full"))

# import IC50 data from other experiments
IC50_data_all <- readxl::read_excel(IC50_path) |>
  dplyr::mutate(
    IC50_nM = as.numeric(IC50_nM),
    linker_length_PEG = as.numeric(linker_length_PEG)
  )
IC50_ABL1_data <- IC50_data_all |> # filter for just ABL1 wt and desired assays
  dplyr::filter(
    variant == "ABL1 wt",
    assay %in% c("SelectScreen", "Kinomescan")
  )

# join runwise thermo to runwise RMSF
runwise_thermo_RMSF_data <- raw_thermo_data |> 
  dplyr::left_join(RMSF_data, by = join_by("kenneth_id", "run")) |> 
  doseplotr::filter_validate_reorder("compound_name_full", compounds) |> 
  dplyr::mutate(entropy_per_mass = entropy / molecular_weight,
                entropy_per_linker_atom = entropy / linker_length_atoms)

# join runwise to IC50 data
thermo_potency_data <- runwise_thermo_RMSF_data |> 
  dplyr::left_join(IC50_ABL1_data, by = join_by("compound_name_full" == "treatment",
                                                "linker_length_PEG" == "linker_length_PEG"))

write_csv(runwise_thermo_RMSF_data, # report processed data
          fs::path(output_dir, str_glue("thermo_RMSF_runwise_data_{get_timestamp()}.csv")))
# set up plot constants-----------------------------------------------
scatter_plot_width <- 7
scatter_plot_height <- 5
# plot run entropy vs linker length-----------------------------------------------
runwise_thermo_RMSF_data |> 
  dplyr::group_by(compound_name_full) |>
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
# plot run entropy vs SelectScreen potency-----------------------------------------------
thermo_potency_data |> 
  dplyr::filter(assay == "SelectScreen") |> 
  ggplot(aes(y = entropy, x = IC50_nM, color = compound_name_full)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    # width = 2,
    alpha = 0.5) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3) +
  scale_x_log10() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  labs(
    y = "compound entropy",
    x = "SelectScreen ABL1 IC50 (nM)"
  ) +
  coord_flip()
ggsave(str_glue(
  "{output_dir}/entropy_vs_IC50_SelectScreen_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)
# plot run entropy per linker length vs Selectscreen potency-----------------
thermo_potency_data |> 
  dplyr::filter(assay == "SelectScreen") |> 
  ggplot(aes(y = entropy_per_linker_atom, x = IC50_nM, color = compound_name_full)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    # width = 2,
    alpha = 0.5) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3) +
  scale_x_log10() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  labs(
    y = "compound entropy per linker atom",
    x = "SelectScreen ABL1 IC50 (nM)"
  ) +
  coord_flip()
ggsave(str_glue(
  "{output_dir}/entropy_per_linker_vs_IC50_SelectScreen_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)
# plot run entropy vs mean linker RMSF-----------------
thermo_potency_data |> 
  dplyr::filter(assay == "SelectScreen") |> 
  ggplot(aes(y = entropy, x = mean_linker_rmsf, color = compound_name_full)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    # width = 2,
    alpha = 0.5) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3) +
  scale_x_log10() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  labs(
    y = "compound entropy per linker atom",
    x = "SelectScreen ABL1 IC50 (nM)"
  ) +
  coord_flip()
ggsave(str_glue(
  "{output_dir}/entropy_per_linker_vs_IC50_SelectScreen_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)