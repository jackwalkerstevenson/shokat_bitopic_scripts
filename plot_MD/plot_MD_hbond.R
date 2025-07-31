# plot hbond data from bitopic MD experiments
# Jack Stevenson started 2025-07-31
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MD_hbond.R"
source(params_path)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MD/plot_MD_hbond.R", output_dir)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# import key file
key_path <- str_glue("{input_dir}{key_filename}")
doseplotr::file_copy_to_dir(key_path, output_dir) # write timestamped file
# import, preprocess and report data-----------------------------------------------
key_data <- readxl::read_excel(key_path)
read_hbond_file <- function(file_path) {
  # extract info from filenames of the form "m19_run1_hbond.dat"
  matches <- stringr::str_match(basename(file_path), "^(m\\d+)_run(\\d+)_hbond_ligprot\\.dat$")
  kenneth_id <- matches[,2]  # "m19"
  run <- matches[,3]       # "1" 
  # read data
  data <- readr::read_table(file_path, comment = "#", col_names = FALSE) |> 
    dplyr::rename(frame = X1, hbonds_lig_to_prot = X2) |> 
    dplyr::mutate(kenneth_id = kenneth_id, # kenneth_id extracted from filename
                  run = run, # run extracted from filename
                  time_ns = frame * .004) # Kenneth uses 4 ps = .004 ns per frame
}

dat_files <- list.files(path = input_dir, pattern = "*_hbond_ligprot\\.dat$", full.names = TRUE)
hbond_data <- dat_files |> 
  purrr::map_dfr(read_hbond_file) |> 
  dplyr::left_join(key_data, by = join_by(kenneth_id))
mean_hbond_data <- hbond_data |> 
  dplyr::group_by(kenneth_id, run) |> 
  dplyr::summarize(mean_hbonds_lig_to_prot = mean(hbonds_lig_to_prot)) |>
  dplyr::left_join(key_data, by = join_by(kenneth_id))

# report processed data
write_csv(mean_hbond_data,
          fs::path(output_dir,
                   str_glue("mean_hbond_data_{get_timestamp()}.csv")))
# set up plot constants-----------------------------------------------
scatter_plot_width <- 7
scatter_plot_height <- 5
# plot mean hbonds vs linker length-----------------------------------------------
mean_hbond_data |> 
  ggplot(aes(x = linker_length_atoms, y = mean_hbonds_lig_to_prot, color = series)) +
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
       y = "ligand-protein H-bonds (time average)",
       title ="H-bonds vs linker length")
ggsave(str_glue(
  "{output_dir}/mean_hbonds_vs_linker_atoms_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)
# plot hbonds over time-----------------------------------------------
sparse_sample_factor <- 500
hbond_data |>
  # sparse sampling
  dplyr::filter(frame %% sparse_sample_factor == 0) |> 
  ggplot(aes(x = time_ns, y = hbonds_lig_to_prot, color = compound_name_full)) +
  # stat_summary(
  #   fun.data = "mean_se",
  #   geom = "errorbar",
  #   width = 2,
  #   alpha = 0.5) +
  # stat_summary(
  #   fun = "mean",
  #   geom = "point",
  #   size = 1) +
  stat_summary(
  fun = "mean",
  geom = "line",
  alpha = 0.6) +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  theme(plot.background = element_blank()) + # transparent
  labs(x = "time (ns)",
       y = "ligand-protein H-bonds",
       title ="H-bonds over time")
ggsave(str_glue(
  "{output_dir}/hbonds_over_time_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)