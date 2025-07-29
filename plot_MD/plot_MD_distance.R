# plot ligand distance data from bitopic MD experiments
# Jack Stevenson started 2025-05-29
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MD_distance.R"
source(params_path)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MD/plot_MD_distance.R", output_dir)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
read_distance_file <- function(file_path) {
  # extract info from filenames of the form "m19_run1_dist.aend.dat"
  matches <- stringr::str_match(basename(file_path), "^(m\\d+)_run(\\d+)_(.+)\\.dat$")
  kenneth_id <- matches[,2]  # e.g. "m19"
  run <- matches[,3]       # e.g. "1" 
  ligand <- matches[,4]    # e.g. "ponatinib"
  # read data
  data <- readr::read_table(file_path, comment = "#", col_names = FALSE) |> 
    dplyr::rename(frame = X1, distance = X2) |> 
    dplyr::mutate(kenneth_id = kenneth_id, # add id column extracted from filename
                  run = run,
                  ligand = ligand,
                  time_ns = frame * .004) # Kenneth uses 4 ps = .004 ns per frame
}

dat_files <- list.files(path = input_dir, pattern = "*_dist\\..*\\.dat$", full.names = TRUE)
all_data <- dat_files |> 
  purrr::map_dfr(read_distance_file) |> 
  dplyr::mutate(
    ligand = case_when(
      ligand == "dist.aend" ~ "asciminib",
      ligand == "dist.pend" ~ "ponatinib"
    )
  )

# report processed data
write_csv(all_data,
          fs::path(output_dir,
                   str_glue("MD_distance_all_data_{get_timestamp()}.csv")))
# plot asciminib distance by kenneth_id and run--------------------------------------
sparse_sample_factor <- 20
all_data |>
  # sparse sampling
  dplyr::filter(frame %% sparse_sample_factor == 0,
                ligand == "asciminib") |> 
  ggplot(aes(x = time_ns, y = distance, color = kenneth_id, linetype = run)) +
  geom_line(alpha = 0.4) +
  geom_smooth(alpha = 0.7) +
  scale_color_manual(values = pals::cols25()) +
  labs(x = "time (ns)",
       y = "distance from reference point (Å)",
       title = "Position of asciminib in binding site") +
  theme_prism()
ggsave(str_glue(
  "{output_dir}/distance_asc_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
# plot ponatinib distance by kenneth_id and run--------------------------------------
sparse_sample_factor <- 20
all_data |>
  # sparse sampling
  dplyr::filter(frame %% sparse_sample_factor == 0,
                ligand == "ponatinib") |> 
  ggplot(aes(x = time_ns, y = distance, color = kenneth_id, linetype = run)) +
  geom_line(alpha = 0.4) +
  geom_smooth(alpha = 0.7) +
  scale_color_manual(values = pals::cols25()) +
  labs(x = "time (ns)",
       y = "distance from reference point (Å)",
       title = "Position of ponatinib in binding site") +
  theme_prism()
ggsave(str_glue(
  "{output_dir}/distance_pon_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
