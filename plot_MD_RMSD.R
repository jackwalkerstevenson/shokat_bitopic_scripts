# plot global RMSD data from bitopic MD experiments
# Jack Stevenson started 2025-04-17
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MD_RMSD.R"
source(params_path)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MD_RMSD.R", output_dir)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
read_dat_with_id <- function(file_path) {
  # extract ID from filename
  id <- stringr::str_extract(base::basename(file_path), "^.*(?=_rms\\.dat$)")
  # read data
  data <- readr::read_table(file_path, comment = "#", col_names = FALSE) |> 
    dplyr::rename(frame = X1, rmsd = X2) |> 
    dplyr::mutate(id = id, # add id column extracted from filename
                  time_ns = frame * .004) # Kenneth uses 4 ps = .004 ns per frame
}

dat_files <- list.files(path = input_dir, pattern = "*_rms\\.dat$", full.names = TRUE)
all_data <- dat_files |> 
  purrr::map_dfr(read_dat_with_id)

# report processed data
write_csv(all_data,
          fs::path(output_dir,
                   str_glue("MD_RMSD_all_data{get_timestamp()}.csv")))
# plot RMSD by compound-----------------------------------------------
sparse_sample_factor <- 10
all_data |>
  # sparse sampling
  dplyr::filter(frame %% sparse_sample_factor == 0) |> 
  ggplot(aes(x = time_ns, y = rmsd, color = id)) +
  geom_line(alpha = 0.4) +
  geom_smooth(alpha = 0.7) +
  # scale_color_viridis(option = "turbo") +
  scale_color_manual(values = pals::cols25()) +
  labs(x = "time (ns)",
       y = "global RMSD (Å)",
       title = "RMSD of system from starting position") +
  theme_prism()
ggsave(str_glue(
  "{output_dir}/RMSD_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
