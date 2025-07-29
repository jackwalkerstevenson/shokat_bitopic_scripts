# plot protein-only RMSD data from bitopic MD experiments
# Jack Stevenson started 2025-07-02
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
doseplotr::file_copy_to_dir("plot_MD/plot_MD_RMSD.R", output_dir)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
read_rms_file <- function(file_path) {
  # extract info from filenames of the form "m19_run1_rms.pro.dat"
  matches <- stringr::str_match(basename(file_path), "^(m\\d+)_run(\\d+)_rms\\.pro\\.dat$")
  kenneth_id <- matches[,2]  # "m19"
  run <- matches[,3]       # "1" 
  # read data
  data <- readr::read_table(file_path, comment = "#", col_names = FALSE) |> 
    dplyr::rename(frame = X1, rmsd = X2) |> 
    dplyr::mutate(kenneth_id = kenneth_id, # kenneth_id extracted from filename
                  run = run, # run extracted from filename
                  time_ns = frame * .004) # Kenneth uses 4 ps = .004 ns per frame
}

dat_files <- list.files(path = input_dir, pattern = "*_rms\\.pro\\.dat$", full.names = TRUE)
all_data <- dat_files |> 
  purrr::map_dfr(read_rms_file)

# report processed data
write_csv(all_data,
          fs::path(output_dir,
                   str_glue("MD_RMSD_protein_all_data_{get_timestamp()}.csv")))
# count and report number of frames for each run--------------------------------
frame_summary <- all_data |>
  dplyr::count(kenneth_id, run)
write_csv(frame_summary,
          fs::path(output_dir,
                   str_glue("MD_RMSD_protein_frames_{get_timestamp()}.csv")))
# plot RMSD by kenneth_id and run-----------------------------------------------
sparse_sample_factor <- 20
all_data |>
  # sparse sampling
  dplyr::filter(frame %% sparse_sample_factor == 0) |> 
  ggplot(aes(x = time_ns, y = rmsd, color = kenneth_id, linetype = run)) +
  geom_line(alpha = 0.4) +
  geom_smooth(alpha = 0.7) +
  # scale_color_viridis(option = "turbo") +
  scale_color_manual(values = pals::cols25()) +
  labs(x = "time (ns)",
       y = "protein RMSD (Å)",
       title = "RMSD of protein from starting position") +
  theme_prism()
ggsave(str_glue(
  "{output_dir}/RMSD_protein_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
