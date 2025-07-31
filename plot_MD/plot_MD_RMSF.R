# plot atomwise RMSF data from bitopic MD experiments
# Jack Stevenson started 2025-03-31
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MD_RMSF.R"
source(params_path)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MD/plot_MD_RMSF.R", output_dir)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)

IC50_path <- str_glue("{input_dir}{IC50_filename}")
# write timestamped IC50 file to output
doseplotr::file_copy_to_dir(IC50_path, output_dir)

key_path <- str_glue("{input_dir}{key_filename}")
# write timestamped key file to output
doseplotr::file_copy_to_dir(key_path, output_dir)
# import, preprocess and report data-----------------------------------------------
IC50_data <- readxl::read_excel(IC50_path) |> # import previously measured IC50 data
  dplyr::mutate(
    IC50_nM = as.numeric(IC50_nM)
  )
IC50_ABL1_data <- IC50_data |> # filter data for just ABL1 wt
  dplyr::filter(
    variant == "ABL1 wt",
    assay %in% c("SelectScreen", "Kinomescan")
  )

# function to read atomwise RMSF data files and get info from filenames
read_rmsf_file <- function(file_path) {
  # extract info from filenames of the form "m19_run1_fluct_lig.apf"
  matches <- stringr::str_match(basename(file_path), "^(m\\d+)_run(\\d+)_fluct_lig\\.apf$")
  kenneth_id <- matches[,2]  # "m19"
  run <- matches[,3]       # "1" 
  readr::read_table(file_path, comment = "#", col_names = FALSE) |> 
    dplyr::rename(atom = X1, rmsf = X2) |>  # replace default read_table colnames
    # add columns with data extracted from filename
    dplyr::mutate(kenneth_id = kenneth_id,
                  run = run)
}

# find all files ending in "_fluct_lig.apf" and read 
rmsf_files <- list.files(path = input_dir, pattern = "*_fluct_lig\\.apf$", full.names = TRUE)
all_data <- rmsf_files |> 
  purrr::map(read_rmsf_file) |> # read each file into a dataframe, making a list
  purrr::list_rbind()  # bind the list by rows into a single big dataframe

key_data <- readxl::read_excel(key_path) |> 
  dplyr::mutate(orthosteric_end_linker_atom_num = as.numeric(orthosteric_end_linker_atom_num),
                allosteric_end_linker_atom_num = as.numeric(allosteric_end_linker_atom_num))
# add compound data from key
all_data <- all_data |> 
  dplyr::left_join(key_data, by = join_by(kenneth_id)) |> 
  doseplotr::filter_validate_reorder("compound_name_full", compounds) |> 
  # filter for atoms in linker range. use pmin to get min/max for current row
  dplyr::filter(atom >= pmin(orthosteric_end_linker_atom_num,
                            allosteric_end_linker_atom_num) & atom <= pmax(allosteric_end_linker_atom_num,
                                                                          orthosteric_end_linker_atom_num))

# true/false reversed linker numbering

atomwise_data <- all_data |>
  dplyr::mutate(
    linker_length_atoms =  abs(allosteric_end_linker_atom_num - orthosteric_end_linker_atom_num),
    middle_relative_linker_atom = linker_length_atoms / 2 |> floor(),
    linker_atom_num = case_when(
      !reversed_linker_numbering ~ atom - orthosteric_end_linker_atom_num,
      reversed_linker_numbering ~ orthosteric_end_linker_atom_num - atom),
    relative_linker_atom_num = linker_atom_num - middle_relative_linker_atom,
    normalized_linker_atom_num = linker_atom_num / linker_length_atoms
  )
compoundwise_data <- atomwise_data |>
  dplyr::group_by(compound_name_full) |>
  dplyr::summarize(
    mean_linker_rmsf = mean(rmsf),
    linker_length_atoms = mean(linker_length_atoms) # dumb workaround
  ) |>
  dplyr::left_join(IC50_ABL1_data,
                   by = join_by("compound_name_full" == "treatment")) |> 
  doseplotr::filter_validate_reorder("compound_name_full", compounds)

# report processed data
write_csv(atomwise_data,
          fs::path(output_dir,
                   str_glue("MD_RMSF_atomwise_data_{get_timestamp()}.csv")))
write_csv(compoundwise_data,
          fs::path(output_dir,
                   str_glue("compoundwise_data_{get_timestamp()}.csv")))
# plot atomwise fluctuation, centered atom position----------------------------
atomwise_data |> 
  ggplot(aes(x = relative_linker_atom_num, y = rmsf, color = compound_name_full)) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    width = 0.9,
    alpha = 0.2) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 2) +
  stat_summary(
    fun = "mean",
    geom = "line") +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  theme(plot.background = element_blank()) + # transparent
  labs(x = "relative linker atom position",
       y = "root mean square fluctuation (Å)",
       title ="Linker fluctuation (heavy atoms)")
ggsave(str_glue(
  "{output_dir}/atomwise_RMSF_relative_position_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
# plot atomwise fluctuation, normalized atom position----------------------------
atomwise_data |> 
  ggplot(aes(x = normalized_linker_atom_num, y = rmsf, color = compound_name_full)) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    width = 0.01,
    alpha = 0.2) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 2) +
  stat_summary(
    fun = "mean",
    geom = "line") +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  theme(plot.background = element_blank()) + # transparent
  labs(x = "normalized linker atom position",
       y = "root mean square fluctuation (Å)",
       title ="Linker fluctuation (heavy atoms)")
ggsave(str_glue(
  "{output_dir}/atomwise_RMSF_normalized_position_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
# plot mean RMSF vs selectscreen IC50-----------------------------------------------
compoundwise_data |> 
  dplyr::filter(assay == "SelectScreen") |> 
  ggplot(aes(x = mean_linker_rmsf, y = IC50_nM, color = compound_name_full)) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  labs(
    x = "mean RMSF of linker (heavy atoms)",
    y = "SelectScreen ABL1 IC50 (nM)"
    )
ggsave(str_glue(
  "{output_dir}/RMSF_vs_IC50_SelectScreen_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
# plot mean RMSF vs Kinomescan IC50-----------------------------------------------
compoundwise_data |> 
  dplyr::filter(assay == "Kinomescan") |> 
  ggplot(aes(x = mean_linker_rmsf, y = IC50_nM, color = compound_name_full)) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  labs(
    x = "mean RMSF of linker (heavy atoms)",
    y = "Kinomescan ABL1 IC50 (nM)"
  )
ggsave(str_glue(
  "{output_dir}/RMSF_vs_IC50_Kinomescan_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
# plot mean RMSF vs linker length-----------------------------------------------
compoundwise_data |> 
  ggplot(aes(x = linker_length_atoms, y = mean_linker_rmsf, color = compound_name_full)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  labs(
    x = "linker length (heavy atom count)",
    y = "mean RMSF of linker (heavy atoms)"
  )
ggsave(str_glue(
  "{output_dir}/linker_length_vs_RMSF_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)