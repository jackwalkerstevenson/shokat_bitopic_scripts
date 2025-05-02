# plot data from bitopic MD experiments
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
input_path <- str_glue("{input_dir}{input_filename}")
IC50_path <- str_glue("{input_dir}{IC50_filename}")
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_dir)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_MD_RMSF.R", output_dir)
# write timestamped input file to output
doseplotr::file_copy_to_dir(input_path, output_dir)
# write timestamped IC50 file to output
doseplotr::file_copy_to_dir(IC50_path, output_dir)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
raw_data <- readxl::read_excel(input_path)
if(exists("treatments")){
  raw_data <- filter_validate_reorder(raw_data, "compound_name_full", treatments)
}
IC50_data <- readxl::read_excel(IC50_path) |> 
  dplyr::mutate(
    IC50_nM = as.numeric(IC50_nM)
  )
IC50_ABL1_data <- IC50_data |> 
  dplyr::filter(
    variant == "ABL1 wt",
    assay %in% c("SelectScreen", "Kinomescan")
  )
    
atomwise_data <- raw_data |> 
  dplyr::mutate(
    linker_length_atoms =  allosteric_end_linker_atom_num - orthosteric_end_linker_atom_num,
    middle_linker_atom = (orthosteric_end_linker_atom_num + allosteric_end_linker_atom_num) / 2 |> floor(),
    linker_atom_num = atom_serial_number - orthosteric_end_linker_atom_num,
    relative_linker_atom_num = atom_serial_number - middle_linker_atom,
    normalized_linker_atom_num = linker_atom_num / linker_length_atoms
  )
compoundwise_data <- atomwise_data |> 
  dplyr::group_by(compound_name_full) |> 
  dplyr::summarize(
    mean_linker_rmsf = mean(rmsf),
    linker_length_atoms = mean(linker_length_atoms) # dumb workaround
  ) |> 
  dplyr::left_join(IC50_ABL1_data, by = join_by("compound_name_full" == "treatment"))

# report processed data
write_csv(atomwise_data,
          fs::path(output_dir,
                   str_glue("MD_RMSF_atomwise_data_{get_timestamp()}.csv")))
write_csv(compoundwise_data,
          fs::path(output_dir,
                   str_glue("MD_RMSF_compoundwise_data_{get_timestamp()}.csv")))
# plot atomwise fluctuation, centered atom position----------------------------
atomwise_data |> 
  dplyr::group_by(compound_name_full) |> 
  ggplot(aes(x = relative_linker_atom_num, y = rmsf, color = compound_name_full)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  theme(plot.background = element_blank()) + # transparent
        # legend.title = element_text()) +
  labs(x = "relative linker atom position",
       y = "root mean square fluctuation (Å)",
       title ="Linker fluctuation (heavy atoms)")
ggsave(str_glue(
  "{output_dir}/atomwise_RMSF_relative_position_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = 9, height = 5)
# plot atomwise fluctuation, normalized atom position----------------------------
atomwise_data |> 
  dplyr::group_by(compound_name_full) |> 
  ggplot(aes(x = normalized_linker_atom_num, y = rmsf, color = compound_name_full)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme_prism() +
  theme(plot.background = element_blank()) + # transparent
  # legend.title = element_text()) +
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