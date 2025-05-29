# visualizing resistance data from saturation mutagenesis
# Jack Stevenson, started 2023-11-15
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
#library(readxl) # for excel file handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis transforms
library(viridis) # for color schemes
library(pals) # color palettes
library(doseplotr) # you bet
# import global parameters and clear environment--------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_sat_mut.R"
source(params_path)
dir.create(output_directory, showWarnings = FALSE)
# import data and report-----------------------------------------------
raw_data <- readr::read_csv(str_glue("{input_directory}/{input_filename}"))
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("pritchard growth screens/plot_sat_mut.R", output_directory)
# write timestamped input files to output
doseplotr::file_copy_to_dir(fs::path(input_directory,input_filename), output_directory)
# write timestamped raw data to output
readr::write_csv(raw_data, fs::path(output_directory, str_glue("raw_data_{get_timestamp()}.csv")))
# process data-----------------------------------------------
# clearer column names
all_data <- raw_data |> 
  dplyr::rename(
    aa_position = pos.aa,
    plasmid_variant_proportion = variant_proportion,
    count_rep_a = ct.bl.a,
    depth_rep_a = depth.bl.a,
    freq_rep_a = freq.bl.a,
    count_rep_b = ct.bl.b,
    depth_rep_b = depth.bl.b,
    freq_rep_b = freq.bl.b,
  )
# summary stats by aa position
summary_data <- all_data |>
  dplyr::group_by(aa_position) |> 
  dplyr::summarize(
    num_plasmid_variants = n(), # each row represents one variant detected in plasmid library
    total_plasmid_proportion = sum(plasmid_variant_proportion), # for normalizing
    total_count_rep_a = sum(count_rep_a, na.rm = TRUE),
    total_count_rep_b = sum(count_rep_b, na.rm = TRUE)
  )
all_data <- all_data |> 
  dplyr::left_join(summary_data) |> # join summary back to parent data
  dplyr::mutate(
    # normalize plasmid proportions to 1. they're already close, just a little noisy
    normalized_plasmid_proportion = plasmid_variant_proportion / total_plasmid_proportion,
    rep_a_cell_proportion = count_rep_a / total_count_rep_a,
    rep_b_cell_proportion = count_rep_b / total_count_rep_b,
    total_cell_count = total_count_rep_a + total_count_rep_b,
    # this sums to more than 1
    #test_cell_proportion = dplyr::case_when(
    #  !is.na(count_rep_a) & !is.na(count_rep_b) ~ (rep_a_cell_proportion + rep_b_cell_proportion) / 2,
    #  !is.na(count_rep_a) ~ rep_a_cell_proportion,
    #  !is.na(count_rep_b) ~ rep_b_cell_proportion,
    #),
    # just total count over total total
    mean_cell_proportion = rowSums(cbind(count_rep_a, count_rep_b), na.rm = TRUE) / total_cell_count
    )
# summary statistics-----------------------------------------------
detected_seq_data <- all_data |> 
  dplyr::filter(!is.na(count_rep_a) | !is.na(count_rep_b))
first_aa <- all_data$aa_position |> min()
last_aa <- all_data$aa_position |> max()
total_theoretical_variants <- (last_aa - first_aa) * 19
total_plasmid_variants <- all_data$mean_cell_proportion |> length()
total_seq_variants <- detected_seq_data$mean_cell_proportion |> length()
plasmid_vs_theoretical_coverage <- total_plasmid_variants / total_theoretical_variants
seq_vs_theoretical_coverage <- total_seq_variants / total_theoretical_variants
seq_vs_plasmid_coverage <- total_seq_variants / total_plasmid_variants

  
# ggplot chunk for aesthetics of sequence plots---------------------------------
seq_plot <- function(){
  list(
    scale_x_continuous(breaks = scales::breaks_width(25),
                       minor_breaks = scales::breaks_width(5),
                       expand = expansion(mult = 0)),
    theme_prism(),
    theme(#panel.grid = element_line(color = "black", linewidth = 0.5),
          #panel.grid.minor = element_line(color = "black",
          #                                linewidth = 0.1,
          #                                linetype = "dotted"),
          plot.background = element_blank()) # for transparent background
  )
}
# plot number of plasmid variants per position------------------------------------------
(all_data |>
   count(aa_position, alt_aa, name = "num_variants_per_position") |>
   ggplot(aes(x = aa_position, y = num_variants_per_position)) +
   geom_bar(stat = "identity", fill = "skyblue3") +
   geom_hline(yintercept = 19, linetype = "dotted", linewidth = 1) + # line of max coverage
   labs(title = "Amino acid variant coverage in plasmid library",
        x = "amino acid position",
        y = "number of amino acid variants") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1),
                      expand = expansion(mult = c(0,0.05))) +
   seq_plot()) |>
  save_plot(str_glue("output/num_plasmid_variants_{get_timestamp()}.{plot_type}"),
            width = 10, height = 4)
# plot number of cell pool variants per position--------------------------------
(all_data |>
   dplyr::filter(mean_cell_proportion > 0) |> 
   count(aa_position, alt_aa, name = "num_variants_per_position") |>
   ggplot(aes(x = aa_position, y = num_variants_per_position)) +
   geom_bar(stat = "identity", fill = "skyblue3") +
   geom_hline(yintercept = 19, linetype = "dotted", linewidth = 1) + # line of max coverage
   labs(title = "Amino acid variant coverage in baseline cell pool",
        x = "amino acid position",
        y = "number of amino acid variants") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1),
                      expand = expansion(mult = c(0,0.05))) +
   seq_plot()) |>
  save_plot(str_glue("output/num_cell_variants_{get_timestamp()}.{plot_type}"),
            width = 10, height = 4)
# plot all plasmid variants at each position proportionally-----------------------------
(all_data |>
   ggplot(aes(x = aa_position, y = normalized_plasmid_proportion, fill = alt_aa)) +
   geom_bar(stat = "identity", position = "stack") +
   labs(title = "variant proportions in plasmid library",
        x = "amino acid position",
        y = "proportion of amino acids at position") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1),
                      expand = expansion(mult = 0)) +
   scale_fill_viridis(discrete = TRUE,
                      option = "turbo") +
   # scale_fill_discrete(type = pals::glasbey()) +
   seq_plot()) |>
  save_plot(str_glue("output/plasmid_variant_proportion_{get_timestamp()}.{plot_type}"),
            width = 10, height = 6)
# plot all cell pool variants at each position proportionally-------------------
(all_data |>
   ggplot(aes(x = aa_position, y = mean_cell_proportion, fill = alt_aa)) +
   geom_bar(stat = "identity", position = "stack") +
   labs(title = "variant proportions in baseline cell pool",
        x = "amino acid position",
        y = "proportion of amino acids at position") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1),
                      expand = expansion(mult = 0)) +
   scale_fill_viridis(discrete = TRUE,
                      option = "turbo") +
   # scale_fill_discrete(type = pals::glasbey()) +
   seq_plot()) |>
  save_plot(str_glue("output/cell_variant_proportion_{get_timestamp()}.{plot_type}"),
            width = 10, height = 5)