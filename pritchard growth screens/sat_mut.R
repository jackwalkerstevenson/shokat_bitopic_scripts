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
all_data <- readr::read_csv(str_glue("{input_directory}/{input_filename}"))
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("pritchard growth screens/plot_sat_mut.R", output_directory)
# write timestamped input files to output
doseplotr::file_copy_to_dir(fs::path(input_directory,input_filename), output_directory)
# write timestamped raw data to output
readr::write_csv(all_data, fs::path(output_directory, str_glue("all_data_{get_timestamp()}.csv")))
# dataframe of manual annotations for sequence position plots-------------------
manual_highlight_residues <- tibble(
  protein_start = c(255, 315, 468),
  residue_name = c("E255", "T315", "V468")
)
# ggplot chunk for aesthetics of sequence plots---------------------------------
seq_plot <- function(){
  list(
    scale_x_continuous(breaks = scales::breaks_width(25),
                       minor_breaks = scales::breaks_width(5),
                       expand = expansion(mult = 0)),
    theme_prism(),
    theme(panel.grid = element_line(color = "black", linewidth = 0.5),
          panel.grid.minor = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = "dotted"),
          plot.background = element_blank()) # for transparent background
  )
}
# test summary by position-----------------------------------------------
summary_data <- all_data |>
  dplyr::group_by(pos.aa) |> 
  dplyr::summarize(
    num_rows = n(),
    total_variant_proportion = sum(variant_proportion)
  )
# plot number of variants per position------------------------------------------
(all_data |>
   count(pos.aa, alt_aa, name = "num_variants_per_position") |>
   ggplot(aes(x = pos.aa, y = num_variants_per_position)) +
   geom_bar(stat = "identity") +
   geom_hline(yintercept = 19, linetype = "dotted") + # line of max coverage
   labs(title = "variant coverage",
        x = "amino acid position",
        y = "number of unique variants") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1),
                      expand = expansion(mult = c(0,0.05))) +
   seq_plot()) |>
  save_plot(str_glue("output/num_variants_{get_timestamp()}.{plot_type}"),
            width = 12, height = 6)
# plot all variants at each position proportionally-----------------------------
(all_data |>
   ggplot(aes(x = pos.aa, y = variant_proportion, fill = alt_aa)) +
   geom_bar(stat = "identity", position = "stack") +
   labs(title = "variant proportions",
        x = "amino acid position",
        y = "proportion of amino acids at position") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1),
                      expand = expansion(mult = 0)) +
   scale_fill_viridis(discrete = TRUE,
                      option = "turbo") +
   # scale_fill_discrete(type = pals::glasbey()) +
   seq_plot()) |>
  save_plot(str_glue("output/variant_proportion_{get_timestamp()}.{plot_type}"),
            width = 12, height = 6)