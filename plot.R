#'plateplotr generates plots with dose-response curves from 96-well plate data.
#'input: data in plater format for tidy import
#'copy raw plate reader data into JWS Excel template and fill in input cells
#'save input file as csv for import
library(drc)  # for dose response curves
library(tidyverse)
library(ggprism)  # for pretty prism-like plots
library(plater)  # for plate data import
# set up input parameters
# input_filename = "JS-B2-79 CTG 1.xlsx"
# reference_symbol = "<>"  # symbol that appears right before plate data
# rows_to_skip = 1  # rows at beginning of plate without data
# rows_per_condition = 3
# conditions_per_plate = 2
# data_rows = rows_per_condition * conditions_per_plate
# # read input data and slice out just the values
# full_sheet <- read_excel(input_filename, col_names = FALSE)
# ref_row <- which(full_sheet$...1==reference_symbol)
# first_row = ref_row+1+rows_to_skip
# final_row = ref_row+1+data_rows
# plate_data <- full_sheet %>%
#   slice(first_row:final_row) %>%
#   select(3:12)
plater_data <- read_plate("JS-B2-79 plater test.csv") %>%
  filter(drug != "N/A") %>%
  # drop 0 values for plotting and curve fitting
  # note this is only OK because normalization is happening in Excel
  filter(concentration != 0) %>%
  mutate(log.conc = log10(concentration))
plater_data.summary <- plater_data %>%
    group_by(cell_line, concentration) %>%
    summarize(
      sem = sd(CTG_normalized, na.rm = TRUE)/sqrt(n()),
      CTG_normalized = mean(CTG_normalized),
      log.conc = log.conc,
      drug = drug
    )

plater_data.summary %>% ggplot(aes(x = log.conc, y = CTG_normalized, color = cell_line))+
    geom_point()+
    geom_errorbar(aes(ymax = CTG_normalized+sem, ymin = CTG_normalized-sem))+
    geom_smooth(method = "drm", method.args = list(fct = L.4()), se = FALSE)+
    theme_prism()+
    labs(x = "Log [compound] (µM)", y = "Relative cell viability (%)")
