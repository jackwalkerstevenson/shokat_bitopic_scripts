#'This is the R script that is used to generate nice plots from plate data.
#'Assumes a 96-well plate in which the outer border of cells are not read
#'Assumes the format output by Tecan Spark:
#'"<>" appears at the corner of row and column names
library(tidyverse)
library(readxl)
input_filename = "JS-B2-79 CTG 1.xlsx"
reference_symbol = "<>"  # symbol that appears right before plate data
rows_per_condition = 3
conditions_per_plate = 2
data_rows = rows_per_condition * conditions_per_plate
rows_to_skip = 1  # rows at beginning of plate without data
full_sheet <- read_excel(input_filename, col_names = FALSE)
ref_row <- which(full_sheet$...1==reference_symbol)
plate_data <- slice(full_sheet,ref_row+1+rows_to_skip:ref_row+1+data_rows-1)
