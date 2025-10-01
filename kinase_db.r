# interpret data from PKIDB
# Jack Stevenson started 2025-09-30
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr)
library(ggprism) # for prism theme
# set up-----------------------------------------------
rm(list = ls()) # clear environment
input_dir <- "input"
output_dir <- "output"
input_filename <- "pkidb_2025-04-15.tsv"
# write timestamped code to output
doseplotr::file_copy_to_dir("kinase_db.R",output_dir)
# import, preprocess and report data--------------------------------------------
all_data <- readr::read_tsv(fs::path(input_dir, input_filename),
                            na = c("", "NA", "nan"))
data_separate_targets <- all_data |> 
  tidyr::separate_longer_delim(cols = Targets, delim = "; ")

targets <- data_separate_targets |> 
  dplyr::count(Targets)
# plot  first few targets by number of drugs-----------------------------------------------
targets |> 
  tidyr::drop_na(Targets) |> 
  dplyr::arrange(desc(n)) |> 
  dplyr::slice_head(n = 30) |> 
  dplyr::mutate(Targets = forcats::fct_inorder(Targets)) |> 
  ggplot(aes(x = Targets, y = n)) +
  geom_col() +
  labs(x = "kinase",
       y = "number of targets",
       title = "Kinases in the PKIDB by number of targets") +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 270))
ggsave(fs::path(output_dir, str_glue("kinases_by_target_{doseplotr::get_timestamp()}.pdf")),
       width = 7, height = 5)