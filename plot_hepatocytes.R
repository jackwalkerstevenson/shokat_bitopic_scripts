#' ---
#'title: "plateplotr hepatocytes"
#'author: "Jack Stevenson"
#' ---
#'copied from microsomes script 2024-02-19
#'expects the following columns:
#'- 'compound'
#'- 'species'
#'- 'half_life' (mean value, minutes)
#'- 'half_life_SE' (standard error)
#'- 'CLint' (mean value, intrinsic clearance, µg/mL/million cells)
#'- 'CLint_SE' (standard error)

# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(scales) # for fancy plotting scales
library(ggprism)  # for pretty prism-like plots
library(assertthat) # for QC assertions
library(doseplotr) # you bet
# set global variables---------------------------------------------------------
input_filename <- "input/hepatocytes.xlsx"
plot_type <- "pdf" # file type for saved output plots
dir.create("output/", showWarnings = FALSE) # silently create output directory
all_data <- readxl::read_excel(input_filename, .name_repair = "universal") %>%
  # get greater/less information from e.g. "<0.2", then strip values to numbers
  mutate(half_life_greater_than = if_else(startsWith(half_life, ">"), TRUE, FALSE)) %>%
  mutate(CLint_less_than = if_else(startsWith(CLint, "<"), TRUE, FALSE)) %>%
  mutate(across(c(half_life, CLint), ~ str_replace_all(.x, "[><]", ""))) %>%
  mutate(across(c(half_life, CLint), as.numeric))
# check that data is from only one species
species <- distinct(all_data["species"])$species %>% str_to_lower()
assert_that(length(species) == 1,
            msg = "Cannot plot data from multiple species at once")
# reorder compound factors so they plot in input order
compounds <- distinct(all_data["compound"])$compound
all_data <- all_data %>% mutate(compound = fct_relevel(compound, compounds))
# global axis labels-----------------------------------------------------------
ylab = "compound"
# plot CLint-------------------------------------------------------------------
CLint_less <- all_data %>% filter(CLint_less_than) # select data to note less than
all_data %>%
  ggplot(aes(y = compound)) +
  scale_y_discrete(limits = rev) + # order top to bottom
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + # remove space from bottom of bar
  geom_bar(aes(x = CLint), stat = "identity", fill = "gray") +
  geom_errorbar(aes(x = CLint,
                    xmin = CLint - CLint_SE,
                    xmax = CLint + CLint_SE),
                width = .5) +
  # label appropriate subset of data with 'less than' symbol
  geom_text(data = CLint_less, aes(x = CLint * 2, label = "<")) +
  labs(
    title = str_glue("Clearance rate in {species} hepatocytes"),
    x = "Intrinsic clearance (µL/min/million cells)",
    y = ylab) +
  theme_prism() +
  theme(plot.background = element_blank()) # need for transparent background

ggsave(str_glue("output/CLint{get_timestamp()}.{plot_type}"),
       bg = "transparent", width = 7, height = 4)
# plot half life---------------------------------------------------------------
half_life_greater <- all_data %>% filter(half_life_greater_than) # select data to note greater than
all_data %>%
  ggplot(aes(y = compound)) +
  scale_y_discrete(limits = rev) + # order top to bottom
  scale_x_continuous(trans = "log10",
                     expand = expansion(mult = c(0, 0.05))) +
  geom_bar(aes(x = half_life), stat = "identity", fill = "gray") +
  geom_errorbar(aes(x = half_life,
                    xmin = half_life - half_life_SE,
                    xmax = half_life + half_life_SE),
                width = .5) +
  # label appropriate subset of data with 'less than' symbol
  geom_text(data = half_life_greater, aes(x = half_life * 1.1, label = ">")) +
  labs(
    title = str_glue("Half-life in {species} hepatocytes"),
    x = "Half life (minutes)",
    y = ylab) +
  theme_prism() +
  theme(plot.background = element_blank()) # need for transparent background

ggsave(str_glue("output/half_life_{get_timestamp()}.{plot_type}"),
       bg = "transparent", width = 7, height = 4)