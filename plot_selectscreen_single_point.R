#' ---
#'title: "plateplotr selectscreen single point"
#'author: "Jack Stevenson"
#'date: "2023-03-06"
#' ---
#'edited version of plateplotr for dealing with SelectScreen data
#'started 2023-03-06
#'# load required libraries-----------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(assertthat) # for QC assertions
# set global variables----------------------------------------------------------
input_filename <- "ZLYTE_compiled_results_single_point.csv"
plot_type <- "pdf" # file type for saved output plots
dir.create("output/", showWarnings = FALSE) # silently create output directory
source("parameters/compounds.R")
source("parameters/targets.R")
source("import_selectscreen.R")
all_data <- import_selectscreen(input_filename, compounds)
# helper summary function-------------------------------------------------------
inhibition_summarize <- function(x){
  summarize(x,
            compound = compound,
            # standard error for error bars = standard deviation / square root of n
            sem = sd(pct_inhibition, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_pct_inhibition = mean(pct_inhibition),
            w = 0.06 * n()
  )
}
# pilot plot--------------------------------------------------------------------
pt_size = 3 # size for all geom_point
alpha_val = 0.6
all_data %>%
  group_by(compound, target) %>%
  inhibition_summarize() %>%
ggplot(aes(x = target, color = compound)) +
  geom_point(aes(y = mean_pct_inhibition), size = pt_size, alpha = alpha_val) +
  # geom_line(aes(y = mean_pct_inhibition)) +
  geom_errorbar(aes(ymax = mean_pct_inhibition+sem, ymin = mean_pct_inhibition-sem, width = w), alpha = alpha_val) +
  scale_color_viridis(discrete = TRUE, begin = .9, end = 0) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_prism() +
  # rotated, right-justified x labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Single-point SelectScreen inhibition",
       caption = "Note: concentrations not equal between kinases",
       x = "target kinase",
       y = "percent inhibition")
  