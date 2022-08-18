#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2022-08-09"
#' ---
#'plateplotr generates dose-response plots from plate reader data.
#' 
#'Input: data formatted as per the plater package. Use the JWS Excel plater template to make this easy.
#' 
#'Variables expected in the import file (JWS Excel template is already set up for this):
#' 
#'- 'compound': name of compound used (including DMSO wells in dilution series)
#'- 'conc': concentration of compound used
#'- 'cell_line': cell line treated
#'- 'readout': raw plate reader data (this is only used to calculate read_norm)
#'- 'read_norm': normalized plate reader data (calculated in the spreadsheet)
#'- 'replicate': replicate of curve (not actually used)
#'
#'How to use plateplotr:
#'
#'1. Make one copy of the JWS Excel plater template for each plate of data you want to import
#'2. Fill out each copy with the compound + concs used and paste in raw plate reader data
#'3. Save input sheets as csv and put them in a folder within your working directory
#'4. Set `plate_directory` to the folder of .csv files (default is "Input CSVs/")
#'  -Maintain format of default, "Folder" and "~/Folder/" will both fail.
#'5. Run plot.R (this file)
#'6. Output plots will be created in a "Plots Output" folder in the working directory

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(patchwork) # for plot organization

# specify names of input files and import data---------------------------------
# note the order compounds are imported is the order they will be plotted
plate_directory <- "Input CSVs/"
plot_type <- "pdf"

plate_filenames <- c(list.files(plate_directory, pattern = "*.csv")) #gathers all .csv in directory
plate_paths <- paste0(plate_directory, plate_filenames)
plate_names <- seq(1,length(plate_filenames))  # create plate IDs
plate_data <- read_plates(plate_paths, plate_names) %>% # import with plater
  filter(compound != "N/A") %>% # drop empty wells
  # drop 0 values before plotting and curve fitting
  # note this is only OK because normalization happens before import
  filter(conc != 0) %>%
  mutate(log.conc = log10(conc/1e6))  # convert conc µM to M and log transform
# generate global parameters for all plots------------------------------------------
all_compounds <- distinct(plate_data["compound"])$compound
all_lines <- distinct(plate_data["cell_line"])$cell_line
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plate_data$log.conc))
x_max <- ceiling(max(plate_data$log.conc))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all compounds
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# set factors so cell lines get plotted and colored in input order
cell_line_factors <- distinct(plate_data, cell_line)$cell_line
compound_factors <- distinct(plate_data, compound)$compound
plate_data <- plate_data %>% 
  mutate(cell_line = fct_relevel(cell_line, cell_line_factors)) %>%
  mutate(compound = fct_relevel(compound, compound_factors))
# helper function for saving plots----------------------------------------------
scale_facet <- 4 # plot width per col/height per row
#todo: allow overriding width and height if provided
save_plot <- function(plot, nrow = 1, ncol = 1, ...){
  ggsave(plot, bg = "transparent",
         width = ncol*scale_facet + 1,
         height = nrow*scale_facet, ...)}
# helper function for summarizing replicate data for plotting------------------
plate_summarize <- function(x){
  summarize(x,
            # standard error for error bars = standard deviation / square root of n
            sem = sd(read_norm, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_read = mean(read_norm),
            w = 0.1 * n() # necessary for consistent error bar widths across plots
  )
}
# helper function to add ggplot objects common to all plots--------------------
plot_global <- function(plot){
  plot +
    scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                       minor_breaks = minor_x) + # manual minor ticks
    scale_y_continuous(guide = "prism_offset",  # end at last tick
                       breaks = c(0,25,50,75,100)) + # manual y axis ticks
    coord_cartesian(xlim = x_limits, # set x axis zoom from global values
                    ylim = c(0,NA)) + # set y axis zoom locally
    theme_prism() + # make it look fancy like prism
    theme(plot.background = element_blank()) + # need for transparent background
    labs(x = "log [compound] (M)",
         y = "relative cell viability (%)")
}
# helper function to plot one compound----------------------------------------
plot_compound <- function(cpd){
  plate_summary <- plate_data %>%
    filter(compound == cpd) %>% # get data from one compound to work with
    group_by(cell_line, log.conc) %>%  # get set of replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = cell_line)) +
      geom_point() +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w)) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, size = 1)} %>%
    plot_global() +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = c("black","darkred")) +
    labs(title = cpd)
}
# plot data for each compound separately----------------------------------------
for (cpd in all_compounds){
  plot_compound(cpd)
  # save plot with manually optimized aspect ratio
  save_plot(str_glue("plots output/{cpd}.{plot_type}"))
}  
# plot data for all compounds in facets----------------------------------
compound_plots = list()
for (cpd in all_compounds){
  compound_plots <- append(compound_plots, list(plot_compound(cpd)))
}

plot_mar <- 15 # margin between wrapped plots, in points
cols = 4
rows = 2
wrap_plots(compound_plots, guides = "collect", ncol = cols, nrow = rows) &
  theme(plot.margin = unit(c(plot_mar,plot_mar,plot_mar,plot_mar), "pt"),
        plot.background = element_blank(),
        legend.text= element_text(face = "bold", size = 16))
save_plot(str_glue("plots output/compound_facets.{plot_type}"), ncol = cols, nrow = rows)
# plot data for each cell line separately-------------------------------------------------
alpha_val <- 1
viridis_start <- .8
viridis_end <- 0
for (c_line in all_lines){
  plate_summary <- plate_data %>%
    filter(cell_line == c_line) %>%
    group_by(compound, log.conc) %>% # group into replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = compound)) +
      geom_point() +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, size = 1, alpha = alpha_val)} %>%
    plot_global() +
    scale_color_viridis(option = "turbo", discrete = TRUE, begin = viridis_start, end = viridis_end) +
    labs(title = c_line)
  save_plot(str_glue("plots output/{c_line}.{plot_type}"))
}
# plot data for all cell lines at once-----------------------------------------
alpha_val <- 1
viridis_start <- .8
viridis_end <- 0
plate_summary <- plate_data %>%
  group_by(cell_line, compound, log.conc) %>% # group into replicates for each condition
  plate_summarize()
{ggplot(plate_summary,aes(x = log.conc, y = mean_read, color = compound)) +
    geom_point() +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = cell_line), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, size = 1, alpha = alpha_val)} %>%
  plot_global() +
  scale_color_viridis(option = "turbo", discrete = TRUE, begin = viridis_start, end = viridis_end) +
  labs(title = "All data")
save_plot(str_glue("Plots Output/all_data.{plot_type}"))