#' ---
#'title: "plateplotr selectscreen"
#'author: "Jack Stevenson"
#'date: "2022-09-09"
#' ---
#'edited version of plateplotr for dealing with SelectScreen data
#'copied from main plot.R 2022-09-09

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(patchwork) # for plot organization
library(assertthat) # for QC assertions

# import and tidy data---------------------------------
# note the order compounds are imported is the order they will be plotted
input_filename <- "ZLYTE_compiled_results_complete.csv"
source("compounds.R")
source("targets.R")
plot_type <- "pdf"
dir.create("output/", showWarnings = FALSE)
plate_data <- read_csv(input_filename) %>%
  rename(compound = Compound) %>%
  rename(target = Kinase) %>%
  # filter for desired compounds
  filter(compound %in% compounds) %>%
  filter(target %in% targets) %>%
  mutate(compound = fct_relevel(compound, compounds)) %>% # order compounds by list
  # tidy by pivoting to one measurement per row. magic column number, watch out!
  pivot_longer(cols = 11:12, names_to = NULL, values_to = "pct_inhibition") %>%
  # wrangle: convert conc and convert inhibition to activity
  mutate(log.conc = log10(Compound_Conc_nM/1e9)) %>%  # convert conc nM to log(M)
  mutate(activity = 100 - pct_inhibition)
# fit models to output EC values------------------------------------------------
# seems like you should be able to just pipe group_by into drm(), but nope, so doing this instead
# helper function for getting EC for one compound, cell line, and EC threshold
source("get_EC.R")
EC_summary <- plate_data %>%
  group_by(compound, target) %>%
  summarize(
    EC50_nM = 10^get_EC(plate_data, compound, target, 50) * 1e9, # convert M to nM
    # for negative-response data like this, the EC75 is the drop to 25%
    EC75_nM = 10^get_EC(plate_data, compound, target, 25) * 1e9
  )
write_csv(EC_summary, "output/EC_summary_selectscreen.csv")
# generate global parameters for all plots------------------------------------------
pt_size = 2 # size for all geom_point
all_compounds <- distinct(plate_data["compound"])$compound
all_targets <- distinct(plate_data["target"])$target
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plate_data$log.conc))
# x_min <- -10
x_max <- ceiling(max(plate_data$log.conc))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all conc plots
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# set factors so targets get plotted and colored in input order
#compound_factors <- distinct(plate_data, compound)$compound
target_factors <- distinct(plate_data, target)$target
plate_data <- plate_data %>% 
  mutate(target = fct_relevel(target, target_factors))# %>%
#  mutate(compound = fct_relevel(compound, compound_factors))
# helper function for saving plots----------------------------------------------
scale_facet <- 4 # plot width per col/height per row
#todo: allow overriding width and height if provided
save_plot <- function(plot, nrow = 1, ncol = 1, ...){
  ggsave(plot, bg = "transparent",
         width = ncol*scale_facet + 2,
         height = nrow*scale_facet, ...)}
# helper function for summarizing replicate data for plotting------------------
plate_summarize <- function(x){
  summarize(x,
            # standard error for error bars = standard deviation / square root of n
            sem = sd(activity, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_read = mean(activity),
            w = 0.06 * n() # necessary for consistent error bar widths across plots
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
         y = "kinase activity (%)")
}
# helper function to plot one compound----------------------------------------
plot_compound <- function(cpd){
  plate_summary <- plate_data %>%
    filter(compound == cpd) %>% # get data from one compound to work with
    group_by(target, log.conc) %>%  # get set of replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = target)) +
      geom_point(size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w)) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1)} %>%
    plot_global() +
    theme(aspect.ratio = 1) +
    scale_color_viridis(discrete = TRUE, begin = 1, end = 0) +
    # scale_color_manual(values = c("black","darkred")) +
    labs(title = cpd)
}
# plot data for each compound separately----------------------------------------
for (cpd in all_compounds){
  plot_compound(cpd)
  # save plot with manually optimized aspect ratio
  save_plot(str_glue("output/{cpd}.{plot_type}"))
}  
# plot data for all compounds in facets----------------------------------
compound_plots = list()
for (cpd in all_compounds){
  compound_plots <- append(compound_plots, list(plot_compound(cpd)))
}

plot_mar <- 15 # margin between wrapped plots, in points
cols = ceiling(sqrt(length(compounds)))
rows = ceiling(length(compounds)/cols)
wrap_plots(compound_plots, guides = "collect", ncol = cols, nrow = rows) &
  theme(plot.margin = unit(c(plot_mar,plot_mar,plot_mar,plot_mar), "pt"),
        plot.background = element_blank(),
        legend.text= element_text(face = "bold", size = 16))
save_plot(str_glue("output/compound_facets.{plot_type}"), ncol = cols, nrow = rows)
# set color parameters for overlaid plots--------------------------------------
alpha_val <- 1
color_scale <- "viridis"
viridis_start <- .95
viridis_end <- 0
grey_start <- .7
grey_end <- 0
# plot data for each target separately------------------------------------------
for (k in all_targets){
  plate_summary <- plate_data %>%
    filter(target == k) %>%
    group_by(compound, log.conc) %>% # group into replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = compound)) +
      geom_point(aes(shape = compound), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
      # use drm method from drc package to fit dose response curve
      geom_line(#aes(linetype = compound),
        stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
    plot_global() +
    #scale_color_grey(start = grey_start, end = grey_end) +
    scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_start, end = viridis_end) +
    labs(title = k)
  save_plot(str_glue("output/{k}.{plot_type}"))
}
# plot data for all targets at once-----------------------------------------
plate_summary <- plate_data %>%
  group_by(target, compound, log.conc) %>% # group into replicates for each condition
  plate_summarize()
{ggplot(plate_summary,aes(x = log.conc, y = mean_read, color = compound)) +
    geom_point(size = pt_size) +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = target), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
  plot_global() +
  scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_start, end = viridis_end) +
  labs(title = "All data")
save_plot(str_glue("output/all_data.{plot_type}"))