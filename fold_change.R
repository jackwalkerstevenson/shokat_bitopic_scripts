# plot fold changes between targets relative to wt for each treatment
# load libraries----------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# import precalculated IC50 table-----------------------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_fold_change.R")
data <- read_csv(input_filename)
# process data---------------------------------------------------------------
if(exists("treatments")){
  data <- filter_validate_reorder(data, "treatment", treatments)
}
if(exists("targets")){
  data <- filter_validate_reorder(data, "target", targets)
}
# extract wt IC50 and calculate fold change of other targets
wt_IC50s <- data |>
  filter(target == wt_target_name) |> 
  group_by(treatment) |> 
  summarize(wt_IC50_nM = IC50_nM)
data <- data |> 
  left_join(wt_IC50s, by = "treatment") |> 
  mutate(fold_vs_wt_IC50 = IC50_nM/wt_IC50_nM)
# write report of fold changes-----------------------------------------------
report_data <- data |>
  mutate(fold_vs_wt_IC50 = fold_vs_wt_IC50 |> signif(3),
         IC50_nM = IC50_nM |> signif(3)) |> 
  dplyr::select(treatment, target, IC50_nM, fold_vs_wt_IC50)
# tidy report-----------------------------------------------
report_data |>
  write_csv(str_glue("output/fold_change_report_tidy_{get_timestamp()}.csv"))
# untidy report-----------------------------------------------
report_data |> 
  tidyr::pivot_wider(names_from = treatment, values_from = c(IC50_nM, fold_vs_wt_IC50)) |> 
  write_csv(str_glue("output/fold_change_report_untidy_{get_timestamp()}.csv"))

# helper functions for renaming-----------------------------------------------
get_target_labels <- function(){
  if(manually_relabel_targets) display_names_targets else{
    # necessary to preserve correct ordering even after changing plotting order
    named_targets <- targets
    names(named_targets) <- targets
    return(named_targets)
  }
}
legend_labels = get_if(display_names_treatments,
                       manually_relabel_treatments,
                       otherwise = ggplot2::waiver())
# bar plot of fold changes------------------------------------------------------
vr <- viridis_range(length(treatments))
viridis_begin <- vr[[1]]
viridis_end <- vr[[2]]
viridis_option <- vr[[3]]
x_min <- floor(min(log10(data$fold_vs_wt_IC50)))
x_max <- ceiling(max(log10(data$fold_vs_wt_IC50)))
breaks_x <- 10^rep(x_min : x_max)
minor_x <- minor_breaks(x_min, x_max, log_units = TRUE)
p <- data |> 
  # don't plot wt or control
  filter(!target %in% c(wt_target_name, control_target_name)) |>
  ggplot(aes(y = target, x = fold_vs_wt_IC50,
             fill = treatment,
             label = signif(fold_vs_wt_IC50, digits = 2))) +
  geom_bar(stat = "identity",
           position = position_dodge2(reverse = TRUE, padding = 0)) +
  geom_text(position = position_dodge2(width = .9, reverse = TRUE),
            hjust = -0.1,
            parse = FALSE,
            size = 5) +
  scale_x_continuous(trans = "log10",
                     expand = expansion(mult = .1),
                     guide = "prism_offset_minor", # end at last tick
                     breaks = breaks_x,
                     labels = label_comma(accuracy = 1, big.mark = ""),
                     minor_breaks = minor_x,) +
  scale_y_discrete(limits = rev, labels = get_target_labels()) +
  {if(manually_recolor_treatments){
    ggplot2::scale_fill_manual(values = color_map_treatments,
                               labels = legend_labels)
  } else{
    scale_fill_viridis(option = viridis_option,
                       discrete = TRUE,
                       begin = viridis_begin, end = viridis_end,
                       labels = legend_labels)
  }} +
  theme_prism() +
  theme(plot.background = element_blank(), # need for transparent background
        legend.title = element_text(face = "plain"),
        legend.title.align = 0) +
  labs(x = "fold change in IC50 vs wt",
       y = "cell line (K562 pUltra BCR-ABL1)")
save_plot(p, str_glue("output/fold_change_bar_{get_timestamp()}.{plot_type}"),
          width = 14, height = .7*length(targets) + 0.75)
# strip plot of raw IC50s-----------------------------------------------
x_min <- floor(min(log10(data$IC50_nM)))
x_max <- ceiling(max(log10(data$IC50_nM)))
breaks_x <- 10^rep(x_min : x_max)
minor_x <- rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9))
p <- data |> 
  ggplot(aes(y = target, x = IC50_nM,
             color = treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_x_continuous(trans = c("log10", "reverse"),
                     guide = "prism_minor", # end at last tick
                     breaks = breaks_x,
                     labels = label_comma(accuracy = 1, big.mark = ""),
                     minor_breaks = minor_x,
                     expand = expansion(mult = .1)) +
  scale_y_discrete(limits = rev, labels = get_target_labels()) +
  {if(manually_recolor_treatments){
    ggplot2::scale_color_manual(values = color_map_treatments,
                                labels = legend_labels)
  } else{
    scale_fill_viridis(option = viridis_option,
                       discrete = TRUE,
                       begin = viridis_begin, end = viridis_end,
                       labels = legend_labels)
  }} +
  theme_prism() +
  theme(panel.grid = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor = element_line(color = "black",
                                        linewidth = 0.1,
                                        linetype = "dotted")) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "IC50 (nM)",
       y = "cell line (K562 pUltra BCR-ABL1)")
save_plot(p, str_glue("output/IC50_dot_{get_timestamp()}.{plot_type}"),
          width = 12, height = .65*length(targets) + 0.75)
