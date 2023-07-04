# plot fold changes between targets relative to wt for each treatment
# load libraries----------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# import precalculated IC50 table-----------------------------------------------
input_filename <- "input/Ivan IC50 data vs 06-21.csv"
source("parameters/treatments.R") # import list of treatments to include in plots
source("parameters/targets.R") # import list of targets to include in plots
data <- read_csv(input_filename) |> 
  filter(treatment %in% treatments) |> # take only specified treatments
  filter(target %in% targets) |> # take only specified targets
  mutate(target = fct_relevel(target, targets)) |> 
  mutate(treatment = fct_relevel(treatment, treatments))
wt_IC50s <- data |>
  filter(target == "K562 pUltra BCR-ABL1 wt") |>
  group_by(treatment) |>
  summarize(wt_IC50_nM_0621 = IC50_nM_0621,
            wt_IC50_nM_0703 = IC50_nM_0703)
data <- data |> 
  left_join(wt_IC50s, by = "treatment") |>
  mutate(IC50_0703_over_0621 = IC50_nM_0703 / IC50_nM_0621,
         IC50_vs_wt_0621 = IC50_nM_0621 / wt_IC50_nM_0621,
         IC50_vs_wt_0703 = IC50_nM_0703 / wt_IC50_nM_0703,
         IC50_vs_wt_0703_over_0621 =
           IC50_vs_wt_0703 / IC50_vs_wt_0621)
# set up variables for plotting-----------------------------------------------
plot_type = "pdf"
color_scale <- "viridis"
vr <- viridis_range(length(treatments))
viridis_begin <- vr[1]
viridis_end <- vr[2]
# helper function to plot a column, not working yet-----------------------------
plot_col <- function(df, col_name){
  p <- df |> 
    ggplot(aes(x = target, y = col_name,
               fill = treatment,
               label = round(col_name, digits = 1))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(position = position_dodge(width = 1), vjust = -0.5) +
    scale_y_continuous(trans = "log10") +
    scale_fill_viridis(option = color_scale,
                       discrete = TRUE,
                       begin = viridis_begin, end = viridis_end) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.background = element_blank()) + # need for transparent background
    labs(y = "fold change in IC50 07-03 vs 06-21")
  ggsave(str_glue("output/fold_change_{col_name}_{get_timestamp()}.{plot_type}"),
         p, bg = "transparent",
         width = 10, height = 10)
  # save_plot(p, str_glue("output/fold_change_date_comparison_{get_timestamp()}.{plot_type}"),
  # width = 12, height = 8)) 
}
# plot IC50 vs wt 06-21---------------------------------------
p <- data |> 
  ggplot(aes(x = target, y = IC50_vs_wt_0621,
             fill = treatment,
             label = round(IC50_vs_wt_0621, digits = 1))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 1), vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis(option = color_scale,
                     discrete = TRUE,
                     begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(y = "fold change in IC50 vs wt 06-21")
ggsave(str_glue("output/IC50_vs_wt_0621_{get_timestamp()}.{plot_type}"),
       bg = "transparent",
       width = 10, height = 10)
# save_plot(p, str_glue("output/fold_change_date_comparison_{get_timestamp()}.{plot_type}"),
# width = 12, height = 8))
# plot IC50s 07-03-------------------------------------------------------------
p <- data |> 
  ggplot(aes(x = target, y = IC50_vs_wt_0703,
             fill = treatment,
             label = round(IC50_vs_wt_0703, digits = 1))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 1), vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis(option = color_scale,
                     discrete = TRUE,
                     begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(y = "fold change in IC50 vs wt 07-03")
ggsave(str_glue("output/IC50_vs_wt_0703_{get_timestamp()}.{plot_type}"),
       bg = "transparent",
       width = 10, height = 10)
# save_plot(p, str_glue("output/fold_change_date_comparison_{get_timestamp()}.{plot_type}"),
# width = 12, height = 8))
# plot fold change in IC50s between dates---------------------------------------
p <- data |> 
  ggplot(aes(x = target, y = IC50_0703_over_0621,
             fill = treatment,
             label = round(IC50_0703_over_0621, digits = 1))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 1), vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis(option = color_scale,
                      discrete = TRUE,
                      begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(y = "fold change in IC50 07-03 vs 06-21")
ggsave(str_glue("output/fold_change_date_comparison_{get_timestamp()}.{plot_type}"),
       bg = "transparent",
       width = 10, height = 10)
  # save_plot(p, str_glue("output/fold_change_date_comparison_{get_timestamp()}.{plot_type}"),
            # width = 12, height = 8))
# plot fold change in IC50s vs wt between dates---------------------------------
p <- data |> 
  ggplot(aes(x = target, y = IC50_vs_wt_0703_over_0621,
             fill = treatment,
             label = round(IC50_vs_wt_0703_over_0621, digits = 1))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 1), vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis(option = color_scale,
                     discrete = TRUE,
                     begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(y = "fold change in IC50 vs wt 07-03 vs 06-21")
ggsave(str_glue("output/fold_change_vs_wt_date_comparison_{get_timestamp()}.{plot_type}"),
       bg = "transparent",
       width = 10, height = 10)
# save_plot(p, str_glue("output/fold_change_vs_wt_dates_{get_timestamp()}.{plot_type}"),
# width = 12, height = 8))
