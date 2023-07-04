# plot fold changes between targets relative to wt for each treatment
# load libraries----------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# import precalculated IC50 table-----------------------------------------------
input_filename <- "input/EC_summary_2023-07-03T204947.csv"
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
  summarize(wt_IC50_nM = IC50_nM)
data <- data |> 
  left_join(wt_IC50s, by = "treatment") |> 
  mutate(fold_vs_wt_IC50 = IC50_nM/wt_IC50_nM)
# plot fold changes-----------------------------------------------
plot_type = "pdf"
color_scale <- "viridis"
vr <- viridis_range(length(treatments))
viridis_begin <- vr[1]
viridis_end <- vr[2]
p <- data |> 
  filter(target != "K562 pUltra BCR-ABL1 wt") |> 
  ggplot(aes(x = target, y = fold_vs_wt_IC50,
             fill = treatment,
             label = round(fold_vs_wt_IC50, digits = 1))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 1), vjust = -0.5) +
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis(option = color_scale,
                      discrete = TRUE,
                      begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(y = "fold change in IC50 vs wt")
ggsave(str_glue("output/fold_change_{get_timestamp()}.{plot_type}"),
       bg = "transparent",
       width = 14, height = 10)
  # save_plot(p, str_glue("output/fold_change_{get_timestamp()}.{plot_type}"),
            # width = 12, height = 8))
