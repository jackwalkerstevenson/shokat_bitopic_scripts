# plot tumor growth inhibition data
# Jack Stevenson started 2024-10-09
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_TGI.R"
source(params_path)
input_path <- str_glue("{input_directory}{input_filename}")
dosing_path <- str_glue("{input_directory}{dosing_filename}")
PK_path <- str_glue("{input_directory}{PK_filename}")
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_in_vivo/plot_TGI.R", output_directory)
# write timestamped input file to output
doseplotr::file_copy_to_dir(input_path, output_directory)
# set up input and output directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
raw_TGI_data <- readxl::read_excel(
  str_glue("{input_directory}/{input_filename}"),
  na = "NA" # interpret existing string "NA" as NA
)

plot_TGI_data <- raw_TGI_data |>
  mutate(volume = as.numeric(volume),
         body_weight_percent = as.numeric(body_weight_percent) * 100,
         day = as.numeric(day),
         treatment = fct_inorder(as.factor(treatment)))
dosing_data <- readxl::read_excel(
  str_glue("{input_directory}/{dosing_filename}"),
)

PK_data <- readxl::read_excel(
  str_glue("{input_directory}/{PK_filename}"),
) |> 
  dplyr::mutate(measurement = fct_inorder(as.factor(measurement)),
                treatment = fct_inorder(as.factor(treatment)))
# report raw TGI data
write_csv(raw_TGI_data,
          fs::path(output_directory,
                   str_glue("raw_data_TGI_{get_timestamp()}.csv")))
# report raw dosing period data
write_csv(dosing_data,
          fs::path(output_directory,
                   str_glue("raw_data_dosing_{get_timestamp()}.csv")))
# report raw PK data
write_csv(PK_data,
          fs::path(output_directory,
                   str_glue("raw_data_PK_{get_timestamp()}.csv")))

# function to plot individual measurements over time----------------------------
plot_measurement <- function(measurement_data,
                             interval_data, 
                             measurement,
                             y_label){
  x_min = min(measurement_data$day)
  x_max = max(measurement_data$day)
  measurement_data |> 
    ggplot(aes(x = day,
               y = .data[[measurement]],
               color = treatment,
               shape = treatment,
               group = treatment)) +
    stat_summary(
      fun.data = "mean_se",
      geom = "errorbar",
      width = 1
    ) +
    stat_summary(
      fun = "mean",
      geom = "point",
      size = 3
    ) +
    stat_summary(
      fun = "mean",
      geom = "line"
    ) +
    # scale_x_continuous(limits = c(x_min, x_max),
    #                    expand = expansion(mult = c(0.05, 0.1))) +
    scale_color_manual(values = color_map_treatments) +
    scale_shape_manual(values = shape_map_treatments) +
    geom_rect(data = interval_data,
              inherit.aes = FALSE,
              aes(xmin = dosing_start,
                  xmax = dosing_end, #todo: pmin(dosing_end, x_max),
                  ymin = -Inf,
                  ymax = Inf),
              fill = "black",
              alpha = 0.2) +
    theme_prism() +
    theme(legend.title = element_text(), # reinstate legend label
          panel.background = element_blank()) + # transparent background
    labs(x = "days since start of dosing",
         y = y_label)
  ggsave(
    str_glue(
      "{output_directory}/plot_TGI_{measurement}_{get_timestamp()}.{plot_type}"
    ),
    width = 8, height = 4,
    bg = "transparent")
}
# plot tumor volume and body weight---------------------------------------------
plot_measurement(measurement_data = plot_TGI_data,
                 interval_data = dosing_data, 
                 measurement = "volume",
                 y_label = "tumor volume (mm^3)")
plot_measurement(measurement_data = plot_TGI_data,
                 interval_data = dosing_data, 
                 measurement = "body_weight_percent",
                 y_label = "body weight (%)")

# plot PK data-----------------------------------------------------------------
PK_data |> 
  group_by(treatment, measurement) |> 
  dplyr::summarize(
    mean_conc_nM = mean(conc_nM),
    # standard error for error bars = standard deviation / square root of n
    sem = sd(conc_nM, na.rm = TRUE)/sqrt(n()),
  ) |> 
  ggplot(aes(x = treatment,
             y = mean_conc_nM,
             fill = treatment,
             group = measurement)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  geom_errorbar(aes(ymin = mean_conc_nM - sem,
                    ymax = mean_conc_nM + sem),
                position = "dodge",
                linewidth = 0.3) +
  geom_label(aes(label = measurement), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            show.legend = FALSE) + 
  scale_fill_manual(values = color_map_treatments) +
  theme_prism() +
  labs(x = "",
       y = "plasma concentration (nM)",
       title = "Peak and trough compound concentrations",
       caption = "after 14 d QD weekday dosing")
ggsave(
  str_glue(
    "{output_directory}/plot_TGI_PK_{get_timestamp()}.{plot_type}"
  ),
  width = 9, height = 6,
  bg = "transparent")
# plot PK data-----------------------------------------------------------------
PK_data |> 
  group_by(treatment, measurement) |> 
  dplyr::summarize(
    mean_conc_nM = mean(conc_nM),
    # standard error for error bars = standard deviation / square root of n
    sem = sd(conc_nM, na.rm = TRUE)/sqrt(n()),
  ) |> 
  ggplot(aes(x = treatment,
             y = mean_conc_nM,
             fill = treatment,
             group = measurement)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  geom_errorbar(aes(ymin = mean_conc_nM - sem,
                    ymax = mean_conc_nM + sem),
                position = "dodge",
                linewidth = 0.3) +
  geom_label(aes(label = measurement), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            show.legend = FALSE) + 
  scale_fill_manual(values = color_map_treatments) +
  theme_prism() +
  labs(x = "",
       y = "plasma concentration (nM)",
       title = "Peak and trough compound concentrations",
       caption = "after 14 d QD weekday dosing")
ggsave(
  str_glue(
    "{output_directory}/plot_TGI_PK_{get_timestamp()}.{plot_type}"
  ),
  width = 9, height = 6,
  bg = "transparent")