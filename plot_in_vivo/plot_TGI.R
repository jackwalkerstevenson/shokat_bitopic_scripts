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
# write timestamped input files to output
doseplotr::file_copy_to_dir(input_path, output_directory)
doseplotr::file_copy_to_dir(dosing_path, output_directory)
doseplotr::file_copy_to_dir(PK_path, output_directory)
# set up input and output directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
all_TGI_data <- readxl::read_excel(
  str_glue("{input_directory}/{input_filename}"),
  # interpret existing string "NA" as NA
  na = "NA") |>
  dplyr::mutate(volume = as.numeric(volume),
         body_weight_percent = as.numeric(body_weight_percent) * 100,
         day = as.numeric(day),
         treatment = fct_inorder(as.factor(treatment)),
         hit_endpoint = tumor_length_mm > 20 | is.na(tumor_length_mm))

plot_TGI_data <- all_TGI_data |>
  # remove timepoints after too many subjects hit endpoint
  # if >1 measurement in a group timepoint hits endpoint, remove the whole timepoint
  dplyr::group_by(treatment, day) |> 
  dplyr::mutate(endpoint_count = sum(hit_endpoint)) |> 
  dplyr::ungroup() |> 
  dplyr::filter(endpoint_count < 2) |> 
  select(-endpoint_count)

dosing_data <- readxl::read_excel(
  str_glue("{input_directory}/{dosing_filename}"),
)

PK_data <- readxl::read_excel(
  str_glue("{input_directory}/{PK_filename}"),
) |> 
  dplyr::mutate(measurement = fct_inorder(as.factor(measurement)),
                treatment = fct_inorder(as.factor(treatment)))
# report raw TGI data
write_csv(all_TGI_data,
          fs::path(output_directory,
                   str_glue("raw_data_TGI_{get_timestamp()}.csv")))

# function to plot mean/SEM of a type of measurement over time----------------------------
plot_measurement <- function(measurement_data,
                             interval_data, 
                             measurement,
                             y_label,
                             plot_title){
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
    # dosing period rectangles
    # todo: truncate dosing period data to available measurement data
    # scale_x_continuous(limits = c(x_min, x_max),
    #                    expand = expansion(mult = c(0.05, 0.1))) +
    # geom_rect(data = interval_data,
    #           inherit.aes = FALSE,
    #           aes(xmin = dosing_start,
    #               xmax = dosing_end, #todo: pmin(dosing_end, x_max),
    #               ymin = -Inf,
    #               ymax = Inf,
    #               fill = "dosing period"),
    #           alpha = 0.1) +
    scale_color_manual(values = color_map_treatments) +
    scale_shape_manual(values = shape_map_treatments) +
    # scale_fill_manual(values = c("dosing period" = "black"), name = NULL) +
    theme_prism() +
    theme(legend.title = element_text(), # reinstate legend label
          plot.background = element_blank()) + # transparent background
    # guides(color = guide_legend(order = 1),
    #        shape = guide_legend(order = 1),
    #        fill = guide_legend(order = 2)) +
    labs(x = "days since start of dosing",
         y = y_label,
         title = plot_title,
         caption = animal_caption)
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
                 y_label = "tumor volume (mm^3)",
                 plot_title = "Tumor volume")
plot_measurement(measurement_data = plot_TGI_data,
                 interval_data = dosing_data, 
                 measurement = "body_weight_percent",
                 y_label = "body weight (%)",
                 plot_title = "Body weight")

# function to plot individual traces of measurement-------------------
plot_measurement_individual <- function(measurement_data,
                                        interval_data, 
                                        measurement,
                                        y_label,
                                        plot_title){
  x_min = min(measurement_data$day)
  x_max = max(measurement_data$day)
  plot <- measurement_data |> 
    ggplot(aes(x = day,
               y = .data[[measurement]],
               color = animal,
               shape = animal,
               group = animal)) +
    geom_point(size = 2) +
    geom_line() +
    # dosing period rectangles
    # geom_rect(data = interval_data,
    #           inherit.aes = FALSE,
    #           aes(xmin = dosing_start,
    #               xmax = dosing_end, #todo: pmin(dosing_end, x_max),
    #               ymin = -Inf,
    #               ymax = Inf,
    #               fill = "dosing period"),
    #           alpha = 0.1) +
    scale_color_manual(values = pals::cols25(8)) +
    scale_shape_manual(values = doseplotr::shape_scale_default()) +
    # dosing period rectangle stuff
    # scale_fill_manual(values = c("dosing period" = "black")) +
    # guides(color = "none", shape = "none",
    #        fill = guide_legend())+
    theme_prism() +
    theme(plot.background = element_blank(), # transparent background
          legend.position = "none") + # remove legend entirely
    labs(x = "days since start of dosing",
         y = y_label,
         title = plot_title,
         caption = animal_caption)
  return(plot)
}
# plot individual tumor vol and body weight-------------------------------------
for (trt in unique(all_TGI_data$treatment)){
  trt_data <- all_TGI_data |> 
    dplyr::filter(treatment == trt)
  trt_name <- janitor::make_clean_names(trt) # remove special chars from trt
  trt_data |> 
    plot_measurement_individual(
                 interval_data = dosing_data, 
                 measurement = "volume",
                 y_label = "tumor volume (mm³)",
                 plot_title = str_glue("Tumor volume, {trt}"))
  ggsave(
    str_glue(
      "{output_directory}/plot_TGI_individual_vol_{trt_name}_{get_timestamp()}.{plot_type}"
    ),
    width = 6, height = 4,
    bg = "transparent")
  
  trt_data |> 
    plot_measurement_individual(
      interval_data = dosing_data, 
      measurement = "body_weight_percent",
      y_label = "body weight (%)",
      plot_title = str_glue("Body weight, {trt}"))
  ggsave(
    str_glue(
      "{output_directory}/plot_TGI_individual_weight_{trt_name}_{get_timestamp()}.{plot_type}"
    ),
    width = 6, height = 4,
    bg = "transparent")
}
# plot Kaplan-Meier survival curve by group-------------------------------------
plot_survival_data <- all_TGI_data |> 
  dplyr::group_by(treatment, day) |> 
  dplyr::summarize(n_living = sum(hit_endpoint == FALSE)) # count living animals
plot_survival_data |> 
  ggplot(aes(x = day, y = n_living, color = treatment)) +
  # horizontal step first means line goes down on day when mouse is missing
  geom_step(direction = "hv",
            linewidth = 1,
            position = position_jitter(width = 1, height = 0)) +
  # dosing interval rectangles
  # geom_rect(data = dosing_data,
  #           inherit.aes = FALSE,
  #           aes(xmin = dosing_start,
  #               xmax = dosing_end, #todo: pmin(dosing_end, x_max),
  #               ymin = -Inf,
  #               ymax = Inf,
  #               fill = "dosing period"),
  #           alpha = 0.1) +
  # start x axis at start of dosing
  scale_x_continuous(limits = c(0, NA)) +
  # start y axis at 0 mice
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = color_map_treatments) +
  # scale_fill_manual(values = c("dosing period" = "black"), name = NULL) +
  theme_prism() +
  theme(legend.title = element_text(), # reinstate legend label
        plot.background = element_blank()) + # for transparent background
  # guides(color = guide_legend(order = 1),
  #        fill = guide_legend(order = 2)) +
  labs(x = "days since start of dosing",
       y = "surviving mice",
       title = "Survival",
       caption = animal_caption)
ggsave(
  str_glue(
    "{output_directory}/plot_TGI_survival_{get_timestamp()}.{plot_type}"
  ),
  width = 7, height = 4,
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
    "{output_directory}/plot_TGI_PK_{get_timestamp()}.{plot_type}",
    bg = "transparent"
  ),
  width = 9, height = 6,
  bg = "transparent")