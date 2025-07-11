# plot tumor growth inhibition data
# Jack Stevenson started 2024-10-09
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(survival) # for model fitting
library(survminer)
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_TGI.R"
source(params_path)
input_path <- str_glue("{input_directory}{input_filename}")
dosing_path <- str_glue("{input_directory}{dosing_filename}")
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_in_vivo/plot_TGI.R", output_directory)
# write timestamped input files to output
doseplotr::file_copy_to_dir(input_path, output_directory)
doseplotr::file_copy_to_dir(dosing_path, output_directory)
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
         animal = fct_inorder(as.factor(animal)),
         hit_endpoint = tumor_length_mm > 20 | is.na(tumor_length_mm))

last_days <- all_TGI_data |> 
  dplyr::group_by(animal, treatment) |> 
  dplyr::summarize(last_day=max(day))

endpoint_days <- all_TGI_data |> 
  dplyr::group_by(animal, treatment) |> 
  dplyr::filter(hit_endpoint == TRUE) |>
  dplyr::summarize(endpoint_day = min(day))

survival_data <- last_days |> 
  dplyr::left_join(endpoint_days, by=c("animal", "treatment")) |> 
  dplyr::mutate(reached_endpoint=!is.na(endpoint_day),
                endpoint_day=dplyr::coalesce(endpoint_day, last_day)) |> 
  dplyr::select(-last_day)

# trim data for plotting: remove group after too many subjects hit endpoint
# if >0 measurement in a group timepoint hits endpoint, remove the whole timepoint
plot_TGI_data <- all_TGI_data |>
  dplyr::group_by(treatment, day) |> 
  dplyr::mutate(endpoint_count = sum(hit_endpoint)) |> 
  dplyr::ungroup() |> 
  dplyr::filter(endpoint_count < 1) |> # choice of how many endpoints can appear
  dplyr::select(-endpoint_count)

dosing_data <- readxl::read_excel(
  str_glue("{input_directory}/{dosing_filename}"),
)

# report raw TGI data
write_csv(all_TGI_data,
          fs::path(output_directory,
                   str_glue("raw_data_TGI_{get_timestamp()}.csv")))
# report survival summary statistics
survival_summary <- survival_data |> 
  dplyr::group_by(treatment) |> 
  summarize(
    n = n(),
    total_reached_endpoint = sum(reached_endpoint),
    event_rate = mean(reached_endpoint),
    min_time_to_endpoint = ifelse(any(reached_endpoint), min(endpoint_day[reached_endpoint], na.rm=TRUE), NA),
    max_time_to_endpoint = ifelse(all(reached_endpoint), max(endpoint_day[reached_endpoint], na.rm=TRUE), NA),
  )
write_csv(survival_summary,
          fs::path(output_directory,
                   str_glue("group_summary_TGI_{get_timestamp()}.csv")))
# experiment with survival modeling-----------------------------------------------
# mostly claude
surv_obj <- Surv(time = survival_data$endpoint_day,
                 event = survival_data$reached_endpoint)
km_fit <- survfit(surv_obj ~ treatment, data = survival_data)
km_plot <- ggsurvplot(
  km_fit,
  data = survival_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Survival probability",
  title = "Time to Endpoint Analysis",
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)
# overall log-rank test
overall_log_rank <- survdiff(surv_obj ~ treatment, data = survival_data)

# pairwise log-rank tests vs vehicle-----------------------------------------------
all_treatments <- levels(all_TGI_data$treatment)
non_control_treatments <- all_treatments[all_treatments != control_treatment]
# get pval from a log-rank test between a treatment and control as a tibble
get_log_rank_pval_vs_ctrl <- function(data, trt){
    # Subset data to just control and specified treatment
    control_and_trt_data <- data |> 
      dplyr::filter(treatment %in% c(control_treatment, trt))
    # run log-rank test
    test_result <- survdiff(
      Surv(time = endpoint_day,
           event = reached_endpoint) ~ treatment,
      data = control_and_trt_data,
      rho = 0) # rho 0 -> log-rank test
    # Extract p-value
    tibble(
      treatment = trt,
      p.value = 1 - pchisq(test_result$chisq, df = 1)
    )
}
# get log-rank p values for all non-control treatments
log_rank_pvals_vs_ctrl <- purrr::map(non_control_treatments,
                                            get_log_rank_pval_vs_ctrl,
                                            data = survival_data) |> 
  purrr::list_rbind() |> # combine into one tibble
  dplyr::rename(raw_logrank_p_vs_ctrl = p.value) |> 
  # adjust p values for multiple hypothesis testing
  dplyr::mutate(adjusted_logrank_p_vs_ctrl = p.adjust(raw_logrank_p_vs_ctrl,
                                                      method = "holm"))

# fully pairwise log-rank comparison that includes comparison between treatments
# pairwise_log_rank <- pairwise_survdiff(Surv(time = endpoint_day,
#                                             event = reached_endpoint,
#                                             type = "right") ~ treatment,
#                                        data = survival_data,
#                                        rho = 0)

# report values
write_csv(log_rank_pvals_vs_ctrl,
          fs::path(output_directory,
                   str_glue("log_rank_pvals_vs_ctrl_TGI_{get_timestamp()}.csv")))

# plot custom Kaplan-Meier survival curve by group------------------------------
plot_survival_data <- all_TGI_data |> 
  dplyr::group_by(treatment, day) |> 
  dplyr::summarize(n_living = sum(hit_endpoint == FALSE)) # count living animals
plot_survival_data |> 
  ggplot(aes(x = day, y = n_living, color = treatment)) +
  # horizontal step first means line goes down on day when mouse is missing
  geom_step(direction = "hv",
            linewidth = 1) +
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
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  # scale_fill_manual(values = c("dosing period" = "black"), name = NULL) +
  theme_prism() +
  theme(legend.title = element_text(), # reinstate legend label
        plot.background = element_blank()) + # for transparent background
  # guides(color = guide_legend(order = 1),
  #        fill = guide_legend(order = 2)) +
  labs(x = "days since start of dosing",
       y = "mice below tumor size limit",
       title = "Time to tumor size limit",
       caption = animal_caption)
ggsave(
  str_glue(
    "{output_directory}/plot_TGI_survival_{get_timestamp()}.{plot_type}"
  ),
  # WARNING manual width adjust 2025-04-17
  width = 9, height = 4,
  bg = "transparent")
# function to plot mean/SEM of a type of measurement over time----------------------------
plot_measurement <- function(measurement_data,
                             interval_data, 
                             measurement,
                             y_label,
                             plot_title,
                             width = 8){
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
    scale_x_continuous(limits = c(x_min, x_max),
                       expand = expansion(mult = c(0.05, 0.1))) +
    geom_rect(data = interval_data,
              inherit.aes = FALSE,
              aes(xmin = dosing_start,
                  xmax = dosing_end, #todo: pmin(dosing_end, x_max),
                  ymin = -Inf,
                  ymax = Inf,
                  fill = "dosing period"),
              alpha = 0.1) +
    scale_color_manual(values = color_map_treatments,
                       labels = display_names_treatments) +
    scale_shape_manual(values = shape_map_treatments,
                       labels = display_names_treatments) +
    scale_fill_manual(values = c("dosing period" = "black"), name = NULL) +
    theme_prism() +
    theme(legend.title = element_text(), # reinstate legend label
          plot.background = element_blank()) + # transparent background
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 1),
           fill = guide_legend(order = 2)) +
    labs(x = "days since start of dosing",
         y = y_label,
         title = plot_title,
         caption = animal_caption)
  ggsave(
    str_glue(
      "{output_directory}/plot_TGI_{measurement}_{get_timestamp()}.{plot_type}"
    ),
    width = width, height = 4,
    bg = "transparent")
}
# plot tumor volume and body weight---------------------------------------------
if(!exists("group_plot_width")){group_plot_width <-  8} # default width
plot_measurement(measurement_data = plot_TGI_data,
                 interval_data = dosing_data, 
                 measurement = "volume",
                 y_label = "tumor volume (mm³)",
                 plot_title = "Tumor volume",
                 width = group_plot_width)
plot_measurement(measurement_data = plot_TGI_data,
                 interval_data = dosing_data, 
                 measurement = "body_weight_percent",
                 y_label = "body weight (%)",
                 plot_title = "Body weight",
                 width = group_plot_width)

# function to plot individual traces of measurement-------------------
plot_measurement_individual <- function(measurement_data,
                                        interval_data, 
                                        measurement,
                                        y_label,
                                        plot_title,
                                        y_limits=NULL){
  x_min = min(measurement_data$day)
  x_max = max(measurement_data$day)
  if(is.null(y_limits)){y_limits <- c(NA,NA)}
  plot <- measurement_data |> 
    ggplot(aes(x = day,
               y = .data[[measurement]],
               color = animal,
               shape = animal,
               group = animal)) +
    geom_point(size = 2) +
    geom_line() +
    # dosing period rectangles
    geom_rect(data = interval_data,
              inherit.aes = FALSE,
              aes(xmin = dosing_start,
                  xmax = dosing_end, #todo: pmin(dosing_end, x_max),
                  ymin = -Inf,
                  ymax = Inf,
                  fill = "dosing period"),
              alpha = 0.1) +
    scale_y_continuous(limits = y_limits) +
    scale_color_manual(values = pals::cols25(8)) +
    scale_shape_manual(values = doseplotr::shape_scale_default()) +
    # dosing period rectangle stuff
    scale_fill_manual(values = c("dosing period" = "black")) +
    guides(color = "none", shape = "none",
           fill = guide_legend())+
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
                 y_limits = c(0,3200), # WARNING MAGIC NUMBERS SPECIFIC TO DATA
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
      y_limits = c(77,125), # WARNING MAGIC NUMBERS SPECIFIC TO DATA
      plot_title = str_glue("Body weight, {trt}"))
  ggsave(
    str_glue(
      "{output_directory}/plot_TGI_individual_weight_{trt_name}_{get_timestamp()}.{plot_type}"
    ),
    width = 6, height = 4,
    bg = "transparent")
}