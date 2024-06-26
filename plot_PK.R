# plot PK data
# Jack Stevenson started 2024-06-25
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
library(pracma) # for integration
# library(ggthemes) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_PK.R"
source(params_path)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# set up input and output directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
raw_data <- readxl::read_excel(
  str_glue("{input_directory}/{input_filename}"),
  na = "BQL" # "below quantitation limit"
  )
plot_data <- raw_data |> 
  mutate(study = fct_inorder(study),
         dose_mg_kg = fct_inseq(as.factor(dose_mg_kg))) |> 
  replace_na(list(conc_nM = 0, conc_ng_mL = 0))
write_csv(plot_data, fs::path(output_directory,
                              str_glue("PK_raw_data_{get_timestamp()}.csv")))
# function for plotting one study-----------------------------------------------
plot_study <- function(data, specific_study){
  num_doses <- (data |> filter(study == specific_study))$dose_mg_kg |> 
    unique() |> 
    length()
  vr <- doseplotr::viridis_range(num_doses)
  data |>
    filter(study == specific_study) |> 
    group_by(dose_mg_kg, time_hr) |> 
    ggplot(
      aes(x = time_hr,
          y = conc_nM,
          color = dose_mg_kg,
          group = dose_mg_kg)
      ) +
    stat_summary(
      fun.data = "mean_se",
      geom = "errorbar",
      width = 2
    ) +
    stat_summary(
      fun = "mean",
      geom = "point",
      size = 2.5
    ) +
    stat_summary(
      fun = "mean",
      geom = "line"
    ) +
    scale_color_viridis(discrete = TRUE, begin = vr[[1]], end = vr[[2]]) +
    labs(
      x = "time (hours)",
      y = "plasma concentration (nM)",
      title = str_glue("Pharmacokinetics study {study}"),
      color = "dose (mg/kg)"
    ) +
    theme_prism() +
    theme(legend.title = element_text(),
          plot.background = element_blank()) # for transparent background
  ggsave(
    str_glue(
      "{output_directory}/plot_PK_{study}_{get_timestamp()}.{plot_type}"
      ),
    width = 8, height = 6,
    bg = "transparent"
    )
}
# plot each study-----------------------------------------------
studies <- as.vector(unique(plot_data$study))
for (study in studies){
  plot_study(plot_data, study)
}
# plot all studies together by linetype-----------------------------------------------
num_doses <- plot_data$dose_mg_kg |>
  unique() |>
  length()
vr <- doseplotr::viridis_range(num_doses)
plot_data |>
  group_by(study, dose_mg_kg, time_hr) |>
  ggplot(
    aes(
      x = time_hr,
      y = conc_nM,
      color = dose_mg_kg,
      group = interaction(study, dose_mg_kg)
    )
  ) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    width = 2
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 2.5
  ) +
  stat_summary(
    aes(linetype = study),
    fun = "mean",
    geom = "line"
  ) +
  scale_color_viridis(discrete = TRUE, begin = 1, end = 0) +
  # scale_color_viridis(discrete = TRUE, begin = vr[[1]], end = vr[[2]]) +
  # scale_color_colorblind() +
  labs(
    x = "time (hours)",
    y = "plasma concentration (nM)",
    title = str_glue("Pharmacokinetics"),
    color = "dose (mg/kg)"
  ) +
  theme_prism() +
  theme(
    legend.title = element_text(),
    plot.background = element_blank()
  ) # for transparent background
ggsave(
  str_glue(
    "{output_directory}/plot_PK_studies_{get_timestamp()}.{plot_type}"
  ),
  width = 8, height = 6,
  bg = "transparent"
)
# calculate AUC -----------------------------------------------
AUC_data <- plot_data |> 
  summarize(
    mean_conc_nM = mean(conc_nM),
    .by = c(study, dose_mg_kg, time_hr)) |> 
  summarize(
    AUC = pracma::trapz(time_hr, mean_conc_nM) |> signif(3),
    .by = c(study, dose_mg_kg)) |> 
  # temporarily treat dose as numeric to calculate AUC per dose
  mutate(AUC_per_dose = AUC / as.numeric(paste(dose_mg_kg)))
# bar plot components-----------------------------------------------
bar_plot <- function() {
  list(
    geom_bar(stat = "identity", position = "dodge"),
    scale_fill_viridis(
      discrete = TRUE,
      begin = vr_studies[[1]],
      end = vr_studies[[2]]
    ),
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))),
    theme_prism(),
    theme(
      legend.title = element_text(),
      plot.background = element_blank()
    )
  )
}
# plot AUC of each dose-----------------------------------------------
vr_studies <- viridis_range(length(studies))
AUC_data |> 
  ggplot(aes(x = dose_mg_kg, y = AUC, fill = study)) +
  bar_plot() +
  labs(
    x = "dose (mg/kg)",
    y = "AUC(nM•hr)",
    title = "AUC")
ggsave(
  str_glue("{output_directory}/AUC_{get_timestamp()}.{plot_type}"),
  width = 6, height = 6,
  bg = "transparent")
# plot AUC per dose-----------------------------------------------
AUC_data |> 
  ggplot(aes(x = dose_mg_kg, y = AUC_per_dose, fill = study)) +
  bar_plot() +
  labs(
    x = "dose (mg/kg)",
    y = "AUC/dose ratio ((nM•hr)/(mg/kg)",
    title = "dose-dependence of AUC"
    )
ggsave(
  str_glue("{output_directory}/AUC_per_dose_{get_timestamp()}.{plot_type}"),
  width = 6, height = 6,
  bg = "transparent")
