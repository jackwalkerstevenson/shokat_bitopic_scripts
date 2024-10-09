# plot PK data post dosing from a chronic-dose experiment
# Jack Stevenson started 2024-09-30
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
library(readxl) # for excel import
# library(ggthemes) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_PK_chronic.R"
source(params_path)
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# set up input and output directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
raw_data <- readxl::read_excel(
  str_glue("{input_directory}/{input_filename}")
  )
plot_data <- raw_data |> 
  mutate(dose_mg_kg = fct_inseq(as.factor(dose_mg_kg))) |> 
  mutate(days_since_dose = hr_since_dose / 24)
write_csv(plot_data, fs::path(output_directory,
                              str_glue("PK_raw_data_{get_timestamp()}.csv")))
# plot PK by timepoint-----------------------------------------------
num_doses <- plot_data$dose_mg_kg |> 
  unique() |> 
  length()
vr <- doseplotr::viridis_range(num_doses)
plot_data |> 
  ggplot(aes(x = days_since_dose,
             y = conc_nM,
             color = dose_mg_kg,
             group = dose_mg_kg)) +
  stat_summary(
    fun.data = "mean_se",
    geom = "errorbar",
    width = .5
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
  scale_x_continuous(breaks = scales::breaks_width(1))+
  scale_y_continuous(transform = "identity", breaks = scales::breaks_width(1000))+
  scale_color_viridis(discrete = TRUE, begin = vr[[1]], end = vr[[2]]) +
  labs(
    x = "time since last dose (days)",
    y = "plasma concentration (nM)",
    title = "PK after chronic dosing",
    color = "dose (mg/kg)",
    caption = "n=3, 14d QD IP"
  ) +
  theme_prism() +
  theme(legend.title = element_text(), # reenable legend title
        panel.grid = element_line(linetype = "dotted", linewidth = 0.4),
        panel.grid.minor = element_line(linetype = "dotted", linewidth = 0.1),
        plot.background = element_blank()) # for transparent background
ggsave(
  str_glue(
    "{output_directory}/plot_PK_chronic_{get_timestamp()}.{plot_type}"
  ),
  width = 8, height = 6,
  bg = "transparent"
)