# simple therapeutic index plots
# load required libraries-------------------------------------------------------
library(tidyverse) # for tidy data handling
library(readxl) # for excel file handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(janitor) # for name repair
library(doseplotr) # you bet
# import potency data-----------------------------------------------------------
potency_data <- read_excel("input/therapeutic window.xlsx")
# relevel factors
treatments = unique(potency_data$treatment)
potency_data <- potency_data |> 
  mutate(treatment = fct_relevel(treatment, treatments))
# calculate therapeutic window--------------------------------------------------
plot_data <- potency_data |> 
  pivot_wider(names_from = target, values_from = IC50_nM) |> 
  mutate(HUVEC_over_T315I = HUVEC/K562_pUltra_T315I)
# bar plot----------------------------------------------------------------------
plot_data |>
  ggplot(aes(y = treatment,
             x = HUVEC_over_T315I,
             label = round(HUVEC_over_T315I, digits = 0))) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(position = position_nudge(x = 2500)) +
  theme_prism() +
  theme(plot.background = element_blank()) + # need for transparent background
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)),
                     label = scales::comma) +
  scale_y_discrete(limits = rev) +
  labs(x = "HUVEC IC50 / K562 pUltra T315I IC50",
       title = "Ratio of cell-killing potency: HUVEC vs. K562")
ggsave("output/therapeutic_window_{get_timestamp()}.pdf",
       bg = "transparent", width = 6, height = 3)