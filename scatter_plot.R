scatter_plot <- function(data,
                         pt_size = 3, alpha = 0.7,
                         viridis_begin = 1, viridis_end = 0,
                         width = 7, height = 7, dodge_width = .2,
                         plot_name = "scatter_plot"){
  data %>%
    inhibition_summarize() %>%
    ggplot(aes(y = target, color = compound)) +
    geom_point(aes(x = mean_pct_inhibition), size = pt_size, alpha = alpha,
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(xmax = mean_pct_inhibition+sem,
                      xmin = mean_pct_inhibition-sem,
                      width = bar_size),
                  position = position_dodge(width = dodge_width),
                  alpha = alpha) +
    scale_color_viridis(discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) + # opaque legend
    theme_prism() +
    # rotated, right-justified x labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Single-point SelectScreen inhibition",
         caption = "Note: compound concentrations not equal between kinases",
         x = "percent inhibition",
         y = "target kinase")
  ggsave(str_glue("output/{plot_name}.{plot_type}"),
         bg = "transparent", width = width, height = height)
}