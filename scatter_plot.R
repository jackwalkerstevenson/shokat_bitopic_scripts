scatter_plot <- function(data, pt_size = 3, alpha = 0.7,
                         viridis_begin = 1, viridis_end = 0,
                         width = 7, height = 7,
                         plot_name = "scatter_plot"){
  data %>%
    group_by(compound, target) %>%
    inhibition_summarize() %>%
    ggplot(aes(x = target, color = compound)) +
    geom_point(aes(y = mean_pct_inhibition), size = pt_size, alpha = alpha) +
    # geom_line(aes(y = mean_pct_inhibition)) +
    geom_errorbar(aes(ymax = mean_pct_inhibition+sem,
                      ymin = mean_pct_inhibition-sem,
                      width = w),
                  alpha = alpha) +
    scale_color_viridis(discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) + # opaque legend
    theme_prism() +
    # rotated, right-justified x labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Single-point SelectScreen inhibition",
         caption = "Note: concentrations not equal between kinases",
         x = "target kinase",
         y = "percent inhibition")
  ggsave(str_glue("output/{plot_name}.{plot_type}"),
         bg = "transparent", width = 7, height = 7)
}