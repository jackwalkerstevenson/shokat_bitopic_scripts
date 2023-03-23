# plot a scatter plot grouped by target
# each group of treatments is dodged top to bottom for visibility without overlap
scatter_plot <- function(data,
                         pt_size = 3, alpha = 0.7,
                         viridis_begin = 1, viridis_end = 0,
                         width = 7, height = 7, dodge_width = -.2,
                         file_name = "scatter_plot",
                         title = "[plot title here]",
                         caption = NULL,
                         xlab = "[x label here]",
                         ylab = "[y label here]"
                         ){
  data %>%
    inhibition_summarize() %>%
    ggplot(aes(y = target, color = treatment)) +
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
    theme(plot.background = element_blank()) + # need for transparent background
    labs(title = title,
         caption = caption,
         x = xlab,
         y = ylab)
  ggsave(str_glue("output/{file_name}.{plot_type}"),
         bg = "transparent", width = width, height = height)
}