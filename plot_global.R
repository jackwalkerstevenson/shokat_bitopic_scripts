plot_global <- function(plot){
  plot +
    scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                       breaks = scales::breaks_width(1),
                       minor_breaks = minor_x) + # manual minor ticks
    scale_y_continuous(guide = "prism_offset",  # end at last tick
                       breaks = c(0,25,50,75,100)) + # manual y axis ticks
    coord_cartesian(xlim = x_limits, # set x axis zoom from global values
                    ylim = c(0,NA)) + # set y axis zoom locally
    theme_prism() + # make it look fancy like prism
    theme(plot.background = element_blank()) + # need for transparent background
    labs(x = "log [compound] (M)")
}