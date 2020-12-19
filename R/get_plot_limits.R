get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(ymin = ymin, ymax = ymax)
}
