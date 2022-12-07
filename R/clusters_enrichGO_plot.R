clusters_enrichGO_plot <- function(go_result,
                                   save = FALSE,
                                   filename = NULL,
                                   width = 160){
  
  # find min and max p.adjust values
  # round to whole number
  min_p <- round(min(go_result$p.adjust) - 0.5, 0)
  max_p <- round(max(go_result$p.adjust) + 0.5, 0)
  # adjust min_p if too low
  if (min_p < 1){
    min_p <- 1
  }
  # define limits and breaks
  limits <- c(1, max_p)
  breaks <- unique(round(seq(min_p+1, max_p-1, (max_p-min_p-2)/3), 0))
  
  # plot
  go_plot <- ggplot2::ggplot(data = go_result,
                             mapping = ggplot2::aes(x = Cluster,
                                                    y = Description,
                                                    colour = p.adjust,
                                                    size = p.adjust,
                                                    group = Description)) +
    ggplot2::geom_point() +
    viridis::scale_colour_viridis(option = "plasma",
                                  discrete = FALSE,
                                  begin = 0,
                                  end = 0.85,
                                  direction = -1,
                                  limits = limits,
                                  breaks = breaks) +
    ggplot2::scale_size(limits = limits,
                        breaks = breaks) +
    ggplot2::scale_x_discrete(expand = c(0, 1)) +
    ggplot2::scale_y_discrete(expand = c(0, 1)) +
    ggplot2::facet_grid(ONTOLOGY ~ .,
                        scales = "free_y",
                        space = "free_y") +
    ggplot2::xlab(label = "Cluster") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(colour = "black"),
                   axis.text.y = ggplot2::element_text(colour = "black",
                                                       size = 8),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank(),
                   legend.position = "right")
  
  # save if appropriate
  if (save == TRUE){
    ggplot2::ggsave(filename = filename,
                    plot = go_plot,
                    device = "svg",
                    height = nrow(go_result)*3.5+20,
                    width = width,
                    units = "mm")
  }
  
  # return plot
  go_plot
}
