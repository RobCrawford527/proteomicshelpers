plot_fuzzy_clusters <- function(eset,
                                clusters,
                                samples,
                                plot_background = TRUE,
                                plot_centres = TRUE,
                                mem_thres = 0.5,
                                lwd = 0.75,
                                major_x = NULL,
                                minor_x = NULL,
                                nrow = NULL,
                                ncol = NULL,
                                ...){
  
  # convert eset to data frame
  exprs <- as.data.frame(Biobase::exprs(eset))
  exprs$protein <- row.names(exprs)
  
  # convert membership values to data frame
  mem <- as.data.frame(clusters[["membership"]])
  mem$protein <- row.names(mem)
  
  # convert cluster centres to data frame
  # keep cluster levels
  cen <- as.data.frame(clusters[["centers"]])
  cen[,"protein"] <- row.names(cen)
  cen[,"cluster"] <- as.factor(as.numeric(row.names(cen)))
  clus <- levels(cen[,"cluster"])
  cen <- tidyr::gather(data = cen,
                       key = "sample",
                       value = "intensity",
                       dplyr::all_of(samples))
  cen[,"sample"] <- factor(x = cen[,"sample"],
                           levels = samples)
  
  # merge expression and membership values
  exprs <- merge(exprs,
                 mem,
                 by = "protein",
                 all = TRUE)
  
  # check number of columns
  # tidy data to allow plotting with ggplot
  if (ncol(exprs) != 1 + length(samples) + length(clus)){
    stop("incorrect number of columns")
  }
  exprs <- tidyr::gather(data = exprs,
                         key = "cluster",
                         value = "membership",
                         dplyr::all_of(clus))
  exprs <- tidyr::gather(data = exprs,
                         key = "sample",
                         value = "intensity",
                         dplyr::all_of(samples))
  
  # convert sample and cluster columns to factor
  exprs[,"sample"] <- factor(x = exprs[,"sample"],
                             levels = samples)
  exprs[,"cluster"] <- as.factor(as.numeric(exprs[,"cluster"]))
  
  # create basic structure for plot
  fuzzy_plot <- ggplot2::ggplot(data = exprs,
                                mapping = ggplot2::aes(x = sample,
                                                       y = intensity,
                                                       colour = membership,
                                                       group = protein))
  
  # plot background data (i.e. proteins that don't belong to a given cluster)
  if (plot_background == TRUE){
    fuzzy_plot <- fuzzy_plot +
      ggplot2::geom_line(data = dplyr::filter(.data = exprs,
                                              membership < mem_thres),
                         lwd = lwd,
                         colour = 'grey75',
                         alpha = 0.15)
  }
  
  # plot main data for each cluster
  # layer up the data in 0.05 increments ...
  # ... so that proteins with the highest membership are on top
  div <- seq(mem_thres, 1, 0.05)
  for (d in div){
    fuzzy_plot <- fuzzy_plot +
      ggplot2::geom_line(data = dplyr::filter(.data = exprs,
                                              membership >= d & membership < d+0.05),
                         lwd = lwd,
                         alpha = 0.9)
  }
  
  # add in additional plot elements
  fuzzy_plot <- fuzzy_plot +
    
    # divisions between groups of samples
    ggplot2::geom_vline(xintercept = minor_x,
                        colour = "grey60",
                        lty = 2) +
    ggplot2::geom_vline(xintercept = major_x,
                        colour = "grey25",
                        lty = 2) +
    # axis labels
    ggplot2::scale_x_discrete(labels = samples) +
    ggplot2::scale_y_continuous(breaks = seq(-4, 4, 1)) +
    # colour scale
    viridis::scale_colour_viridis(begin = 0.2, end = 0.9,
                                  direction = -1,
                                  option = "inferno",
                                  breaks = seq(mem_thres, 1, 0.1),
                                  limits = c(mem_thres, 1)) +
    # facetting
    ggplot2::facet_wrap(cluster ~ .,
                        as.table = FALSE,
                        nrow = nrow,
                        ncol = ncol) +
    # theme
    ggplot2::theme_classic() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = 0.5))
  
  # plot cluster centres
  if (plot_centres == TRUE){
    fuzzy_plot <- fuzzy_plot +
      ggplot2::geom_line(data = cen,
                         colour = "black",
                         lwd = lwd*1.5)
  }
  
  # return plot
  fuzzy_plot
}
