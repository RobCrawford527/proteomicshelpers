#' Calculate Pairwise Correlations Between Replicates
#'
#' @param data A data frame in long format
#' @param samples List of samples to assess (defaults to all)
#' @param protein_id Column containing protein names
#' @param sample Column containing samples
#' @param value Column containing intensity values
#' @param replicate Column containing replicates
#' @param nrow Number of rows for facet_wrap
#' @param ncol Number of columns for facet_wrap
#' @param limits Limits for colour scale
#'
#' @return A plot showing pairwise linear correlations between replicates
#'     of the same sample
#' @export
#'
#' @examples
#'
replicate_correlations <- function(data,
                                   samples = NULL,
                                   protein_id,
                                   sample,
                                   value,
                                   replicate,
                                   nrow = NULL,
                                   ncol = NULL,
                                   limits = NULL){

  # create list of all samples if not defined already
  if (is.null(samples)){
    samples <- unique(dplyr::pull(data, {{ sample }} ))
  }

  # define replicates
  replicates <- unique(dplyr::pull(data, {{ replicate }} ))

  # keep only columns of interest
  data <- dplyr::select(data,
                        {{ protein_id }}, {{ sample }}, {{ replicate }}, {{ value }})

  # create output data frame
  correlations <- as.data.frame(tidyr::expand_grid(samples, replicates, replicates))
  colnames(correlations) <- c("sample", "rep_x", "rep_y")
  correlations <- dplyr::filter(correlations,
                                rep_x < rep_y)
  correlations <- dplyr::mutate(correlations,
                                adj_rsq = NA,
                                shared = NA)

  # calculate correlations between replicates
  for (i in 1:nrow(correlations)){

    # determine parameters for row of interest
    s <- dplyr::pull(correlations, sample)[i]
    rx <- dplyr::pull(correlations, rep_x)[i]
    ry <- dplyr::pull(correlations, rep_y)[i]

    # filter data
    data_rx <- dplyr::filter(data,
                             {{ sample }} == s & {{ replicate }} == rx)
    data_ry <- dplyr::filter(data,
                             {{ sample }} == s & {{ replicate }} == ry)

    # check that data from both replicates are present
    # if not, proceed directly to writing into table
    if (nrow(data_rx) > 0 & nrow(data_ry) > 0){

      # combine data_rx and data_ry
      # spread data
      # remove proteins with missing values
      data_i <- rbind.data.frame(data_rx, data_ry)
      data_i <- tidyr::spread(data_i,
                              key = {{ replicate }},
                              value = {{ value }})
      data_i <- na.omit(data_i[,c(rx, ry)])
      colnames(data_i) <- c("rx", "ry")

      # check that at least one protein with values for both replicates remains
      # if not, proceed directly to writing into table
      if (nrow(data_i) > 0){
        # calculate adjusted r-squared
        # write into results table
        correlations[i,"adj_rsq"] <- summary(stats::lm(ry ~ rx, data_i))["adj.r.squared"]
        correlations[i,"shared"] <- nrow(data_i)

      } else if (nrow(data_i) == 0){
        # write into results table
        correlations[i,"adj_rsq"] <- NA
        correlations[i,"shared"] <- 0
      }

    } else {
      # write into results table
      correlations[i,"adj_rsq"] <- NA
      correlations[i,"shared"] <- NA
    }
  }

  # format output data frame
  correlations[,"rep_x"] <- as.factor(correlations[,"rep_x"])
  correlations[,"rep_y"] <- as.factor(correlations[,"rep_y"])

  # define scale limits if not already
  if (is.null(limits)){
    limits <- c(round(min(correlations[,"adj_rsq"]) - 0.05, 1), 1.0)
  }

  # plot correlations
  plot <- ggplot2::ggplot(data = correlations,
                          mapping = ggplot2::aes(x = rep_x,
                                                 y = rep_y,
                                                 fill = adj_rsq)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(mapping = ggplot2::aes(label = signif(adj_rsq, 2),
                                              fill = NULL),
                       nudge_y = 0.15,
                       colour = "white") +
    ggplot2::geom_text(mapping = ggplot2::aes(label = shared,
                                              fill = NULL),
                       nudge_y = -0.15,
                       colour = "white") +
    viridis::scale_fill_viridis(discrete = FALSE,
                                begin = 0.15,
                                end = 0.85,
                                direction = -1,
                                option = "inferno",
                                limits = limits,
                                breaks = seq(0, 1, 0.1)) +
    ggplot2::facet_wrap(dplyr::vars(sample), nrow = nrow, ncol = ncol) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
                   strip.background = ggplot2::element_blank())

  # return plot
  plot
}
