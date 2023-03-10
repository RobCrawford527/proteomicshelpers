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
    samples <- unique(data[[sample]])
  }

  # define replicates
  replicates <- unique(data[[replicate]])

  # spread data
  data <- tidyr::spread(data = data[,c( {{ protein_id }}, {{ sample }}, {{ replicate }}, {{ value}} )],
                        key = {{ replicate }},
                        value = {{ value }} )

  # create output data frame
  correlations <- tidyr::expand_grid(samples, replicates, replicates)
  rownames(correlations) <- c("sample", "rep_x", "rep_y")
  correlations <- dplyr::filter(correlations,
                                rep_x < rep_y)
  correlations <- dplyr::mutate(correlations,
                                adj_rsq = NA,
                                shared = NA)

  # calculate correlations between replicates
  for (i in 1:nrow(correlations)){

    # determine parameters for row of interest
    s <- correlations[i, "sample"]
    rx <- as.character(correlations[i, "rep_x"])
    ry <- as.character(correlations[i, "rep_y"])

    # filter data
    data_i <- dplyr::filter(data,
                            sample == s)
    data_i <- na.omit(data_i)
    colnames(data_i) <- c("rx", "ry")

    if (nrow(data_i) > 0){
      # calculate adjusted r-squared
      # write into results table
      correlations[i, "adj_rsq"] <- summary(stats::lm( {{ ry }} ~ {{ rx }} , data_i))["adj.r.squared"]
      correlations[i, "shared"] <- nrow(data_i)

    } else if (nrow(data_i) == 0){
      # write into results table
      correlations[i, "adj_rsq"] <- NA
      correlations[i, "shared"] <- 0
    }
  }

  # format output data frame
  correlations[,"rep_x"] <- as.factor(correlations[,"rep_x"])
  correlations[,"rep_y"] <- as.factor(correlations[,"rep_y"])

  # define scale limits if not already
  if (is.null(limits)){
    limits <- c(round(min(correlations[["adj_rsq"]]) - 0.05, 1), 1.0)
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
