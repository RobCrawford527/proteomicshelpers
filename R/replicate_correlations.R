#' Calculate Pairwise Correlations Between Replicates
#'
#' @param data A data frame in long format
#' @param samples List of samples to assess (defaults to all)
#' @param sam_col Column containing samples
#' @param val_col Column containing intensity values
#' @param rep_col Column containing replicates
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
                                   names_col,
                                   sam_col,
                                   val_col,
                                   rep_col,
                                   nrow = NULL,
                                   ncol = NULL){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == val_col] <- "val"
  colnames(data)[colnames(data) == rep_col] <- "rep"

  # create samples if not defined already
  if (is.null(samples)){
    samples <- unique(data$sam)
  }

  # define replicates
  replicates <- levels(as.factor(data$rep))

  # spread data
  data <- tidyr::spread(data = data[,c("name", "sam", "rep", "val")],
                        key = rep,
                        value = val)

  # create output data frame
  correlations <- data.frame()
  for (s in samples){
    for (r in replicates){
      correlations <- rbind.data.frame(correlations,
                                       data.frame(sample = s,
                                                  RepX = r,
                                                  RepY = replicates,
                                                  adjrsq = NA,
                                                  shared = NA,
                                                  limits = NULL))
    }
  }
  correlations <- dplyr::filter(correlations, RepX < RepY)

  # calculate correlations between replicates
  for (i in 1:nrow(correlations)){

    # determine parameters for row of interest
    s <- correlations[i, "sample"]
    rx <- as.character(correlations[i, "RepX"])
    ry <- as.character(correlations[i, "RepY"])

    # filter data
    data_i <- dplyr::filter(data, sam == s)
    data_i <- na.omit(data_i[,c(rx, ry)])
    colnames(data_i) <- c("rx", "ry")

    if (nrow(data_i) > 0){
      # calculate adjusted r-squared
      # write into results table
      correlations[i, "adjrsq"] <- summary(stats::lm(ry ~ rx, data = data_i))["adj.r.squared"]
      correlations[i, "shared"] <- nrow(data_i)

    } else if (nrow(data_i) == 0){
      # write into results table
      correlations[i, "adjrsq"] <- NA
      correlations[i, "shared"] <- 0
    }
  }

  # format output data frame
  correlations[,"RepX"] <- as.factor(correlations[,"RepX"])
  correlations[,"RepY"] <- as.factor(correlations[,"RepY"])

  # define scale limits if not already
  if (is.null(limits)){
    limits <- c(round(min(correlations$adjrsq)-0.05, 1), 1.0)
  }

  # plot correlations
  plot <- ggplot2::ggplot(data = correlations,
                          mapping = ggplot2::aes(x = RepX,
                                                 y = RepY,
                                                 fill = adjrsq)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(mapping = ggplot2::aes(label = signif(adjrsq,2),
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
