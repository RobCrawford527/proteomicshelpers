#' Determine & Plot Number of Proteins Per Sample
#'
#' @param data A data frame in long format
#' @param samples List of samples to count (defaults to all)
#' @param names_col Column containing protein names
#' @param sam_col Column containing samples
#' @param val_col Column containing intensity values
#' @param plot Whether to the output should be a plot or not (default FALSE)
#' @param format Format of the sample names, as a vector
#' @param fill Column to use for colouring plot (must be within format)
#'
#' @return A plot and/or data frame showing the number of proteins detected
#'     in each sample
#' @export
#'
#' @examples
#'
protein_count <- function(data,
                          samples = NULL,
                          names_col,
                          sam_col,
                          val_col,
                          plot = FALSE,
                          format = NULL,
                          fill = NULL){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "names"
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == val_col] <- "val"

  # create samples if not defined already
  if (is.null(samples)){
    samples <- unique(data$sam)
  }

  # summarise protein counts
  summary <- data.frame(sam = samples,
                        count = NA)
  for (s in 1:nrow(summary)){
    summary[s, "count"] <- nrow(dplyr::filter(data, sam == summary[s, "sam"] & !is.na(val)))
  }

  # plot if plot == TRUE
  if (plot == TRUE){

    # create copy of sample column
    # split into individual parts
    # define fill column
    summary[,"fill"] <- summary[,"sam"]
    if (!is.null(fill)){
      summary <- tidyr::separate(summary, col = fill, into = format, sep="_")
      colnames(summary)[colnames(summary) == fill] <- "fill"
    }

    # plot
    plot <- ggplot2::ggplot(data = summary,
                            mapping = ggplot2::aes(x = sam,
                                                   y = count,
                                                   fill = fill)) +
      ggplot2::geom_col() +
      ggplot2::geom_text(ggplot2::aes(label = count),
                         angle = 90,
                         nudge_y = -75,
                         vjust = 0.5,
                         hjust = 0.5) +
      ggplot2::scale_y_continuous(limits = c(0, NA),
                                  breaks = seq(0, 5000, 200),
                                  expand = c(0, 0)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                     strip.background = ggplot2::element_blank(), legend.position = 'none')

    # print summary and return plot
    print(summary)
    plot

  } else {

    # return summary
    summary
  }
}
