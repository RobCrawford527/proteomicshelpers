#' Determine & Plot Number of Proteins Per Sample
#'
#' @param data A data frame in long format
#' @param samples List of samples to count (defaults to all)
#' @param protein_id Column containing protein names
#' @param sample Column containing samples
#' @param value Column containing intensity values
#' @param plot Whether to the output should be a plot or not (default FALSE)
#' @param fill Column to use for colouring plot
#'
#' @return A plot and/or data frame showing the number of proteins detected
#'     in each sample
#' @export
#'
#' @examples
#'
protein_count <- function(data,
                          samples = NULL,
                          protein_id,
                          sample,
                          value,
                          plot = FALSE,
                          fill){

  # filter data frame for samples of interest (if appropriate)
  if (!is.null(samples)){
    data <- dplyr::filter(data,
                          {{ sample }} %in% samples)
  }

  # group data frame
  # summarise protein counts
  counts <- dplyr::group_by(data, {{ sample }} )
  counts <- dplyr::summarise(counts, count = dplyr::n() )

  # plot if appropriate
  if (plot == TRUE){

    # add fill column
    counts <- dplyr::left_join(counts,
                               dplyr::distinct(dplyr::select(data,
                                                             {{ sample }} , {{ fill }} )),
                               by = dplyr::join_by( {{ sample }} ),
                               keep = FALSE)

    # create plot
    plot <- ggplot2::ggplot(data = counts,
                            mapping = ggplot2::aes(x = {{ sample }} ,
                                                   y = count,
                                                   fill = {{ fill }} )) +
      ggplot2::geom_col() +
      ggplot2::geom_text(ggplot2::aes(label = count),
                         angle = 90,
                         nudge_y = -75,
                         vjust = 0.5,
                         hjust = 0.5) +
      ggplot2::scale_y_continuous(limits = c(0, NA),
                                  breaks = seq(0, 10000, 200),
                                  expand = c(0, 0)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                     strip.background = ggplot2::element_blank(), legend.position = "none")

    # print summary and return plot
    print(counts)
    plot

  } else {

    # return summary
    counts
  }
}
