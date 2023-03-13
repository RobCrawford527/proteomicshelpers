#' Produce PCA Plots to Compare Proteomics Samples
#'
#' @param data A data frame in long format
#' @param samples List of samples to count (defaults to all)
#' @param protein_id Column containing protein names
#' @param sample Column containing samples
#' @param value Column containing intensity values
#' @param colour Column to use for colouring points (defaults to sample,
#'     otherwise must be within format)
#' @param shape Column to use for determining shape of points (defaults
#'     to NA, otherwise must be within format)
#' @param label Column to use for determining label of points (defaults
#'     to sample, otherwise must be within format)
#' @param normalise Logical indicating whether or not to normalise values for plotting
#' @param plot_label Logical indicating whether or not to add labels to plot
#'
#' @return A plot showing the first two PCs
#' @export
#'
#' @examples
#'
pca_enhanced <- function(data,
                         samples = NULL,
                         protein_id,
                         sample,
                         value,
                         colour = NULL,
                         shape = NULL,
                         label = NULL,
                         normalise = FALSE,
                         plot_label = TRUE){

  # create list of all samples if not defined already
  if (is.null(samples)){
    samples <- unique(dplyr::pull(data, {{ sample }} ))
  }

  # keep only columns of interest
  # spread data to wide format
  # name rows with protein names
  data_pc <- dplyr::select(data,
                           {{ protein_id }} , {{ sample }} , {{ value }} )
  data_pc <- tidyr::spread(data = data_pc,
                           key = {{ sample }} ,
                           value = {{ value }} )

  # transpose and remove missing values
  # print number of proteins included
  data_pc <- t(stats::na.omit(data_pc[,samples]))
  print(paste(ncol(data_pc), " protein(s) included", sep = ""))

  # calculate principle components
  pr <- stats::prcomp(x = data_pc,
                      scale = TRUE,
                      retx = TRUE,
                      center = TRUE)

  # calculate percentage of variance explained by PC1 and PC2
  var1 <- paste("PC1 (",
                round((pr$sdev[1] ^ 2) / sum(pr$sdev ^ 2) * 100, 1),
                "% of variance)",
                sep = "")
  var2 <- paste("PC2 (",
                round((pr$sdev[2] ^ 2) / sum(pr$sdev ^ 2) * 100, 1),
                "% of variance)",
                sep = "")

  # extract pr data frame
  # normalise values if appropriate
  pr <- as.data.frame(pr$x)
  pr[,"sample"] <- rownames(pr)
  if (normalise == TRUE){
    pr[,c("PC1", "PC2")] <- apply(pr[,c("PC1", "PC2")], 2, function(x) (x - mean(x)) / sd(x))
  }

  # add colour, shape and label columns
  if (plot_label == TRUE){
    pr <- dplyr::left_join(pr,
                           dplyr::distinct(dplyr::select(data,
                                                         {{ sample }} , {{ colour }} , {{ shape }} , {{ label }} )),
                           by = dplyr::join_by( sample == {{ sample }} ),
                           keep = FALSE)
  } else {
    pr <- dplyr::left_join(pr,
                           dplyr::distinct(dplyr::select(data,
                                                         {{ sample }} , {{ colour }} , {{ shape }} )),
                           by = dplyr::join_by( sample == {{ sample }} ),
                           keep = FALSE)
  }

  # plot
  plot <- ggplot2::ggplot(data = pr,
                          mapping = ggplot2::aes(x = PC1,
                                                 y = PC2,
                                                 colour = {{ colour }} ,
                                                 shape = {{ shape }} )) +
    ggplot2::geom_point(size = 5, alpha = 0.8) +
    ggplot2::coord_equal() +
    ggplot2::xlab(var1) +
    ggplot2::ylab(var2) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")

  # add labels to plot if appropriate
  if (plot_label == TRUE){
    plot <- plot +
      ggplot2::geom_text(ggplot2::aes(label = {{ label }} ),
                         colour = "black")
  }

  # return plot
  plot
}
