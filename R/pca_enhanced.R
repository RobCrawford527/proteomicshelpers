#' Produce PCA Plots to Compare Proteomics Samples
#'
#' @param data A data frame in long format
#' @param samples List of samples to count (defaults to all)
#' @param names_col Column containing protein names
#' @param sam_col Column containing samples
#' @param val_col Column containing intensity values
#' @param format Format of the sample names, as a vector
#' @param colour_col Column to use for colouring points (defaults to sample,
#'     otherwise must be within format)
#' @param shape_col Column to use for determining shape of points (defaults
#'     to NA, otherwise must be within format)
#' @param label_col Column to use for determining label of points (defaults
#'     to sample, otherwise must be within format)
#'
#' @return A plot showing the first two PCs
#' @export
#'
#' @examples
#'
pca_enhanced <- function(data,
                         samples = NULL,
                         names_col,
                         sam_col,
                         val_col,
                         format = NULL,
                         colour_col = NULL,
                         shape_col = NULL,
                         label_col = NULL){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "names"
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == val_col] <- "val"

  # create samples if not defined already
  if (is.null(samples)){
    samples <- unique(data$sam)
  }

  # keep only columns of interest
  # spread data to wide format
  # name rows with protein names
  data <- data[,c("names", "sam", "val")]
  data <- tidyr::spread(data = data,
                        key = sam,
                        value = val)
  rownames(data) <- data$names

  # transpose and remove missing values
  # print number of proteins included
  data <- t(stats::na.omit(data[,samples]))
  print(paste(ncol(data), " protein(s) included", sep=""))

  # calculate principle components
  pr <- stats::prcomp(x = data,
                      scale = TRUE,
                      retx = TRUE,
                      center = TRUE)

  # calculate percentage of variance explained by PC1 and PC2
  var1 <- paste("PC1 (",
                round((pr$sdev[1]^2) / sum(pr$sdev^2) * 100, 1),
                "% of variance)",
                sep="")
  var2 <- paste("PC2 (",
                round((pr$sdev[2]^2) / sum(pr$sdev^2) * 100, 1),
                "% of variance)",
                sep="")

  # add sample data to pr
  # create colour and shape columns
  pr <- as.data.frame(pr$x)
  pr[,c("sam", "colour", "label")] <- rownames(pr)
  pr[,"shape"] <- as.factor("s")

  # split into individual parts
  # define colour column
  # define shape column
  if (!is.null(colour_col)){
    pr <- tidyr::separate(data = pr,
                          col = colour,
                          into = format,
                          sep="_")
    colnames(pr)[colnames(pr) == colour_col] <- "colour"

    if (!is.null(shape_col)){
      pr$shape <- pr[,shape_col]
    }
    if (!is.null(label_col)){
      pr$label <- pr[,label_col]
    }
  }

  # plot
  plot <- ggplot2::ggplot(data = pr,
                          mapping = ggplot2::aes(x = PC1,
                                                 y = PC2,
                                                 colour = colour,
                                                 shape = shape)) +
    ggplot2::geom_point(size = 5, alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       colour = "black") +
    ggplot2::coord_equal() +
    ggplot2::xlab(var1) +
    ggplot2::ylab(var2) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'right')

  # return plot
  plot
}
