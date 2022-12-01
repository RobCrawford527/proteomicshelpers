#' Sum Up Fractions to Calculate Ribosome Engagement
#'
#' @param data A data frame in long format. The "fractions" must be present
#'     in a separate column to the other sample information.
#' @param val_col Column containing intensity values
#' @param key_col Column containing fractions
#' @param ribosome_fractions The ribosome fractions to sum up
#' @param total The "total" fraction to use for normalisation
#' @param normalise Whether or not to normalise the ribosome_engagement score
#'     to the total (default is TRUE)
#'
#' @return A data frame containing (normalised) ribosome engagement scores
#'     and totals
#' @export
#'
#' @examples
#'
ribosome_engagement <- function(data,
                                val_col,
                                key_col,
                                ribosome_fractions,
                                total = NULL,
                                normalise = TRUE){

  # change column names
  colnames(data)[colnames(data) == val_col] <- "val"
  colnames(data)[colnames(data) == key_col] <- "key"

  # define all values of key column
  # spread data using key column
  # reverse log2 transform
  fractions <- unique(data$key)
  data <- tidyr::spread(data = data,
                        key = key,
                        value = val)
  data[,ribosome_fractions] <- 2 ^ data[,ribosome_fractions]

  # calculate ribosome engagement for each protein
  data$ribosome.engagement <- rowSums(data[,ribosome_fractions], na.rm = TRUE)

  # log2 transform
  # replace infinite values with NAs
  data[,c(ribosome_fractions, "ribosome.engagement")] <- log2(data[,c(ribosome_fractions, "ribosome.engagement")])
  data$ribosome.engagement[is.infinite(data$ribosome.engagement)] <- NA

  # normalise to totals
  if (normalise == TRUE){
    data[,c(ribosome_fractions, "ribosome.engagement")] <- data[,c(ribosome_fractions, "ribosome.engagement")] - data[,total]
  }

  # filter to keep only total and ribosome engagement columns
  columns <- colnames(data)[!(colnames(data) %in% ribosome_fractions)]
  data <- data[,columns]

  # revert column names
  colnames(data)[colnames(data) == "val"] <- val_col
  colnames(data)[colnames(data) == "key"] <- key_col
  colnames(data)[colnames(data) == total] <- "total"

  # return output data frame
  data
}
