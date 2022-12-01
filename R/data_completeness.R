#' Evaluate Completeness of Proteomics Data
#'
#' @param data A data frame in long format
#' @param samples List of samples to assess (defaults to all)
#' @param names_col Column containing protein names
#' @param sam_col Column containing samples (including replicates)
#' @param imp_col Column containing whether value is imputed or not
#' @param thresholds Thresholds defining quality categories (must be three
#'     values; default 0.9, 0.8, 0.6)
#'
#' @return A data frame, indicating the proportion of sample values that
#'     are imputed for each protein
#' @export
#'
#' @examples
#'
data_completeness <- function(data,
                              samples = NULL,
                              names_col,
                              sam_col,
                              imp_col,
                              thresholds = c(0.9, 0.8, 0.6)){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "names"
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == imp_col] <- "imp"

  # create samples if not defined already
  if (is.null(samples)){
    samples <- unique(data$sam)
  }

  # create copy of data frame and keep relevant columns
  # spread data frame - wide format
  quality <- data[,c("names", "sam", "imp")]
  quality <- tidyr::spread(quality,
                           key = sam,
                           value = imp)

  # determine how many values are present for protein
  # determine how many values are real (not imputed) for protein
  # calculate proportion of real values
  quality[,"values"] <- rowSums(!is.na(quality[,samples]), na.rm = TRUE)
  quality[,"not_imputed"] <- rowSums(quality[,samples] == FALSE, na.rm = TRUE)
  quality[,"completeness"] <- round(quality[,"not_imputed"] / quality[,"values"], 3)

  # assign quality rating
  quality <- dplyr::mutate(quality,
                           quality = dplyr::case_when(completeness == 1 ~ "complete",
                                                      completeness >= thresholds[1] ~ "fantastic",
                                                      completeness >= thresholds[2] ~ "good",
                                                      completeness >= thresholds[3] ~ "ok",
                                                      TRUE ~ "poor"))

  # keep columns of interest
  # revert column names
  quality <- quality[,c("names", "completeness", "quality")]
  colnames(quality)[colnames(quality) == "names"] <- names_col

  # return quality data frame
  quality
}
