#' Prepare Data for Fuzzy Clustering
#'
#' @param data A data frame in wide format. Columns should be samples, and rows
#'     should be proteins.
#' @param thres The proportion of missing values that is allowed (default 0)
#' @param replace Whether NAs should be replaced
#' @param mode The method for replacing NAs (default "mean")
#' @param ...
#'
#' @return A standardised expression set ready for fuzzy clustering using Mfuzz
#' @export
#'
#' @examples
#'
mfuzz_prep <- function(data,
                       thres = 0,
                       replace = FALSE,
                       mode = "mean",
                       ...){

  # convert data to expression set for Mfuzz analysis
  eset <- Biobase::ExpressionSet(assayData = as.matrix(data))

  # filter for and replace missing data
  eset <- Mfuzz::filter.NA(eset,
                           thres = thres)
  if (replace == TRUE){
    eset <- Mfuzz::fill.NA(eset,
                           mode = mode)
  }

  # standardise data (mean = 0, sd = 1)
  eset <- Mfuzz::standardise(eset)

  # return standardised expression set
  eset
}
