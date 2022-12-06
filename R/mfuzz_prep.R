mfuzz_prep <- function(data,
                       thres = 0,
                       replace = FALSE,
                       mode = 'mean',
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
