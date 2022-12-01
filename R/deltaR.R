deltaR <- function(data,
                   val_col,
                   key_col,
                   delta = NULL){
  
  # change column names
  colnames(data)[colnames(data) == val_col] <- "val"
  colnames(data)[colnames(data) == key_col] <- "key"
  
  # define levels of factor of interest
  if (is.null(delta)){
    delta <- unique(data$key)
  }
  
  # spread data using key column
  data <- tidyr::spread(data = data,
                        key = key,
                        value = val)
  
  # create vector to store comparison names
  # set index
  comparisons <- vector(mode = "character")
  index <- 1
  
  # perform pairwise comparisons
  # only in one direction
  for (i in 1:(length(delta)-1)){
    for (j in 2:length(delta)) {
      if (i < j){
        # calculate deltaR
        data[,"delta"] <- data[,delta[j]] - data[,delta[i]]
        
        # rename delta column
        # write comparison into vector
        # iterate index
        colnames(data)[colnames(data) == "delta"] <- paste(delta[j], "vs", delta[i], sep="_")
        comparisons[index] <- paste(delta[j], "vs", delta[i], sep="_")
        index <- index + 1
      }
    }
  }
  
  # keep comparison columns only
  # gather by comparison
  # create new column to separate sample and reference
  data <- data[,colnames(data)[!(colnames(data) %in% delta)]]
  data <- tidyr::gather(data = data,
                        key = "comparison",
                        value = "deltaR",
                        dplyr::all_of(comparisons))
  
  # return output data frame
  data
}
