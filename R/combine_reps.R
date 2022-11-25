combine_reps <- function(data,
                         samples = NULL,
                         key_col,
                         value_col,
                         sample_format,
                         min_reps = 1){

  # create samples vector if not defined already
  # create list of sample columns
  # define columns
  if (is.null(samples)){
    samples <- unique(data[,key_col])
  }
  columns <- strsplit(sample_format, split="_")[[1]]

  # rename columns
  colnames(data)[colnames(data) == key_col] <- "key_col"
  colnames(data)[colnames(data) == value_col] <- "value_col"

  # format data frame
  output <- separate(data,col=key_col,into=columns,sep="_")
  output <- unite(output,col=key_col,columns[columns!="Replicate"],sep="_")

  # define replicates
  replicates <- levels(as.factor(output[,"Replicate"]))

  # calculate mean for each protein
  output <- spread(output,key=Replicate,value=value_col)
  output$mean <- rowMeans(output[,replicates],na.rm=TRUE)

  # identify number of reps each protein is present in
  output$count <- rowSums(!is.na(output[,replicates]))

  # remove any proteins present in fewer than minimum reps
  output[,"mean"][output[,"count"] < min_reps] <- NA

  # remove count column and format data frame
  # add mean to other replicates
  output <- subset(output,select=-count)
  replicates <- c(replicates,"mean")
  output <- gather(output,key="Replicate",value=value_col,all_of(replicates))
  output <- separate(output,col=key_col,into=columns[columns!="Replicate"],sep="_")
  output <- unite(output,col=key_col,columns,sep="_")

  # rename columns
  colnames(output)[colnames(output)=="key_col"] <- key_col
  colnames(output)[colnames(output)=="value_col"] <- value_col

  # return output data frame
  output
}
