combine_reps <- function(data,
                         samples = NULL,
                         sam_col,
                         rep_col,
                         val_col,
                         min_reps = 1){

  # change column names
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == rep_col] <- "rep"
  colnames(data)[colnames(data) == val_col] <- "val"

  # create samples if not defined already
  if (is.null(samples)){
    samples <- unique(data$sam)
  }

  # define replicates
  replicates <- levels(as.factor(data$rep))

  # calculate mean for each protein
  data <- spread(data, key = rep, value = val)
  data$mean <- rowMeans(data[,replicates], na.rm = TRUE)

  # identify number of reps each protein is present in
  # remove mean values for proteins present in fewer than min_reps
  data$mean[rowSums(!is.na(data[,replicates])) < min_reps] <- NA

  # add mean to replicates
  # gather replicates
  # keep only mean values (and discard NAs)
  replicates <- c(replicates, "mean")
  data <- tidyr::gather(data,
                        key = rep,
                        value = val,
                        tidyr::all_of(replicates))
  data <- dplyr::filter(data, rep == "mean" & !is.na(val))

  # revert column names
  colnames(data)[colnames(data)=="sam"] <- sam_col
  colnames(data)[colnames(data)=="rep"] <- rep_col
  colnames(data)[colnames(data)=="val"] <- val_col

  # return output data frame: mean values only
  data
}
