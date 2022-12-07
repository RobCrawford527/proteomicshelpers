alpha_core_modified <- function(eset,
                                clusters,
                                min.acore = 0.5c,
                                reference,
                                names_col){
  
  # extract alpha cores
  cores <- Mfuzz::acore(eset = eset,
                        cl = clusters,
                        min.acore = min.acore)
  
  # add row names to cluster cores and append into single data frame
  cores_combined <- data.frame()
  for (i in 1:length(cores)){
    cores[[i]][["Cluster"]] <- i
    cores_combined <- rbind(cores_combined,
                            cores[[i]])
  }
  
  # merge cores_combined with reference
  cores_combined <- merge(reference,
                          cores_combined,
                          by.x = names_col,
                          by.y = "NAME",
                          all.y = TRUE)
  
  # return list of cluster cores
  cores_combined
}
