#' Extract Cluster Cores From Fuzzy Clustering
#'
#' @param eset An expression set used for fuzzy clustering
#' @param clusters A set of clusters from fuzzy clustering, with protein names
#'     and membership values
#' @param min.acore Cluster membership threshold, i.e. the minimum membership
#'     value for proteins to be included in the alpha core for each cluster
#' @param reference Reference table for converting protein names
#' @param names_col Column in reference table to use for merging
#'
#' @return A data frame containing the proteins assigned to the alpha core
#'     for each cluster
#' @export
#'
#' @examples
#'
alpha_core_modified <- function(eset,
                                clusters,
                                min.acore = 0.5,
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
