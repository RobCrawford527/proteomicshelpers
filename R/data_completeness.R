data_completeness <- function(data,
                              samples = NULL,
                              names_col,
                              sam_rep_col,
                              sam_col,
                              comp_col,
                              imp_col,
                              ref,
                              ref_comp,
                              thresholds = c(0.9, 0.8, 0.6)){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "names"
  colnames(data)[colnames(data) == sam_rep_col] <- "sam_rep"
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == comp_col] <- "comp"
  colnames(data)[colnames(data) == imp_col] <- "imp"

  # create samples if not defined already
  if (is.null(samples)){
    samples <- unique(data$sam_rep)
  }

  # filter so rows are unique (remove extra reference rows)
  # keep relevant columns
  # spread data frame - wide format
  quality <- dplyr::filter(.data = data,
                           !(sam == ref & comp != ref_comp))
  quality <- quality[,c("names", "sam_rep", "imp")]
  quality <- tidyr::spread(quality,
                           key = sam_rep,
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
