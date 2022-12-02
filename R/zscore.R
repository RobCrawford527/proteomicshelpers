#' Calculate Z-Scores & P-Values Based on a Reference Population
#'
#' @param data A data frame in long format
#' @param samples Vector containing the samples to test (if only a subset of
#'     the possible values)
#' @param comparisons Vector containing the comparisons to test (if only a
#'     subset of the possible values)
#' @param sam_col Column containing samples
#' @param val_col Column containing delta values
#' @param key_col Column containing the comparisons to test
#' @param ref_col Column containing the protein type (used to determine
#'     the reference population)
#' @param reference Indicator of which proteins to use as the reference
#' @param alpha Significance level (default 0.05)
#'
#' @return A data frame containing z-score and p-value columns in addition
#'     to the input columns
#' @export
#'
#' @examples
#'
zscore <- function(data,
                   samples = NULL,
                   comparisons = NULL,
                   sam_col,
                   val_col,
                   key_col,
                   ref_col,
                   reference,
                   alpha = 0.05){

  # change column names
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == val_col] <- "val"
  colnames(data)[colnames(data) == key_col] <- "key"
  colnames(data)[colnames(data) == ref_col] <- "ref"

  # define levels of factors of interest
  samples <- unique(data$sam)
  comparisons <- unique(data$key)

  # create empty output data frame
  output <- data.frame()

  for (s in samples){
    for (c in comparisons){
      # filter for only comparison of interest
      # filter for reference population
      comp <- dplyr::filter(.data = data,
                            sam == s & key == c)
      ref_comp <- dplyr::filter(.data = data,
                                sam == s & key == c & ref == reference)

      # define mean and sd of reference population
      ref_mean <- mean(ref_comp$val, na.rm = TRUE)
      ref_sd <- sd(ref_comp$val, na.rm = TRUE)

      # calculate z-scores
      # calculate raw p-values
      comp[,"zscore"] <- (comp[,"val"] - ref_mean) / ref_sd
      comp[,"pval"] <- ifelse(stats::pnorm(q = comp[,"val"], mean = ref_mean, sd = ref_sd) < 0.5,
                              stats::pnorm(q = comp[,"val"], mean = ref_mean, sd = ref_sd) * 2,
                              (1 - stats::pnorm(q = comp[,"val"], mean = ref_mean, sd = ref_sd)) * 2)

      # FDR correction
      comp <- comp[order(abs(comp$zscore), decreasing = TRUE),]
      comp[,"rank"] <- order(abs(comp$zscore), decreasing = TRUE)
      comp[,"BH"] <- comp[,"rank"] / nrow(comp) * alpha
      comp[,"BH"] <- ifelse(comp[,"pval"] < comp[,"BH"],
                            TRUE,
                            FALSE)

      # adjust BH column
      max <- max(dplyr::filter(.data = comp, BH == TRUE)[,"rank"])
      comp[,"BH"] <- ifelse(comp[,"rank"] <= max,
                            TRUE,
                            FALSE)

      # combine with output data frame
      output <- rbind.data.frame(output, comp)
    }
  }

  # revert column names to original
  colnames(output)[colnames(output) == "sam"] <- sam_col
  colnames(output)[colnames(output) == "val"] <- val_col
  colnames(output)[colnames(output) == "key"] <- key_col
  colnames(output)[colnames(output) == "ref"] <- ref_col

  # return data frame
  output
}
