#' Assign Missingness to Proteomics Samples
#'
#' @param data A data frame in long format
#' @param names_col Column containing protein names
#' @param sam_col Column containing samples
#' @param val_col Column containing intensity values
#' @param rep_col Column containing replicates
#' @param ref Reference sample
#' @param MAR_thres Minimum samples for MAR (default 2)
#' @param MNAR_thres Maximum samples for MNAR (default 0)
#'
#' @return A data frame in long format, with missingness assigned for
#'     each comparison
#' @export
#'
#' @examples
#'
assign_missingness <- function(data,
                               names_col,
                               sam_col,
                               val_col,
                               rep_col,
                               ref,
                               MAR_thres = 2,
                               MNAR_thres = 0){

  # change column names
  colnames(data)[colnames(data) == names_col] <- "names"
  colnames(data)[colnames(data) == sam_col] <- "sam"
  colnames(data)[colnames(data) == val_col] <- "val"
  colnames(data)[colnames(data) == rep_col] <- "rep"

  # define replicates
  # define total number of replicates
  replicates <- levels(as.factor(data$rep))
  full <- length(replicates)

  # create copy of data
  # count number of replicates each protein found in
  count <- data
  count <- tidyr::spread(data = count,
                         key = rep,
                         value = val)
  count$count <- rowSums(!is.na(count[,replicates]))

  # keep only count column
  # spread using sample column
  count <- count[,colnames(count)[!(colnames(count) %in% replicates)]]
  samples <- unique(count$sam)
  count <- tidyr::spread(data = count,
                         key = sam,
                         value = count)

  # create vector to store comparison names
  comparisons <- vector(mode = "character")

  # determine missingness
  for (sample_of_interest in samples[samples != ref]){

    # rename columns
    colnames(count)[colnames(count) == sample_of_interest] <- "s"
    colnames(count)[colnames(count) == ref] <- "ref"

    # create new column with outcome of missingness assignment
    count <- dplyr::mutate(.data = count,
                           missingness = dplyr::case_when(s == full & ref == full ~ "complete",
                                                          s < MAR_thres & ref < MAR_thres ~ "missing",
                                                          ((s >= MAR_thres & ref >= MAR_thres) & (s + ref >= full)) ~ "MAR",
                                                          ((s == full & ref > MNAR_thres) | (s > MNAR_thres & ref == full)) ~ "MAR",
                                                          ((s == full & ref <= MNAR_thres) | (s <= MNAR_thres & ref == full)) ~ "MNAR"))

    # revert columns to original names
    # rename missingness column as comparison
    colnames(count)[colnames(count) == "s"] <- sample_of_interest
    colnames(count)[colnames(count) == "missingness"] <- paste(sample_of_interest, "vs", ref, sep="_")
    colnames(count)[colnames(count) == "ref"] <- ref

    # write comparison into vector
    comparisons[sample_of_interest] <- paste(sample_of_interest, "vs", ref, sep="_")
  }

  # keep comparison columns only
  # gather by comparison
  # create new column to separate sample and reference
  count <- count[,colnames(count)[!(colnames(count) %in% samples)]]
  count <- tidyr::gather(data = count,
                         key = "comparison",
                         value = "missingness",
                         all_of(comparisons))
  count[,"sample"] <- count[,"comparison"]
  count <- tidyr::separate(data = count,
                           col = sample,
                           into = c("sample", "reference"),
                           sep="_vs_")
  count <- tidyr::gather(data = count,
                         key = "type",
                         value = "sam",
                         c("sample", "reference"))
  count <- count[,c("names", "comparison", "missingness", "type", "sam")]

  # merge original data with missingness data frame
  # remove comparisons where all values are missing
  data <- merge(data,
                count,
                by = c("names", "sam"),
                all.x = TRUE)
  data <- dplyr::filter(.data = data,
                        missingness != "missing")

  # revert column names
  colnames(data)[colnames(data) == "names"] <- names_col
  colnames(data)[colnames(data) == "sam"] <- sam_col
  colnames(data)[colnames(data) == "rep"] <- rep_col
  colnames(data)[colnames(data) == "val"] <- val_col

  # return output data frame
  data
}
