#' Assign Missingness to Proteomics Samples
#'
#' @param data A data frame in long format
#' @param name Column containing protein names
#' @param sample Column containing samples
#' @param experiment Column containing experiments
#' @param replicate Column containing replicates
#' @param value Column containing values
#' @param reference Reference sample
#' @param MAR_thres Minimum samples for MAR (default 2)
#' @param MNAR_thres Maximum samples for MNAR (default 0)
#'
#' @return A data frame in long format, with missingness assigned for each comparison
#' @export
#'
assign_missingness <- function(data,
                               name,
                               sample,
                               experiment,
                               replicate,
                               value,
                               reference,
                               MAR_thres = 2,
                               MNAR_thres = 0){

  # define replicates
  # define total number of replicates
  data <- dplyr::ungroup(data)
  replicates <- dplyr::select(data, {{ replicate }} )
  replicates <- unlist(unique(replicates))
  full <- length(replicates)

  # create copy of data
  # count number of replicates each protein found in
  count <- dplyr::group_by(data, {{ name }} , {{ sample }} )
  count <- dplyr::summarise(count,
                            n = dplyr::n(),
                            .groups = "drop")

  count <- dplyr::group_by(count, {{ name }} )
  count_ref <- dplyr::filter(count, {{ sample }} == reference)
  count_ref <- dplyr::select(count_ref, {{ name }}, ref_n = n)

  count <- dplyr::left_join(count, count_ref, by = dplyr::join_by( {{ name }} == {{ name }} ))

  # assign missingness for each protein
  count <- dplyr::mutate(count,
                         comparison = paste( {{ sample }} , reference, sep = "_vs_"),
                         missingness = dplyr::case_when(n == full & ref_n == full ~ "complete",
                                                        (n < MAR_thres & ref_n < MAR_thres) | (n != full & is.na(ref_n)) ~ NA,
                                                        (((n >= MAR_thres & ref_n >= MAR_thres) & (n + ref_n >= full)) | (n == full & ref_n > MNAR_thres) | (n > MNAR_thres & ref_n == full)) ~ "MAR",
                                                        ((n == full & ref_n <= MNAR_thres) | (n <= MNAR_thres & ref_n == full) | (n == full & is.na(ref_n))) ~ "MNAR",
                                                        TRUE ~ NA))

  # join missingness assingments into original data frame
  output <- dplyr::left_join(data, count, by = dplyr::join_by( {{ name }} == {{ name }} ,
                                                               {{ sample }} == {{ sample }} ))


  # reference data
  comparisons <- dplyr::filter(output, {{ sample }} != reference)
  comparisons <- dplyr::select(comparisons,
                               -c( {{ sample }}, {{ experiment }}, {{ replicate }}, {{ value }}, n))

  reference_data <- dplyr::filter(output, {{ sample }} == reference)
  reference_data <- dplyr::select(reference_data,
                                  {{ name }}, {{ sample }}, {{ experiment }}, {{ replicate }}, {{ value }}, n)

  reference_data <- dplyr::left_join(comparisons, reference_data, by = dplyr::join_by( {{ name }} == {{ name }} ), relationship = "many-to-many")
  reference_data <- dplyr::filter(reference_data, !is.na( {{ sample }} ))


  output <- dplyr::filter(output, {{ sample }} != reference)
  output <- dplyr::bind_rows(output, reference_data)

  # return output data frame
  output
}
