names_conversion <- function(data,
                             names_col,
                             fromType = "UNIPROT",
                             toType = c("GENENAME", "ORF"),
                             OrgDb = org.Sc.sgd.db::org.Sc.sgd.db){

  # extract names column from data frame
  names <- as.data.frame(data[,names_col])
  colnames(names) <- names_col

  # convert Uniprot IDs to other name types
  translated_names <- clusterProfiler::bitr(names[,names_col],
                                            fromType = fromType,
                                            toType = toType,
                                            OrgDb = OrgDb,
                                            drop = FALSE)

  # identify and flag duplicated names
  # only keep one version
  translated_names[,"Unique"] <- !(duplicated(translated_names[,fromType], fromLast = FALSE) | duplicated(translated_names[,fromType], fromLast = TRUE))
  translated_names <- translated_names[!duplicated(translated_names[,fromType]), c(fromType, toType, "Unique")]

  # replace missing names with original IDs
  for (c in toType){
    for (r in 1:nrow(translated_names)){
      if (is.na(translated_names[r, c])){
        translated_names[r, c] <- translated_names[r, fromType]
      }
    }
  }

  # return translated names data frame
  translated_names
}
