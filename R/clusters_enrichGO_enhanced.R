clusters_enrichGO_enhanced <- function(data,
                                       OrgDb,
                                       keyType = "ORF",
                                       names_col = keyType,
                                       clus_col,
                                       ont = "all",
                                       reference,
                                       simplify = TRUE,
                                       evaluate = TRUE,
                                       convert = TRUE,
                                       toType = "GENENAME"){
  
  # change column names
  colnames(data)[colnames(data) == names_col] <- "name"
  colnames(data)[colnames(data) == clus_col] <- "cluster"
  
  # calculate GO enrichment for current category
  # simplify and format result
  go <- clusterProfiler::compareCluster(name ~ cluster,
                                        fun = "enrichGO",
                                        data = data,
                                        ont = ont,
                                        OrgDb = OrgDb,
                                        keyType = keyType,
                                        universe = reference[,keyType],
                                        pAdjustMethod = "BH",
                                        pool = FALSE,
                                        pvalueCutoff = 0.05)
  
  # simplify result
  if (simplify == TRUE){
    go <- clusterProfiler::simplify(go)
  }
  
  # extract result data frame
  go <- go@compareClusterResult
  
  # evaluate columns to convert to numeric
  if (evaluate == TRUE){
    go$GeneRatio <- apply(as.matrix(go$GeneRatio),
                          1,
                          function(x) eval(parse(text = x)))
    go$BgRatio <- apply(as.matrix(go$BgRatio),
                        1,
                        function(x) eval(parse(text = x)))
    go$p.adjust <- -log10(go$p.adjust)
    go$Description <- as.factor(go$Description)
  }
  
  # convert names to different format
  if (convert == TRUE){
    # create new names column
    go[,"new_names"] <- go[,"geneID"]
    for (r in 1:nrow(go)){
      # select names from row of interest
      # split string into individual names
      names <- go[r, "new_names"]
      names <- strsplit(names, split="/")[[1]]
      
      # convert each name
      for (p in 1:length(names)){
        names[p] <- reference[reference[,keyType] == names[p], toType]
      }
      
      # paste into original format
      # replace column with new names
      names <- paste(names, collapse="/")
      go[r, "new_names"] <- names
    }
  }
  
  # return output data frame
  go
}
