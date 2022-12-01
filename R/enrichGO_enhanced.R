#' Perform Gene Ontology Analysis
#'
#' @param genes The vector of gene names to be analysed
#' @param OrgDb Organism database to use
#' @param keyType Starting format of the names (default UNIPROT), as
#'     defined by clusterProfiler
#' @param ont The ontologies to analyse (BP, CC and/or MF; defaults to "all")
#' @param reference Data frame containing the background list in multiple
#'     formats
#' @param simplify Whether to simplify the GO result
#' @param evaluate Whether to evaluate the GO result
#' @param convert Whether to convert the names in the GO result
#' @param toType Format to convert the names to (if convert = TRUE)
#'
#' @return A data frame containing the enriched GO terms
#' @export
#'
#' @examples
#'
enrichGO_enhanced <- function(genes,
                              OrgDb,
                              keyType = "ORF",
                              ont = "all",
                              reference,
                              simplify = TRUE,
                              evaluate = TRUE,
                              convert = TRUE,
                              toType = "GENENAME"){

  # GO enrichment analysis on gene list
  go <- clusterProfiler::enrichGO(gene = genes,
                                  OrgDb = OrgDb,
                                  keyType = keyType,
                                  ont = ont,
                                  universe = reference[,keyType])

  # simplify result
  if (simplify == TRUE){
    go <- clusterProfiler::simplify(go)
  }

  # extract result data frame
  go <- go@result

  # evaluate columns to convert to numeric
  if (evaluate == TRUE){
    go$GeneRatio <- apply(as.matrix(go$GeneRatio),
                          1,
                          function(x) eval(parse(text = x)))
    go$BgRatio <- apply(as.matrix(go$BgRatio),
                        1,
                        function(x) eval(parse(text = x)))
    go$p.adjust <- -log10(go$p.adjust)
    go$Description <- factor(go$Description,
                             levels = go$Description[order(go$GeneRatio,
                                                           decreasing = FALSE)])
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
