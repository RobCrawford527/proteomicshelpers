pca_function <- function(data,
                         names_col,
                         samples,
                         sample_format,
                         colour_col = "combined",
                         shape_col = NA){

  # columns
  sample_format <- strsplit(sample_format,split="_")[[1]]

  # select value columns and remove missing values
  # print number of proteins
  rownames(data) <- data[,names_col]
  pca <- t(na.omit(data[,samples]))
  print(ncol(pca))

  # calculate principle components
  pr <- prcomp(pca,scale=TRUE,retx=TRUE,center=TRUE)

  # calculate percentage of variance explained by PC1 and PC2
  var_pc1 <- round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1)
  var_pc2 <- round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1)

  # add sample data to pr
  pr <- as.data.frame(pr$x)
  pr[,c("sample","sample2")] <- rownames(pr)
  pr <- separate(pr,col=sample2,into=sample_format,sep="_")

  # new columns for colour and shape
  # if option is "combined", unite columns to form "condition" column
  if (is.na(colour_col)){
    pr[,"colour_col"] <- "colour"
  } else if (colour_col == "combined"){
    pr <- unite(pr,col="colour_col",sample_format[sample_format!="Replicate"],sep="_")
  } else {
    pr[,"colour_col"] <- pr[,colour_col]
  }
  if (is.na(shape_col)){
    pr[,"shape_col"] <- "shape"
  } else if (colour_col == "combined"){
    pr[,"shape_col"] <- "shape"
  } else {
    pr[,"shape_col"] <- pr[,shape_col]
  }

  # plot PCA result
  plot <- ggplot(pr,aes(x=PC1,y=PC2,colour=colour_col,shape=shape_col)) +
    geom_point(size=5,alpha=0.8) +
    geom_text(aes(label=Replicate),colour="black") +
    coord_equal() +
    xlab(paste("PC1 (",var_pc1,"% of variance)",sep="")) +
    ylab(paste("PC2 (",var_pc2,"% of variance)",sep="")) +
    theme_classic() +
    theme(legend.position='right')

  # return plot
  plot
}
