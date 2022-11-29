enrichGO_plot <- function(go_result, save=FALSE, filename, width=160){

  # find min and max p.adjust values
  # round to whole number
  min_p <- round(min(go_result$p.adjust)-0.5,0)
  max_p <- round(max(go_result$p.adjust)+0.5,0)
  # adjust min_p if too low
  if (min_p < 1){
    min_p <- 1
  }
  # define limits and breaks
  limits <- c(1,max_p)
  breaks <- unique(round(seq(min_p+1,max_p-1,(max_p-min_p-2)/3),0))

  # plot
  go_plot <- ggplot(go_result,aes(x=GeneRatio,y=Description,colour=p.adjust,size=p.adjust,group=ID)) +
    geom_point() +
    scale_colour_viridis(option='plasma',discrete=FALSE,begin=0,end=0.85,direction=-1,
                         limits=limits,breaks=breaks) +
    scale_size(limits=limits,breaks=breaks) +
    scale_x_continuous(limits=c(0,max(go_result[,"GeneRatio"])+0.05),
                       breaks=seq(0,1,0.1),expand=c(0,0)) +
    scale_y_discrete(expand=c(0,1)) +
    facet_grid(ONTOLOGY ~ .,scales='free_y',space='free_y') +
    xlab(label="Gene ratio") +
    theme_classic() +
    theme(axis.title.y=element_blank(),axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black",size=8),strip.background=element_blank(),
          strip.text=element_blank(),legend.position='right')

  # save if appropriate
  if (save == TRUE){
    ggsave(filename=filename,plot=go_plot,device='svg',
           height=nrow(go_result)*3.5+20,width=width,units='mm')
  }

  # return plot
  go_plot
}
