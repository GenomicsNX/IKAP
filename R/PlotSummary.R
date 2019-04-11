
#' Plot the IKAP result
#'
#' This function plots a heat map indicating the increase of gap statistics as the number of cluster (k)
#' increases by varying the number of top principal components (nPC) used. The best set and alternative
#' sets of cell groups are indicated by 'B' and 'X' respectively. The best set is determined based on
#' classification error summarized in DecisionTree().
#' @param gap.gain a numeric matrix for the increase of gap statistic from k-1 to k. Columns are
#' the number of clusters (k = 2 to k.max) and rows are the number of top principal components (nPC).
#' For example, the first column is the increase of gap statistic from k=1 to k=2 for each nPC (row).
#' @param summary.rpart a numberic matrix with overall classification error summarized for each n split (column)
#' and each candidate set (row). This is usually the output return from DecisionTree().
#' @param markers.all a list of marker gene tables for candidate sets. This is usually the list returned from
#' ComputeMarkers().
#' @param out.dir the path for output directory
#' @keywords IKAP
#' @export
#' @examples
#' PlotSummary(gap.gain, summary.rpart, markers.all, out.dir)

PlotSummary <- function(gap.gain, summary.rpart, markers.all, out.dir){
  
  data.plot <- melt(as.matrix(gap.gain))
  colnames(data.plot) <- c("PC","k","gap.gain")
  data.plot$PC <- factor(data.plot$PC, levels = rev(rownames(gap.gain)))
  data.plot$k <- factor(data.plot$k, levels = colnames(gap.gain))
  
  temp <- str_split_fixed(gsub("^PC","",rownames(summary.rpart)),"K",2)
  best.ind <- if(nrow(summary.rpart) > 1) which.min(apply(summary.rpart[,5:15], 1, mean)) else 1
  
  data.plot$label <- ""
  for(i in 1:nrow(temp)) data.plot$label[which(data.plot$PC == temp[i,1] & data.plot$k == temp[i,2])] <- "X"
  data.plot$label[which(data.plot$PC == temp[best.ind,1] & data.plot$k == temp[best.ind,2])] <- "B"
  
  p <- ggplot(data = data.plot, aes(x = k, y=PC, fill = gap.gain, label = label)) + geom_tile() + scale_x_discrete(position = "top") + geom_text() +
    scale_fill_gradient(low = "darkgrey", high = "darkred") + labs(x="# of clusters (k)",y="# of top PCs")
  
  ggsave(p, filename = paste0(out.dir,"/PC_K.pdf"))
}