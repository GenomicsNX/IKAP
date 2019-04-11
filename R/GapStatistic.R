
#' Compute gap statistic
#'
#' This function computes gap statistics for K clusterings with the number of clusters = 1 to K. 
#' @param x a numeric matrix: Columns are samples and rows are dimensions.
#' @param clusterings a list of vectors, each of which indicates the clustering label for each data point
#' @keywords IKAP
#' @export
#' @examples
#' GapStatistic(x, clusterings)

GapStatistic <- function(x, clusterings){
  gap.stat <- data.frame(matrix(NA, ncol = 3, nrow = length(clusterings)))
  colnames(gap.stat) <- c("E.log.W","log.W","gap")
  
  for(i in 1:length(clusterings)){
    gap.stat[i,"E.log.W"] <- ExpectedLogW(apply(x, 2, range), clusterings[[i]])
    gap.stat[i,"log.W"] <- ObservedLogW(x, clusterings[[i]])
  }
  
  gap.stat[,"gap"] <- gap.stat[,"E.log.W"] - gap.stat[,"log.W"]
  
  return(gap.stat)
}
