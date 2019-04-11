
#' Report the two clusters with the nearest centers
#'
#' This function reports the two clusters with the nearest centers.
#' @param x a numeric matrix. Columns are samples and rows are dimensions.
#' @param clustering a vector indicating the clustering label for each data point. 
#' @keywords IKAP
#' @export
#' @examples
#' NearestCluster(x, clustering)

NearestCluster <- function(x, clustering){
  clust.center <- data.frame(matrix(NA, ncol = length(unique(clustering)), nrow = ncol(x)))
  colnames(clust.center) <- unique(clustering)
  for(clust in unique(clustering)){
    clust.center[,clust] <- if(length(which(clustering == clust))==1) as.numeric(x[which(clustering == clust),]) else
      as.numeric(apply(x[which(clustering == clust),], 2, mean))
  }
  
  dist.mat <- as.matrix(dist(t(clust.center)))
  diag(dist.mat) <- Inf
  
  cluster.nearest <- colnames(clust.center)[which(upper.tri(dist.mat)*(dist.mat == min(dist.mat)) == 1, arr.ind = T)[1,]]
  
  return(cluster.nearest)
}