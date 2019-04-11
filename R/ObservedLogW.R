
#' Compute log observed sum of pairwise distances for a given grouping
#'
#' This function computes log observed sum of pairwise distances for a given grouping. It first computes the sum of
#' pairwise distances for each group (cluster) and then takes log for the sum of sums acros all groups.
#' @param X a numeric matrix. Columns are samples and rows are dimensions.
#' @param one.clustering a vector indicating the clustering membership for each data point
#' @keywords IKAP
#' @export
#' @examples
#' ObservedLogW(pc.ranges, one.clustering)

ObservedLogW <- function(X, one.clustering){
  if(is.data.frame(X)) X <- as.matrix(X)
  dist.sum.observed <- log(0.5 * sum(vapply(split(seq_len(nrow(X)), one.clustering), function(I) {
    xs <- X[I, , drop = FALSE]
    sum(dist(xs)/nrow(xs))
  }, 0)))
  
  return(dist.sum.observed)
}