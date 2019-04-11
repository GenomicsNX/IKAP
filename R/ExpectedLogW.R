
#' Compute log expected sum of pairwise distances for a given grouping
#'
#' This function computes log expected sum of pairwise distances for a given grouping. It assumes data points are
#' uniformly distributed within a range of values in every dimension. It first computes the expected sum of
#' pairwise distances for each group (cluster) and then takes log for the sum of expected sums acros all groups.
#' @param pc.ranges a numeric matrix: the minimum (1st row) and the maximum (2nd row) for each dimension (column)
#' @param one.clustering a vector indicating the clustering membership for each data point
#' @keywords IKAP
#' @export
#' @examples
#' ExpectedLogW(pc.ranges, one.clustering)

ExpectedLogW <- function(pc.ranges, one.clustering){
  dist.expected <- sqrt(sum(apply(pc.ranges, 2, function(r){(r[1]-r[2])^2/6})))
  dist.sum.expected <- 0.5*sum(sapply(table(one.clustering), function(n){dist.expected*(n-1)}))
  return(log(dist.sum.expected))
}