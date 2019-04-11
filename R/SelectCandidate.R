
#' Select the candidate sets of cell groups
#'
#' This function selects the candidate sets of cell groups based on the changes of gap statistics.
#' @param gap.gain a numeric matrix for the increase of gap statistic from k-1 to k. Columns are
#' the number of clusters (k = 2 to k.max) and rows are the number of top principal components (nPC).
#' For example, the first column is the increase of gap statistic from k=1 to k=2 for each nPC (row).
#' @keywords IKAP
#' @export
#' @examples
#' SelectCandidate(gap.gain)

SelectCandidate <- function(gap.gain){
  gg.mean <- mean(as.matrix(gap.gain))
  gg.sd <- sd(as.matrix(gap.gain))
  
  gap.gain.candidates <- sort(apply(gap.gain*(gap.gain > gg.mean + gg.sd), 2, max), decreasing = T)
  gap.gain.candidates <- gap.gain.candidates[which(gap.gain.candidates > 0)]
  
  k.candidates <- c()
  pc.candidates <- c()
  
  if(length(gap.gain.candidates) == 0){
    max.ind <- which(gap.gain == max(gap.gain), arr.ind = T)
    k.candidates <- c(max.ind[1,2]+1)
    pc.candidates <- c(as.integer(rownames(gap.gain)[max.ind[1,1]]))
  } else {
    for(i in gap.gain.candidates){
      ind <- which(gap.gain == i, arr.ind = T)
      if(length(k.candidates) == 0){
        k.candidates <- c(ind[1,2]+1)
        pc.candidates <- c(as.integer(rownames(gap.gain)[ind[1,1]]))
      } else {
        k.temp <- ind[1,2]+1
        pc.temp <- as.integer(rownames(gap.gain)[ind[1,1]])
        
        if(all(k.temp > k.candidates) && all(pc.temp > pc.candidates)){
          k.candidates <- c(k.candidates, k.temp)
          pc.candidates <- c(pc.candidates, pc.temp)
        }
      }
    }
  }
  
  return(list(k = k.candidates, pc = pc.candidates))
}