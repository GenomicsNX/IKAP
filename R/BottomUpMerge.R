
#' Generate clusterings by iteratively merging two nearest clusters for a Seurat object
#'
#' This function generates clusterings by iteratively merging two nearest clusters for a Seurat object. It generates an
#' initial clustering with the number of clusters (k) >= k.max using Seurat SNN clustering. Then iteratively merge
#' two clusters with the nearest centers to generate the clustering with k-1 clusters. Finally, it returns k.max
#' clusterings with k = 1 to k.max.
#' @param sobj Seurat object
#' @param k.max the maximal number of clusters 
#' @param npc the number of top principal components used for Seurat SNN clustering
#' @param random.seed the random seed set for Seurat SNN clustering
#' @keywords IKAP
#' @export
#' @examples
#' BottomUpMerge(sobj, k.max, npc, random.seed)

BottomUpMerge <- function(sobj, k.max, npc, random.seed){
  k.clustering <- 0
  clusterings <- list()
  
  clust.r <- 1.0
  sobj <- FindClusters(object = sobj, reduction.type = "pca", dims.use = 1:npc, resolution = clust.r, print.output = 0, save.SNN = TRUE,
                       random.seed = random.seed)
  
  cat("Iteration for nPC =",npc,", r = 1.0")
  while(length(unique(sobj@ident)) < k.max){
    clust.r <- clust.r + 0.2
    cat(",", clust.r)
    sobj <- FindClusters(object = sobj, reduction.type = "pca", dims.use = 1:npc, resolution = clust.r, print.output = 0, save.SNN = TRUE,
                         random.seed = random.seed)
  }
  cat("\n")
  clusterings[[(k.clustering <- length(unique(sobj@ident)))]] <- as.character(sobj@ident)
  
  # Merging clusters by nearest centers
  while(k.clustering > 2){
    merged <- NearestCluster(sobj@dr$pca@cell.embeddings[,1:npc], clusterings[[k.clustering]])
    clustering.merged <- clusterings[[k.clustering]]
    clustering.merged[which(clustering.merged == merged[1])] <- merged[2]
    clusterings[[k.clustering-1]] <- as.character(as.integer(as.factor(clustering.merged)))
    
    k.clustering <- k.clustering - 1
  }
  
  clusterings[[1]] <- rep("1", nrow(sobj@dr$pca@cell.embeddings))
  
  return(clusterings[1:k.max])
}