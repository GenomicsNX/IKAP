
#' Build decision trees for each cell group in the candidate sets
#'
#' This function builds decision trees for each cell group each candidate set using the R package rpart.
#' The overall classification error is summarized for each candidate set of cell groups.
#' @param sobj the Seurat object
#' @param markers a list of marker gene tables for candidate sets. This is usually the list returned from
#' ComputeMarkers().
#' @param out.dir the path for output directory
#' @param plot.decision.tree whether to plot decision trees?
#' @keywords IKAP
#' @export
#' @examples
#' DecisionTree(sobj, markers, out.dir, plot.decision.tree)


DecisionTree <- function(sobj, markers, out.dir, plot.decision.tree){
  data.rpart <- list()
  
  for(candidate in names(markers)){
    data.rpart[[candidate]] <- list()
    genes.candidate <- unique(markers[[candidate]]$gene[which(markers[[candidate]]$p_val_adj < 0.01)])
    for(clust in as.character(unique(markers[[candidate]]$cluster))){
      if(length(genes.candidate) == 0){
        next
      } else if(length(genes.candidate) == 1){
        data <- data.frame(as.factor(sobj@meta.data[[candidate]] == clust), as.numeric(sobj@data[genes.candidate,]))
        colnames(data) <- c("label", genes.candidate)
      } else {
        data <- as.data.frame(t(as.matrix(sobj@data[genes.candidate,])))
        data$label <- as.factor(sobj@meta.data[[candidate]] == clust)
      }
      data.rpart[[candidate]][[clust]] <- rpart(label ~., data = data)
    }
  }
  
  summary.rpart <- data.frame(matrix(NA, nrow = length(markers), ncol = 20))
  rownames(summary.rpart) <- names(markers)
  colnames(summary.rpart) <- paste0("s",1:20)
  
  for(candidate in names(data.rpart)){
    for(nsplit in 1:20){
      err.rate <- 0
      for(clust in unique(as.character(sobj@meta.data[[candidate]]))){
        if(is.null(data.rpart[[candidate]][[clust]])){
          err.rate <- err.rate + 1.0
        } else {
          err.rate <- err.rate + 
            data.rpart[[candidate]][[clust]]$cptable[max(which(data.rpart[[candidate]][[clust]]$cptable[,"nsplit"] <= nsplit)),"rel error"]*
            length(which(sobj@meta.data[[candidate]] == clust))/nrow(sobj@meta.data)
        }
      }
      summary.rpart[candidate, nsplit] <- err.rate
    }
  }
  
  # plot decision tree
  pdf(paste0(out.dir,"/DT_plot.pdf"))
  if(plot.decision.tree){
    for(candidate in names(data.rpart)){
      clust.sorted <- unique(as.character(sobj@meta.data[[candidate]]))
      if(!any(is.na(suppressWarnings(as.integer(clust.sorted))))) clust.sorted <- as.character(sort(as.integer(clust.sorted)))
      for(clust in clust.sorted){
        if(!is.null(data.rpart[[candidate]][[clust]])){
          rpart.plot(data.rpart[[candidate]][[clust]], roundint = F, main = paste0(candidate,":",clust))
        }
      }
    }
  }
  dev.off()
  
  saveRDS(data.rpart, file = paste0(out.dir,"/DT.rds"))
  saveRDS(summary.rpart, file = paste0(out.dir,"/DT_summary.rds"))
  
  return(summary.rpart)
}
