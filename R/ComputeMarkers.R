
#' Compute upregulated marker genes for each cell group in the Seurat object
#'
#' This function computes the list of upregulated marker genes for each cell group in the Seurat object 
#' using Seurat::FindAllMarkers. It also computes the area under the ROC curve (AUROC) for each marker gene.
#' Finally it saves all relevant in Excel sheets.
#' @param sobj the Seurat object
#' @param gap.gain a numeric matrix for the increase of gap statistic from k-1 to k. Columns are
#' the number of clusters (k = 2 to k.max) and rows are the number of top principal components (nPC).
#' For example, the first column is the increase of gap statistic from k=1 to k=2 for each nPC (row).
#' @param candidates the list of candidates returned from SelectCandidate().
#' @param out.dir the path for output directory
#' @keywords IKAP
#' @export
#' @examples
#' ComputeMarkers(sobj, gap.gain, candidates, out.dir)

ComputeMarkers <- function(sobj, gap.gain, candidates, out.dir){
  
  markers.all <- list()
  out.xls <- list()
  out.xls$gap.gain <- cbind(PC_K = rownames(gap.gain), gap.gain)
  
  for(i in 1:length(candidates$k)){
    clustering.label <- paste0("PC",candidates$pc[i],"K",candidates$k[i])
    sobj <- SetAllIdent(sobj, id = clustering.label)
    sobj <- RunTSNE(sobj, dims.use = 1:candidates$pc[i])
    ggsave(TSNEPlot(sobj, do.return = T), filename = paste0(out.dir,"/",clustering.label,"_tSNE.pdf"))
    
    sobj.markers <- FindAllMarkers(object= sobj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
    
    sobj.markers$AUROC <- NA
    for(j in 1:nrow(sobj.markers)){
      sobj.markers$AUROC[j] <- roc.curve(scores.class0 = sobj@data[sobj.markers$gene[j],], weights.class0 = sobj@ident == sobj.markers$cluster[j])$auc
    }
    
    top.10 <- sobj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
    ggsave(DoHeatmap(object = sobj, genes.use = top.10$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row = 7),
           filename = paste0(out.dir,"/",clustering.label,"_DE_genes_LCF.png"), units = "in", width = 12, height = 8)
    
    out.xls[[clustering.label]] <- sobj.markers[,c("gene","p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","AUROC")]
    markers.all[[clustering.label]] <- sobj.markers
  }
  
  WriteXLS(out.xls, ExcelFileName = paste0(out.dir,"/data.xls"))
  saveRDS(markers.all, file = paste0(out.dir,"/markers.all.rds"))
  
  return(markers.all)
}
