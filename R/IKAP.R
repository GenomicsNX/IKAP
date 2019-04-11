
#' IKAP: Identifying K mAjor cell Population group in single cell analysis
#'
#' IKAP identifies candidate set(s) of major cell groups using the single cell analysis R package
#' Seurat by evaluating sets of possible cell groups generated using different parameters in Seurat
#' SNN clustering (i.e. resulotion r and the number of top principal components (nPC)). The results
#' (tables and plots) are saved in the output directory. A Seurat object is returned with all sets of
#' evaluated cell groups saved in the metadata data frame. 
#' @param sobj a Seurat object with cell expression normalized
#' @param pcs the list of principal components used for clustering. default is NA (to be determined by IKAP; recommended)
#' @param pc.range the range of nPCs. default is 20
#' @param k.max the maximal number of clusters. default is NA (to be determined by IKAP; recommended)
#' @param r.kmax.est resolution for IKAP determining k.max using Seurat SNN clustering. default is 1.5
#' @param out.dir the path for output directory
#' @param scale.data whether scale the data using Seurat::ScaleData. default is TRUE (recommended)
#' @param confounders a vector of confounders that need to be regressed out in Seurat::ScaleData.
#' default is c('nUMI','percent.mito') (see Seurat tutorial: https://satijalab.org/seurat/pbmc3k_tutorial.html)
#' @param plot.decision.tree whether to plot decision trees that classify the cell groups. default is TRUE
#' @param random.seed random seed
#' @keywords IKAP
#' @export
#' @examples
#' sobj.new <- IKAP(sobj, out.dir = "./IKAP")
#' 
#' saveRDS(sobj.new, file = "./IKAP/sobj.new.rds")

IKAP <- function(sobj, pcs = NA, pc.range = 20, k.max = NA, r.kmax.est = 1.5, out.dir = "./IKAP", scale.data = TRUE,
                 confounders = c('nUMI','percent.mito'), plot.decision.tree = TRUE, random.seed = 0){
  
  dir.create(out.dir, recursive = T)
  
  if(scale.data){
    if(!all(confounders %in% colnames(sobj@meta.data))){
      warning(confounders[which(!confounders %in% colnames(sobj@meta.data))],"not in Seurat metadata: skipped for regression.\n")
    }
    
    confounders <- intersect(confounders, colnames(sobj@meta.data))
    
    if(length(confounders) > 0) sobj <- ScaleData(sobj, vars.to.regress = confounders)
    else sobj <- ScaleData(sobj)
  }
  
  cat("Finding variable genes for clustering ... \n")
  sobj <- FindVariableGenes(sobj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
  
  cat("Running PCA ... \n")
  if(is.na(pcs)){
    sobj <- RunPCA(sobj, pcs.compute = 50, do.print = F)
    pc.change <- which(abs(diff(sobj@dr$pca@sdev)/sobj@dr$pca@sdev[2:length(sobj@dr$pca@sdev)]) > 0.1)
    while(length(pc.change) > 0 && max(pc.change)+pc.range+2 > length(sobj@dr$pca@sdev)){
      sobj <- RunPCA(sobj, pcs.compute = max(pc.change)+pc.range+2+10, do.print = F)
      pc.change <- which(abs(diff(sobj@dr$pca@sdev)/sobj@dr$pca@sdev[2:length(sobj@dr$pca@sdev)]) > 0.1)
    }
    pcs <- if(length(pc.change) == 0) 2:(pc.range+2) else (max(pc.change)+2):(max(pc.change)+pc.range+2)
  } else {
    sobj <- RunPCA(sobj, pcs.compute = max(pcs))
  }
  
  if(is.na(k.max)){
    cat("Determine k.max.\n")
    
    sobj <- FindClusters(object = sobj, reduction.type = "pca", dims.use = 1:min(pcs), resolution = r.kmax.est, print.output = 0, save.SNN = TRUE,
                         random.seed = random.seed)
    k.min.pc <- length(unique(sobj@ident))
    sobj <- FindClusters(object = sobj, reduction.type = "pca", dims.use = 1:max(pcs), resolution = r.kmax.est, print.output = 0, save.SNN = TRUE,
                         random.seed = random.seed)
    k.max.pc <- length(unique(sobj@ident))
    
    k.max <- as.integer((k.min.pc + k.max.pc)/2)
    
    cat("k.max =", k.max, "\n")
  }
  
  
  gap.gain <- data.frame(matrix(NA, ncol = k.max - 1, nrow = length(pcs)))
  colnames(gap.gain) <- as.character(2:k.max)
  rownames(gap.gain) <- paste0(pcs)
  
  cat("Perform clustering for every nPC:\n")
  for(npc in pcs){
    clusterings <- BottomUpMerge(sobj, k.max, npc, random.seed)
    gap.stat <- GapStatistic(sobj@dr$pca@cell.embeddings[,1:npc], clusterings)
    
    names(clusterings) <- paste0("PC",npc,"K",1:k.max)
    sobj@meta.data <- cbind(sobj@meta.data, as.data.frame(clusterings)[,2:k.max])
    gap.gain[as.character(npc),] <- diff(gap.stat$gap)
  }
  
  candidates <- SelectCandidate(gap.gain)
  
  cat("Compute marker gene lists ... \n")
  markers.all <- ComputeMarkers(sobj, gap.gain, candidates, out.dir)
  
  cat("Build decision tree ... \n")
  summary.rpart <- DecisionTree(sobj, markers.all, out.dir, plot.decision.tree)
  
  cat("Plotting summary ... \n")
  PlotSummary(gap.gain, summary.rpart, markers.all, out.dir)
  
  return(sobj)
}