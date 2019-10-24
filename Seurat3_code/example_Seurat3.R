
# In this example, suppose your working directory is ./
# put IKAP_Seurat3.R in the working directory and raw data matrix (such as 10x data) under ./raw/

# read IKAP functions
source("IKAP_Seurat3.R")

# create Seurat object from raw data. Here I use 10x data as an example.
temp.obj <- Read10X(data.dir = "./raw/")
sobj <- CreateSeuratObject(counts = temp.obj, min.cells = 3, min.features = 200, project = "project")

# compute percentage of mitochondrial genes (optional)
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

# filter cells (optional)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & percent.mt < 5)

# normalize UMI counts
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

# run IKAP
sobj <- IKAP(sobj, out.dir = "./IKAP")

# save the Seurat object with IKAP results
saveRDS(sobj,"./IKAP/sobj.rds")


###### Example code for rerunning IKAP for a selected major group ######

sobj <- readRDS("./IKAP/sobj.rds")

# suppose you want to rerun IKAP for cluster 2 using clustering PC10K8
# subset the data to select cells in cluster 2 using clustering PC10K8
sobj.PC10K8C2 <- subset(sobj, subset = PC10K8 == "2")

# remove all clustering generated from the previous IKAP run in the Seurat metadata
sobj.PC10K8C2@meta.data <- sobj.PC10K8C2@meta.data[,-grep("^PC",colnames(sobj.PC10K8C2@meta.data))]

# remove variable genes identified in the previous IKAP run
VariableFeatures(sobj.PC10K8C2) <- c()

# run IKAP for the subsetted data
sobj.PC10K8C2 <- IKAP(sobj.PC10K8C2, out.dir = "./IKAP/PC10K8C2")

# save the Seurat object
saveRDS(sobj.PC10K8C2,"./IKAP/PC10K8C2/sobj.rds")

