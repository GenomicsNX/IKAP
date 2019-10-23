
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