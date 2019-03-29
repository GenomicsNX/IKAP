IKAP – "Identifying K mAjor cell Population groups in single-cell RNA-seq analysis"
==================================================================================

This R script is accompanying software for the paper: IKAP – "Identifying K mAjor cell Population groups in single-cell RNA-seq analysis", with authors Yun-Ching Chen, Abhilash Suresh, Chingiz Underbayev, Clare Sun, Komudi Singh, Fayaz Seifuddin, Adrian Wiestner, and Mehdi Pirooznia.

The main function, IKAP, takes a Seurat object with the normalized expression matrix and other parameters set by default values if not specified. IKAP explores sets of cell groups (clustering) by varying resolution (r) and the number of top principal components (nPC) for Seurat SNN clustering and picks a few candidate sets among all explored sets with one marked as the best that likely produces distinguishing marker genes.

Note: IKAP will, by default, regress out the percentage of mitochondrial gene counts and total UMI counts and scale the expression matrix using Seurat ScaleData function. These two values should be save in Seurat metadata with column names 'percent.mito' and 'nUMI' respectively. If you want to regress out different confounding variables or use different column names, please save these variables in Seurat metadata and set 'confounders' (an IKAP parameter) as their column names in the Seurat metadata data frame.

Please install the following R libraries for running IKAP:
Seurat, dplyr, reshape2, PRROC, WriteXLS, rpart, stringr, and rpart.plot


Usage:
-------

Seurat_obj <- IKAP(Seurat_obj, out.dir = "./IKAP")

Returned data and output files (saved in the output directory, default = ./IKAP/):

Seurat object:
IKAP returns a Seurat object with all explored sets in the metadata data frame.

- PC_K.pdf:
The heatmap shows the statistics for every combination of r and nPC explored. Candidate sets are marked as 'X' with the best marked as 'B'. The corresponding cell membership can be found in the metadata of the returned Seurat object with column name 'PC?K?'. For example, if 'B' (the best set) is marked at nPC = 20 and k = 8, the corresponding cell membership is stored in column 'PC20K8' in the metadata.

- data.xls and markers.all.rds:
It saves the statistics (plotted in PC_K.pdf) for determining candidate sets in the first sheet. The other sheets display the (upregulated) marker genes for candidate sets. The R object, markers.all.rds, contains a data frame of marker genes for every candidate set.

- *.png:
Heatmaps show expression of top 10 (ranked by expression fold change) marker genes from each cell group for candidate sets. They are plotted using Seurat DoHeatmap function.

DT_plot.pdf, DT_summary.rds, and DT.rds:
Decision tree output files. A decision tree is built using marker genes for every cell group in every candidate set using R package rpart. All decision trees are plotted in DT_plot.pdf. Classification errors are summarized in the R object DT_summary.rds. DT.rds is the output object from rpart.

- *_tSNE.pdf:
tSNE plots for candidate sets.




Functions in the R script:
--------------------------

- IKAP:
The main function runs the following steps: (1) regress out confounding variables and scale data using Seurat::ScaleData; (2) find variable genes for principal component analysis (PCA) using Seurat::FindVariableGenes; (3) perform PCA using Seurat::RunPCA; (4) estimate k.max; (5) explore ranges of k and nPC and compute gap statistics; (6) select candidate sets; (7) compute marker genes using Seurat::FindAllMarkers; (8) build decision trees; and (9) plot tSNE plots and PC_K.pdf.

- GapStatistic, ObservedLogW, and ExpectedLogW (5):
Compute gap statistics given a data matrix (used for computing data point Euclidean distances) and K sets of clusters with k = 1 … K. GapStatistic calls ObservedLogW and ExpectedLogW to compute sum of within-group distances for observed data and random data respectively.

- BottomUpMerge and NearestCluster (5):
Generate sets of cell groups by exploring ranges of k and nPC. BottomUpMerge finds k.max groups using Seurat::FindClusters and gradually merges two nearest clusters measured by NearestCluster.

- SelectCandidate (6):
Select candidate sets based on gap statistics.

- ComputeMarkers (7):
Compute marker genes for all cell groups in all candidate sets using Seurat::FindAllMarkers. In addition, compute Area Under the ROC curve (AUROC) for each marker genes using the R package PRROC. Plot marker gene heatmap(s) using Seurat::DoHeatmap.

- DecisionTree (8):
Build decision trees for all cell groups in all candidate sets using the R package rpart and compute the classification error for each candidate set.

- PlotSummary (9):
Mark the best set based on classification error and plot PC_K.pdf.




Contact
--------
If you have any question, please contact: yun-ching.chen@nih.gov

