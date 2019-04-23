IKAP – Identifying K mAjor cell Population groups in single-cell RNA-seq analysis
------------------------------

<br>

> <div>BioRxiv preprint:<br><b>IKAP - Identifying K mAjor cell Population groups in single-cell RNA-seq analysis</b>
>  <br>Yun-Ching Chen, Abhilash Suresh, Chingiz Underbayev, Clare Sun, Komudi Singh, Fayaz Seifuddin, Adrian Wiestner, Mehdi Pirooznia. 
  >  <b>doi:</b> <a href="https://doi.org/10.1101/596817" class="" classname="" target="_blank" name="">https://doi.org/10.1101/596817</a></div> 

<br>


Installation
-----------

Please install the following R libraries before installing IKAP: <br>
[Seurat](https://satijalab.org/seurat/install.html), [dplyr](https://cloud.r-project.org/web/packages/dplyr), [reshape2](https://cran.r-project.org/web/packages/reshape2), [PRROC](https://cran.r-project.org/web/packages/PRROC/), [WriteXLS](https://cran.r-project.org/web/packages/WriteXLS/), [rpart](https://cran.r-project.org/web/packages/rpart/), [stringr](https://cran.r-project.org/web/packages/stringr), and [rpart.plot](https://cran.r-project.org/web/packages/rpart.plot) 

<br>
<h4>IKAP installation:</h4>

<ol>
  <li>
    <p>First, you need to install the
<a href="https://github.com/hadley/devtools">devtools</a> package. You can do
this from <a href="https://cran.r-project.org">CRAN</a>. Invoke R and then type</p>

```{r, eval = FALSE}
install.packages("devtools")
```
  </li>
  <li>
    <p>Load the devtools package.</p>

```{r, eval = FALSE}
library(devtools)
```
  </li>
  <li>
    <p>Install IKAP </p>

```{r, eval = FALSE}
devtools::install_github("NHLBI-BCB/IKAP")
```
  </li>
</ol>








<br>
The main function, IKAP, takes a Seurat object with the normalized expression matrix and other parameters set by default values if not specified. IKAP explores sets of cell groups (clustering) by varying resolution (r) and the number of top principal components (nPC) for Seurat SNN clustering and picks a few candidate sets among all explored sets with one marked as the best that likely produces distinguishing marker genes.

Note: IKAP will, by default, regress out the percentage of mitochondrial gene counts and total UMI counts and scale the expression matrix using Seurat ScaleData function. These two values should be save in Seurat metadata with column names 'percent.mito' and 'nUMI' respectively. If you want to regress out different confounding variables or use different column names, please save these variables in Seurat metadata and set 'confounders' (an IKAP parameter) as their column names in the Seurat metadata data frame.


<div class="paragraph"><br><br></div>

IKAP Workflow 
--------
<p align="center">
  <img src="Figure 1.png" width="933" height="1050" title="IKAP workflow">
</p>


Usage:
-------

```
Seurat_obj <- IKAP(Seurat_obj, out.dir = "./IKAP")
```

Returned data and output files (saved in the output directory, default = ./IKAP/):

Seurat object:
IKAP returns a Seurat object with all explored sets in the metadata data frame.

- **_PC_K.pdf_**:

The heatmap shows the statistics for every combination of r and nPC explored. Candidate sets are marked as 'X' with the best marked as 'B'. The corresponding cell membership can be found in the metadata of the returned Seurat object with column name 'PC?K?'. For example, if 'B' (the best set) is marked at nPC = 20 and k = 8, the corresponding cell membership is stored in column 'PC20K8' in the metadata.

- **_data.xls_** and **_markers.all.rds_**:

It saves the statistics (plotted in PC_K.pdf) for determining candidate sets in the first sheet. The other sheets display the (upregulated) marker genes for candidate sets. The R object, markers.all.rds, contains a data frame of marker genes for every candidate set.

- **_*.png_**:

Heatmaps show expression of top 10 (ranked by expression fold change) marker genes from each cell group for candidate sets. They are plotted using Seurat DoHeatmap function.

- **_DT_plot.pdf_**, **_DT_summary.rds_**, and **_DT.rds_**:

Decision tree output files. A decision tree is built using marker genes for every cell group in every candidate set using R package rpart. All decision trees are plotted in DT_plot.pdf. Classification errors are summarized in the R object DT_summary.rds. DT.rds is the output object from rpart.

- **_*tSNE.pdf_**:

tSNE plots for candidate sets.

<div class="paragraph"><br><br></div>



Functions in the R script:
--------------------------

- IKAP:
The main function runs the following steps: 
    - (1) regress out confounding variables and scale data using Seurat::ScaleData; 
    - (2) find variable genes for principal component analysis (PCA) using Seurat::FindVariableGenes; 
    - (3) perform PCA using Seurat::RunPCA; 
    - (4) estimate k.max; 
    - (5) explore ranges of k and nPC and compute gap statistics; 
        - GapStatistic, ObservedLogW, and ExpectedLogW: <br>Compute gap statistics given a data matrix (used for computing data point Euclidean distances) and K sets of clusters with k = 1 … K. GapStatistic calls ObservedLogW and ExpectedLogW to compute sum of within-group distances for observed data and random data respectively.
        - BottomUpMerge and NearestCluster (5): Generate sets of cell groups by exploring ranges of k and nPC. BottomUpMerge finds k.max groups using Seurat::FindClusters and gradually merges two nearest clusters measured by NearestCluster.
    - (6) select candidate sets; 
        - SelectCandidate: <br>Select candidate sets based on gap statistics.
    - (7) compute marker genes using Seurat::FindAllMarkers; 
        - ComputeMarkers: <br>Compute marker genes for all cell groups in all candidate sets using Seurat::FindAllMarkers. In addition, compute Area Under the ROC curve (AUROC) for each marker genes using the R package PRROC. Plot marker gene heatmap(s) using Seurat::DoHeatmap.
    - (8) build decision trees; 
        - DecisionTree: <br>Build decision trees for all cell groups in all candidate sets using the R package rpart and compute the classification error for each candidate set.
    - (9) plot tSNE plots and PC_K.pdf 
        - PlotSummary: <br>Mark the best set based on classification error and plot PC_K.pdf

 

<div class="paragraph"><br><br></div>


 
 
 
 
License
--------
MIT license: https://opensource.org/licenses/MIT 
 
<br><br> 
Contact
--------
If you have any question, please contact: yun-ching.chen@nih.gov

