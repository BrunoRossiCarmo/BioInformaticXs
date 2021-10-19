#scRNA-seq - packages
# https://satijalab.org/seurat/articles/get_started.html

# Enter commands in R (or R studio, if installed)
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('Seurat')

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
#https://satijalab.org/seurat/vignettes.html


?Read10X
# Inform the directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X to create the data.
data_10x <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
head(data_10x)
data_10x
# Lets examine a few genes in the first thirty cells
data_10x[c("CD3D","TCL1A","MS4A1"), 1:30]
dim(data_10x)
#The `.` values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0,  Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.
head(as.matrix(data_10x[1:3,1:3]))
dense.size <- object.size(as.matrix(data_10x))
dense.size


head(data_10x[1:3,1:3])
sparse.size <- object.size(data_10x)
sparse.size

dense.size / sparse.size


?CreateSeuratObject
seurat_obj <- CreateSeuratObject(counts = data_10x,
                                 project = "pbmc3k",
                                 min.cells = 3, # genes expresso em no min 3 cells
                                 min.features = 200) # manter cells com no min 200 genes
dim(seurat_obj)
head(seurat_obj[["RNA"]]@counts)
head(seurat_obj@meta.data)
head(seurat_obj$orig.ident)
head(seurat_obj$nCount_RNA)
head(seurat_obj$nFeature_RNA)

## Standard pre-processing workflow
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/

#QC and selecting cells for further analysis
#Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include:
# - The number of unique genes detected in each cell.
# - Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# - The percentage of reads that map to the mitochondrial genome

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

#In order to detect mitochondrial genes, we need to tell Seurat how to distinguish these genes. We do this using a regular expression as in “mito.genes <- grep(pattern = "^MT-". If your mitochondrial genes are named differently, then you will need to adjust this pattern accordingly (e.g. "mt-", “mt.”, or “MT_” , "MTN".)
grep("^MT-", rownames(seurat_obj), v=T)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
head(seurat_obj$percent.mt)
head(seurat_obj[["percent.mt"]])
dim(seurat_obj)

# Show QC metrics for the first 5 cells
head(seurat_obj@meta.data, 5)
seurat_obj@meta.data[1,2] # Cell 1 and nCount
seurat_obj@meta.data[1,3] # Cell 1 and nFeature
seurat_obj@meta.data[1,4] # mitochondrial percentage
seurat_obj@meta.data[1,2] / seurat_obj@meta.data[1,3] #counts per cell

grep("^MT-", rownames(seurat_obj), v=T)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_obj), value = TRUE)
percent.mito <- (Matrix::colSums(seurat_obj[mito.genes, ])/Matrix::colSums(seurat_obj))*100
head(percent.mito)

# Visualize QC metrics as a violin plot
?VlnPlot
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
?FetchData

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
?FeatureScatter
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#We filter out cells that have unique gene counts over 2,500 or less than 200 and percent.mt less than 5
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & nCount_RNA < 10000)
# & nCount_RNA < 10000
dim(seurat_obj)

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = "blue")

plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                        feature2 = "percent.mt",
                        cols = "Blue")
plot4 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA",
                        cols = "Blue")
plot3 + plot4

## Normalizing the data
#After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data.

?NormalizeData
head(seurat_obj[["RNA"]]@data)
head(colSums(seurat_obj))
summary(colSums(seurat_obj))

seurat_obj <- NormalizeData(seurat_obj,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)
head(seurat_obj[["RNA"]]@data)
head(as.matrix(seurat_obj[["RNA"]]@data))[1:3, 1:3]

head(colSums(seurat_obj))
summary(colSums(seurat_obj))

## Identification of highly variable features (feature selection)
#We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
#https://www.nature.com/articles/nmeth.2645

?FindVariableFeatures
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
top10
# plot variable features with and without labels
?VariableFeaturePlot
plot5 <- VariableFeaturePlot(seurat_obj)
plot5
?LabelPoints
plot6 <- LabelPoints(plot = plot5, points = top10, repel = TRUE)
plot6
seurat_obj

## Scaling the data
#Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The `ScaleData` function:
  
#  * Shifts the expression of each gene, so that the mean expression across cells is 0
# * Scales the expression of each gene, so that the variance across cells is 1
#   + This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# * The results of this are stored in `pbmc[["RNA"]]@scale.data`
all.genes <- rownames(seurat_obj)
?ScaleData
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
head(seurat_obj[["RNA"]]@scale.data[1:3,1:3])

#leandro- cod
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("percent.mt"), features = all.genes)

#How can I remove unwanted sources of variation
# we could ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination.
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt") #pay attention!!!

## Perform linear dimensional reduction
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
?RunPCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 10)

?VizDimLoadings
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")

?DimPlot
DimPlot(seurat_obj, reduction = "pca")

?DimHeatmap
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset

#To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset. However, how many componenets should we choose to include? 10? 20? 100?
  
#In Macosko et al, we implemented a resampling test inspired by the JackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
?JackStraw
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100) # take a long time!!!!

?ScoreJackStraw
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)

#The JackStrawPlot function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.
JackStrawPlot(seurat_obj, dims = 1:15)

#An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs
?ElbowPlot
ElbowPlot(seurat_obj)

# ElbowPlot - Harvard
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

#Cluster the cells
?FindNeighbors
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)

?FindClusters
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) #CHANGE resolution!!!
seurat_obj$RNA_snn_res.0.5 #information about the cluster's number that cell belongs
seurat_obj@meta.data
table(seurat_obj@meta.data$RNA_snn_res.0.5 == seurat_obj@meta.data$seurat_clusters)
# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
#Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
?RunUMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
seurat_obj$seurat_clusters

?RunTSNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
?DimPlot
DimPlot(seurat_obj, reduction = "umap")
DimPlot(seurat_obj, reduction = "tsne")

#You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(seurat_obj, file = "pbmc_tutorial.rds")

## Finding differentially expressed features (cluster biomarkers)
#Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

#The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed features will likely still rise to the top.

# find all markers of cluster 1
?FindMarkers
Idents(seurat_obj)
cluster1.markers <- FindMarkers(seurat_obj, ident.1 = 1, min.pct = 0.25, only.pos = F)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_obj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(seurat_obj, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
cluster1.markers_roc <- FindMarkers(seurat_obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#We include several tools for visualizing marker expression. VlnPlot (shows expression probability distributions across clusters), and FeaturePlot (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring RidgePlot, CellScatter, and DotPlot as additional methods to view your dataset.
?VlnPlot
plot7 <- VlnPlot(seurat_obj, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
plot8 <- VlnPlot(seurat_obj, features = c("MS4A1", "CD79A"), slot = "counts", log = TRUE)
plot7 + plot8
VlnPlot(seurat_obj, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

?FeaturePlot
FeaturePlot(seurat_obj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
FeaturePlot(seurat_obj,
            features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),
            reduction = "tsne")

#DoHeatmap generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
?DoHeatmap
DoHeatmap(seurat_obj, features = top10$gene, size = 2.5) + NoLegend()

?RidgePlot
RidgePlot(object = seurat_obj, features = 'PC_1')
?CellScatter
CellScatter(object = seurat_obj, cell1 = 'AAACATACAACCAC-1', cell2 = 'AAACATTGAGCTAC-1')
head(colnames(seurat_obj))
?DotPlot()
cd_genes <- c("CD247", "CD3E", "CD9")
DotPlot(object = seurat_obj, features = cd_genes)

## Assigning cell type identity to clusters
# Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat_obj)
?RenameIdents
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = 0.5) 

FeaturePlot(seurat_obj, features = c("IL7R","CCR7","S100A4","CD14","LYZ","MS4A1","CD8A","FCGR3A", "MS4A7","GNLY","NKG7", "FCER1A", "CST3", "PPBP"))

VlnPlot(seurat_obj, features = c("IL7R", "CCR7"), log = TRUE) #Naive CD4+ T
VlnPlot(seurat_obj, features = c("IL7R", "S100A4"), log = TRUE) #Memory CD4+
VlnPlot(seurat_obj, features = c("CD14", "LYZ"), log = TRUE) #CD14+ Mono
VlnPlot(seurat_obj, features = c("MS4A1"), log = TRUE) #B
VlnPlot(seurat_obj, features = c("CD8A"), log = TRUE) #CD8+ T
VlnPlot(seurat_obj, features = c("FCGR3A", "MS4A7"), log = TRUE) #FCGR3A+ Mono
VlnPlot(seurat_obj, features = c("GNLY","NKG7"), log = TRUE) #NK
VlnPlot(seurat_obj, features = c("FCER1A", "CST3"), log = TRUE) #Dendritic Cells
VlnPlot(seurat_obj, features = c("PPBP"), log = TRUE) #Platelet

top10[top10$gene %in% c("IL7R","CCR7", "CD14", "LYZ", "MS4A1", "CD8A","FCGR3A", "MS4A7","GNLY","NKG7", "FCER1A", "CST3", "PPBP"),]


norm_data <-as.matrix(seurat_obj[["RNA"]]@data)
dim(norm_data)
meta.data.cluster <- seurat_obj@meta.data
dim(meta.data.cluster)
table(meta.data.cluster$seurat_clusters == 0)
rownames(subset(meta.data.cluster, seurat_clusters = 0))
head(norm_data)[1:3, 1:3]
cl_0 <-  norm_data[,colnames(norm_data) %in% rownames(subset(meta.data.cluster, seurat_clusters == 0))]
dim(cl_0)
avr_cl_0 <- as.data.frame(rowMeans(cl_0)) 


?WhichCells

for(group in meta.data.cluster) {
  group.cells <- WhichCells(object = seurat_obj, subset.name = "seurat_clusters" , accept.value = group)
  data_to_write_out <- as.data.frame(x = as.matrix(x = seurat_obj[["RNA"]]@data)[, group.cells])
  write.csv(x = data_to_write_out, row.names = TRUE, file = paste0(group, "_Cluster_outfile.csv"))
}


saveRDS(seurat_obj, file = "pbmc3k_final.rds")
