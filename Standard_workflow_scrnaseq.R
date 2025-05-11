## Standard workflow steps to analyse single cell RNAseq data

# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0

# setwd("F:/scRNAseq")

# Install Libraries
#install.packages("Seurat")
#install.packages("tidyverse")
#install.packages("hdf5r")

# load libraries
#library(Seurat)
#library(tidyverse)
#library(hdf5r)

# Load the NSCLC dataset
#nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
#str(nsclc.sparse.m)
#cts <-  nsclc.sparse.m$`Gene Expression`

# Initialize the Seurat object with the raw (non-normalized data).
#nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
#str(nsclc.seurat.obj)
# The Seurat object has 29552 features across 42081 samples

# 1. QC ############################
#View(nsclc.seurat.obj@meta.data)

# % MT reads
#nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
#View(nsclc.seurat.obj@meta.data)

# Create a violin plot
#pdf("Violin_plot_1.pdf", width = 7, height = 5) 
#VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

#Create scatter plot
#pdf("scatter_plot_2.pdf", width = 7, height = 5)
#FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
#geom_smooth(method = 'lm')
#dev.off()

# 2. Filtering ##########################
#nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
#                             percent.mt < 5)

# 3. Normalize data ######################
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
#str(nsclc.seurat.obj)

# 4. Identify highly variable features ##################
#nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)


# plot variable features with and without labels
#pdf("variable_features_3.pdf", width = 7, height = 5) 
#plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
#LabelPoints(plot = plot1, points = top10, repel = TRUE)
#dev.off()

# 5. Scaling #######################################
#all.genes <- rownames(nsclc.seurat.obj)
#nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
#str(nsclc.seurat.obj)


# 6. Perform Linear dimensionality reduction ##############################
#nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))


# visualize PCA results
#print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
#pdf("DimHeatmap_4.pdf", width = 7, height = 7)
#DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
#dev.off()


# determine dimensionality of the data
 
#pdf("Elbow_plot_5.pdf", width = 7, height = 5) 
#ElbowPlot(nsclc.seurat.obj) 
#dev.off()

# 7. Clustering ------------
#nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
#nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
#View(nsclc.seurat.obj@meta.data)

#pdf("Dimplot_6.pdf",width = 7, height = 5)
#DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
#dev.off()

# setting identity of clusters
#Idents(nsclc.seurat.obj)
#Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
#Idents(nsclc.seurat.obj)


#8. non-linear dimensionality reduction ################################

# install umap
#install.packages("reticulate")
#reticulate::install_miniconda()
#reticulate::conda_create("r-reticulate", packages = "python=3.8")
#reticulate::use_condaenv("r-reticulate", required = TRUE)
#reticulate::py_install("umap-learn", envname = "r-reticulate")
#reticulate::py_install(packages = 'umap-learn')
#reticulate::py_install("umap-learn")
# Reload Seurat and reticulate
#library(Seurat)
#library(reticulate)

#nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

#pdf("Dimplot_7.pdf",width = 7, height = 5)
#DimPlot(nsclc.seurat.obj, reduction = "umap")
#dev.off()





