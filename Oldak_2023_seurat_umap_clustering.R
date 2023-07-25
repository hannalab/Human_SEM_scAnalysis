######################################################
# Single cell RNA-seq analysis with Seurat R package. 
# Analysis that was used in the publication: 
# 
# Written by: Dr. Noa Novershtern
# Runs on R.4.2.2
######################################################
library(cowplot)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(pheatmap)
library(data.table)
set.seed(1234)

##### Load SEM experiments data and create Seurat object #####

SEM.data <- Read10X(data.dir = "<path_to_CellRanger_counts>/filtered_feature_bc_matrix/")
SEM.data <- CreateSeuratObject(counts = SEM.data, project = "SynHum", min.cells = 3, min.features = 200)

#####Add sample names ######
groups<-as.numeric(str_extract(colnames(SEM.data),"\\d+"))
unique(groups)
SEM.data$sample<-groups
type<- SEM.data$sample
type<-ifelse(groups==1,type[groups]<-"W3_SEM_Day4",type[groups]<-type)
type<-ifelse(groups==2,type[groups]<-"W3_SEM_Day6",type[groups]<-type)
type<-ifelse(groups==3,type[groups]<-"W3_SEM_Day8",type[groups]<-type)
unique(type)
SEM.data$sample<-type


###### QC samples ###### 
SEM.data[["percent.mt"]] <- PercentageFeatureSet(SEM.data, pattern = "^MT-")
VlnPlot(SEM.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,same.y.lims=F,log = T,split.by = "sample",group.by="sample",cols =c("coral1","limegreen","steelblue2"))+theme(legend.position = "none",plot.title = element_blank())+xlab("")

#### QC scatter plots ####
plot1 <- FeatureScatter(SEM.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SEM.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

##### Filter cells ######
SEM.data$percent.mt
SEM.data_sub <- subset(SEM.data, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & percent.mt < 15)

###### Normalize each sample #######
new_sub.list <- SplitObject(SEM.data_sub, split.by = "sample")
new_sub.list <- lapply(X = new_sub.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

##### Integrate the samples ######
dims=1:10
# The following step may take a few minutes
anchors <- FindIntegrationAnchors(object.list = new_sub.list, dims = dims,anchor.features = 2000)
sem <- IntegrateData(anchorset = anchors, dims = dims)
DefaultAssay(sem) <- "integrated"
# Run the standard workflow for visualization and clustering
sem <- ScaleData(sem, verbose = FALSE)
sem <- RunPCA(sem, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
sem <- RunUMAP(sem, reduction = "pca", dims = dims)
sem <- FindNeighbors(sem, reduction = "pca", dims = dims)
sem <- FindClusters(sem, resolution = 0.5)# seq(0,1.2,0.1))

ElbowPlot(sem)

##### A recommended step: save sem.rds, so the processed object can be read for future analysis #####
#saveRDS(sem, file = "sem.rds")

###### Visualization #######
DimPlot(sem, reduction = "umap", group.by = "sample",label=F,pt.size = 0.1)
DimPlot(sem, reduction = "umap", split.by="sample",label=F,pt.size = 0.1)
DimPlot(sem, reduction = "umap", group.by="seurat_clusters",label=F,pt.size = 0.1)

# Number of counts per cluster
VlnPlot(sem, features = c("nCount_RNA"), ncol = 3,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters") +theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(sem,features="nFeature_RNA",order = T) & NoAxes()

########## Save cluster's markers #############
SEM_markers = FindAllMarkers(sem, group.by="seurat_clusters",only.pos=T, logfc.threshold = 0.25)
write.csv(SEM_markers,"SEM_markers.csv")

