######################################################
# Single cell RNA-seq analysis with Seurat R package. 
# Analysis that was used in the publication: 
# 
# Written by: Dr. Noa Novershtern
# Runs on R.4.2.2
######################################################
library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(terra)
library(monocle3)
library(ggplot2)
library(pheatmap)

set.seed(1234)

### read sem seurat object
sem<-readRDS("sem.rds")

### subset Epiblast clusters
Idents(sem) <- "seurat_clusters"
epi <- subset(x = sem, idents = c("1","4","7","11"))
epi$orig_cluster <- Idents(epi)
my_colors <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#E6AB02")

DimPlot(epi, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 1)
dims=1:15
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 20)
# t-SNE and Clustering
epi <- RunUMAP(epi, reduction = "pca", dims = dims)
epi <- FindNeighbors(epi, reduction = "pca", dims = dims)
epi <- FindClusters(epi, resolution = 0.2)# seq(0,1.2,0.1))
DimPlot(epi, reduction = "umap", group.by="seurat_clusters",label=F,pt.size = 0.5,cols = my_colors)
epi$orig.ident <- epi$seurat_clusters

#####
cds <- as.cell_data_set(epi)
cds <- cluster_cells(cds,resolution=0.001) # "leiden" clustering was used in the paper
cds <- cluster_cells(cds,cluster_method = "louvain",k = 200) # alternative clustering with louvain

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE,close_loop = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, group_label_size = 7,
           reduction_method="UMAP")

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == c(6)]),reduction_method = "UMAP")

plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black",
           reduction_method="UMAP")


##### print heatmap 
epi$orig.ident <- epi$seurat_clusters

committed_markers <- FindMarkers(epi, ident.1 = c("1"), ident.2 = c("0","2","4"), only.pos = TRUE)
posterior_markers <- FindMarkers(epi, ident.1 = c("3"), ident.2 = c("0","2","4"), only.pos = TRUE)

commit_features = c(rownames(committed_markers[order(committed_markers$avg_log2FC,decreasing=TRUE),][1:31,]))
commit_features_clean <- commit_features[!commit_features == "MALAT1"]
post_features = c(rownames(posterior_markers[order(posterior_markers$avg_log2FC,decreasing=TRUE),][1:30,]))


cell_annot <- as.data.frame(epi$orig.ident)
clusters = c("0","2","4","1")
epi$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime

epi_commited <- subset(epi, idents=c("0","4","2","1"))
epi_post <- subset(epi, idents=c("0","4","2","3"))

commit_cell_annot <- as.data.frame(epi_commited$seurat_clusters)
post_cell_annot <- as.data.frame(epi_post$seurat_clusters)

my_commit_colors <- list(`epi_commited$seurat_clusters` = c("0" = "#1B9E77", "1" = "#D95F02", "2" = "#7570B3","3" = "#E7298A","4" = "#E6AB02"))
my_post_colors <- list(`epi_post$seurat_clusters` = c("0" = "#1B9E77", "1" = "#D95F02", "2" = "#7570B3","3" = "#E7298A","4" = "#E6AB02"))

pheatmap(epi_commited@assays$RNA[commit_features_clean,names(sort(epi_commited$pseudotime))],cluster_cols = FALSE,fontsize = 6,annotation_col = commit_cell_annot,annotation_names_col = F,annotation_colors = my_commit_colors)
pheatmap(epi_post@assays$RNA[post_features,names(sort(epi_post$pseudotime))],cluster_cols = FALSE,fontsize = 6,annotation_col = post_cell_annot,annotation_names_col = FALSE,clustering_distance_rows = "minkowski",annotation_colors = my_post_colors)
