######################################################
# Single cell RNA-seq analysis with Seurat R package. 
# Analysis that was used in the publication: 
# 
# Written by: Dr. Noa Novershtern
# Runs on R.4.2.2
######################################################
library(enrichR)
library(cowplot)
library(dplyr)
library(patchwork)
library(erer)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(pheatmap)
library(data.table)
library(limma)

set.seed(1234)

##### Read seurat object #####
sem<-readRDS("sem.rds")

##### Annotate seurat clusters #####
sem$Epiblast <- ifelse((sem$seurat_clusters=="7" |sem$seurat_clusters=="1" | sem$seurat_clusters=="4"| sem$seurat_clusters=="11"),"Epi","rest")
sem$primed <- ifelse((sem$seurat_clusters=="7" ),"Primed","rest")
sem$PostEpi <- ifelse((sem$seurat_clusters=="4" ),"PostEpi","rest")
sem$Amnion <- ifelse(sem$seurat_clusters=="10","Amnion","rest")
sem$STB<- ifelse(sem$seurat_clusters=="12","STB","rest")
sem$ExEM <- ifelse((sem$seurat_clusters=="0" |sem$seurat_clusters=="6" | sem$seurat_clusters=="2"|sem$seurat_clusters=="8"),"ExEM","rest")
sem$YS <- ifelse((sem$seurat_clusters=="9" | sem$seurat_clusters=="3"| sem$seurat_clusters=="5"),"PrE","rest")
sem$SYS <- ifelse((sem$seurat_clusters=="9"),"SYS","rest")

DimPlot(sem, reduction = "umap", group.by="Epiblast",label=F,pt.size = 0.1,cols=c("brown1","grey")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="primed",label=F,pt.size = 0.1,cols=c("brown1","grey")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="PostEpi",label=F,pt.size = 0.1,cols=c("brown1","grey")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="Amnion",label=F,pt.size = 0.1,cols=c("brown1","grey")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="STB",label=F,pt.size = 0.1,cols=c("grey","brown1")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="YS",label=F,pt.size = 0.1,cols=c("brown1","grey")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="SYS",label=F,pt.size = 0.1,cols=c("grey","brown1")) & NoLegend() & NoAxes() 
DimPlot(sem, reduction = "umap", group.by="ExEM",label=F,pt.size = 0.1,cols=c("brown1","grey")) & NoLegend() & NoAxes() 

#saveRDS(sem, file = "sem.rds")


######### Define cell type features ############
epiblast_features=c("POU5F1","NANOG","SOX2","GDF3","SFRP2") 
com_epi_features = c("OTX2","ZIC2","ZEB2","VIM") 
ps_epi_features = c("TBXT","MIXL1", "MESP1", "EOMES", "WNT8A","SNAI2")

ys_ExEM_features = c("GATA6", "GATA4", "PDGFRA") 
ys_features = c("LINC00261", "SOX17","APOA1","OTX2") 
sys_features = c("TTR", "APOB","GSTA1")
ave_features = c("CER1","DKK1","LHX1") 

amnion_features = c("ISL1","GABRP","VTCN1","TFAP2A","GATA3","PRKD1","IGFBP3","KCNMA1","TFAP2C") 
exmc_features = c("FOXF1","VIM","POSTN","BST2","GATA6", "HAND2", "TBX4", "HGF" ) 

stb_features = c("GATA3","SDC1","CPM","RAI14","GADD45G","DEPP1","TREM1","HOPX","HEXIM1","HSPB8","ERVW-1","CGB5","CGB8","ADAM12","SORBS1","CSGALNACT1")
ctb_features = c("GABRP","PAGE4","LGALS3","OVOL1","SIGLEC6","NR2F2","KRT7","CCKBR")

####### Features plot ########

FeaturePlot(sem,features=epiblast_features,order = T,min.cutoff = 0,max.cutoff = 100)  & NoAxes()
FeaturePlot(sem,features=com_epi_features,order = T,min.cutoff = 0,max.cutoff = 10) & NoAxes()
FeaturePlot(sem,features=ps_epi_features,order = T,min.cutoff = 0,max.cutoff = 100) & NoAxes()

FeaturePlot(sem,features=ys_ExEM_features,order = T,min.cutoff = 0,max.cutoff = 100) & NoAxes()
FeaturePlot(sem,features=ys_features,order = T,min.cutoff = 0,max.cutoff = 100) & NoAxes()
FeaturePlot(sem,features=sys_features,order = T,min.cutoff = 0,max.cutoff = 100) & NoAxes()
FeaturePlot(sem,features=ave_features,order = T,min.cutoff = 0,max.cutoff = 100)  & NoAxes()

FeaturePlot(sem,features=amnion_features,order = T,min.cutoff = 0,max.cutoff = 100) & NoAxes()
FeaturePlot(sem,features=exmc_features,order = T,min.cutoff = 0,max.cutoff = 100)  & NoAxes()

FeaturePlot(sem,features=stb_features,order = T,min.cutoff = 0,max.cutoff = 100)  & NoAxes() 
FeaturePlot(sem,features=ctb_features,order = T,min.cutoff = 0,max.cutoff = 100) & NoLegend() & NoAxes()

####### Dot plot ########

###### Read cluster order file ########  
cluster_order<- read.table("SEM_order_new.csv",header=T)
sem$cluster_order <- cluster_order$order

## prepare dot plot (step 1)
dot_epi = DotPlot(sem,features = epiblast_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #6
dot_stb = DotPlot(sem,features = stb_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #7
dot_com_epi = DotPlot(sem,features = com_epi_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #3
dot_ys = DotPlot(sem,features = ys_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #7
dot_sys = DotPlot(sem,features = sys_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #4
dot_ave = DotPlot(sem,features = ave_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #3
dot_amnion = DotPlot(sem,features = amnion_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #3
dot_exmc = DotPlot(sem,features = exmc_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #8
dot_ps_epi = DotPlot(sem,features = ps_epi_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #4
dot_ctb = DotPlot(sem,features = ctb_features,group.by = "cluster_order",cluster.idents = F,scale =T,scale.min =0,assay="RNA") #6

## prepare dot plot (step 2)
dot_epi = dot_epi + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5)) + theme(aspect.ratio = 2.4) #dot epiblast
dot_ps_epi = dot_ps_epi + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 3.6) #dot sys
dot_com_epi = dot_com_epi + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 4.8) #dot primed
dot_ys = dot_ys + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 2.4) #dot ys
dot_sys = dot_sys + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 3.6) #dot sys
dot_ave = dot_ave + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+theme(aspect.ratio = 4.8) #dot ave
dot_exmc = dot_exmc + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 2) #dot exmc
dot_amnion = dot_amnion + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 2.4) #dot amnion
dot_stb = dot_stb + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5)) + theme(aspect.ratio = 1.2) #dot epiblast
dot_ctb = dot_ctb + theme(legend.text = element_text(family = 'Helvetica', size = 8)) + NoLegend() +   theme(axis.text.x = element_text(size=14, angle=90,hjust = 1,vjust = 0.5))+ theme(aspect.ratio = 2.4) #dot epiblast

## print dot plot
dot_epi
dot_ps_epi
dot_com_epi
dot_ys
dot_sys
dot_ave
dot_exmc
dot_amnion
dot_stb
dot_ctb

####### Find and show DEGs between ExEM subtypes ########
library(ggrepel)
DEGs = FindMarkers(sem,ident.1 = c("2","8") ,ident.2 = c("6","0"))

write.csv((DEGs),"Early_vs_Late_ExEM_markers.csv")
DEGs[order(DEGs$avg_log2FC, decreasing = T),] 

DEGs$sig = abs(DEGs$avg_log2FC)>1&DEGs$p_val_adj<0.05

DEGs$pval_c <- as.numeric(ifelse(DEGs$p_val_adj ==0,"1e-300",DEGs$p_val_adj))

ggplot(DEGs,aes(x=avg_log2FC,y=-log10(pval_c))) + 
  geom_point(size=1)+
  geom_text_repel(aes(label=ifelse(sig,rownames(DEGs),'')),segment.size = 0.1,box.padding = 0.2,force = 2,size=2)+
  geom_point(aes(color=sig))+
  scale_color_manual(values = c("black","red"))+
  theme_classic()




