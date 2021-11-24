#############################################################
### sub-type marker identification
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
#library(clusterProfiler)
##library(SingleR)


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################

outFolder="./5_harmonySubTypesDGE/"
system(paste0("mkdir -p ", outFolder))



anno<-read_rds("3_MergeDemux_Output/anno.rds")
sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")


dim(sc)

table(sc$Library)

table(sc$Location) 

system(paste0("mkdir -p ", outFolder))


clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

names(clust2Names)<-as.character(c(0:23))
cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte", "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



SMC<-names(clust2Names)[which(clust2Names %in% c("Smooth muscle cells-1","Smooth muscle cells-2","Smooth muscle cells-3"))]
Macrophages<-names(clust2Names)[which(clust2Names %in% c("Macrophage-1","Macrophage-2","Macrophage-3","Macrophage-4"))]
Macrophages<-names(clust2Names)[which(clust2Names %in% c("Macrophage-1","Macrophage-2","Macrophage-3","Macrophage-4"))]
Stromal<-names(clust2Names)[which(clust2Names %in% c("Stromal-1","Stromal-2","Myofibroblast"))]

#sc1 <- subset(sc, subset = seurat_clusters %in% SMC)
#sc1 <- subset(sc, subset = seurat_clusters %in% Macrophages)
sc1 <- subset(sc, subset = seurat_clusters %in% Stromal)

# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers$Celltype<-as.character(clust2Names[as.character(markers$cluster)])


m2 <- markers %>% left_join( dplyr::select(anno,gene=kbid,symbol=gene_name,rs)) 

m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)


top100 <- m2 %>% top_n(n = 100, wt = avg_log2FC)

fname=paste0(outFolder,"ClusterDEG.tsv");
write_tsv(m2,fname)
write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))


 fname=paste0(outFolder,"ClusterDEGtop100.tsv");
 write_tsv(top100,fname)
 write.csv(top100,file=paste0(outFolder,"ClusterDEGtop100.csv"))
 

#################
 
 outFolder="./5_harmonySubTypes_v2_DGE/"
 system(paste0("mkdir -p ", outFolder))
 
 
subtypes<-c("SMC","Macrophage","Stromal")
allmarkers<-read.csv(file=paste0("5_harmonyClustersDGE/ClusterDEG.csv"),stringsAsFactors = FALSE)
allmarkers<-allmarkers %>%filter(p_val_adj<0.05)
allmarkers$Cell_type<-clust2Names[as.character(allmarkers$cluster)]

SMC<-names(clust2Names)[which(clust2Names %in% c("Smooth muscle cells-1","Smooth muscle cells-2","Smooth muscle cells-3"))]
Macrophage<-names(clust2Names)[which(clust2Names %in% c("Macrophage-1","Macrophage-2","Macrophage-3","Macrophage-4"))]
Stromal<-names(clust2Names)[which(clust2Names %in% c("Stromal-1","Stromal-2","Myofibroblast"))]

subtypes_list<-list(SMC,Macrophage,Stromal)

 sapply(1:length(subtypes),function(x){
   
   print(subtypes[x])
   sc1 <- subset(sc, subset = seurat_clusters %in% subtypes_list[[x]])
   markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   markers$Celltype<-as.character(clust2Names[as.character(markers$cluster)])
   m2 <- markers %>% left_join( dplyr::select(anno,gene=kbid,symbol=gene_name,rs)) 
   m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)
   
   outFolder<-paste0(outFolder,subtypes[x],"/")
   system(paste0("mkdir -p ", outFolder,subtypes[x],"/"))
   allmarkers_subtypes<-allmarkers %>% filter(!cluster %in% unique(m2$cluster)) 
   
   
   m2<-m2 %>% filter(!gene %in% allmarkers_subtypes$gene)
   
   fname=paste0(outFolder,"ClusterDEG.tsv");
   write_tsv(m2,fname)
   write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))
   
   print(dim(m2))
   top100 <- m2 %>% top_n(n = 100, wt = avg_log2FC)
   
   fname=paste0(outFolder,"ClusterDEGtop100.tsv");
   write_tsv(top100,fname)
   write.csv(top100,file=paste0(outFolder,"ClusterDEGtop100.csv"))
 })
 

 