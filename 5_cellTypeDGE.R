#############################################################
### cluster marker identification
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
#library(clusterProfiler)
##library(SingleR)

clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################

outFolder="./5_harmonyClustersDGE/"
system(paste0("mkdir -p ", outFolder))



anno<-read_rds("3_MergeDemux_Output/anno.rds")
sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")


dim(sc)

table(sc$Library)

table(sc$Location) 


# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


m2 <- markers %>% left_join(select(anno,gene=kbid,symbol=gene_name,rs)) 

m2 <- m2 %>% arrange(cluster,-avg_logFC) %>% group_by(cluster)

fname=paste0(outFolder,"ClusterDEG.tsv");
write_tsv(m2,fname)
m2<-read_tsv(fname)
write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))


m2$celltype<-clust2Names[as.character(m2$cluster)]
write.csv(m2,file=paste0(outFolder,"ClusterDEG_withcelltypes.csv"))




top20 <- m2 %>% top_n(n = 20, wt = avg_logFC)


# fname=paste0(outFolder,"ClusterDEGtop20.tsv");
# write_tsv(top20,fname)
top20<-read_tsv(fname)
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop20.csv"))

top20$celltype<-clust2Names[as.character(top20$cluster)]
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop20_withcelltypes.csv"))


