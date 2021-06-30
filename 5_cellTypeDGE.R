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


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################

outFolder="./5_harmonyClustersDGE/"
system(paste0("mkdir -p ", outFolder))



load("3_MergeDemux_Output/annoOnly.Rdata")
sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")


dim(sc)

table(sc$Library)

table(sc$Location) 


# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


m2 <- markers %>% left_join(select(anno,gene=kbid,symbol=gene_name,rs)) 

m2 <- m2 %>% arrange(cluster,-avg_logFC) %>% group_by(cluster)


top20 <- m2 %>% top_n(n = 20, wt = avg_logFC)

fname=paste0(outFolder,"ClusterDEG.tsv");
write_tsv(m2,fname)
m2<-read_tsv(fname)
write.csv(m2,file=paste0(outFolder,"ClusterDEG.csv"))

fname=paste0(outFolder,"ClusterDEGtop20.tsv");
write_tsv(top20,fname)
top20<-read_tsv(fname)
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop20.csv"))



