## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
##library(SingleR)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

load("3_MergeDemux_Output/annoOnly.Rdata")

sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")



outFolder="./5_harmonyClustersDGE/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 


## DGE.
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


## Try to do the cluster profiler plot. It would be good to do it for cell-type thing to. 
## http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Single_cell_markers.txt

##cell_markers <- read_tsv('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') 
cell_markers <- read_tsv('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Single_cell_markers.txt')

cm <- cell_markers %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
cm

library(clusterProfiler)

##genelist = m2 %>% filter(cluster==4) %>% ungroup() %>% dplyr::select(symbol) %>% unlist()
##genelist2 <- genelist %>%  bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
##y <- enricher(genelist2$ENTREZID, TERM2GENE=cm, minGSSize=1)

m3 <- left_join(m2,bitr(m2$symbol,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% dplyr::select(symbol=SYMBOL,ENTREZID)) %>%
    filter(!is.na(ENTREZID))

## cluster profiler plot?
cres <- compareCluster(ENTREZID ~ cluster, data=m3, fun="enricher", TERM2GENE=cm, minGSSize=1)

fname=paste0(outFolder,"dotplot.pdf");
pdf(fname,width=15,height=15)
dotplot(cres)
dev.off()


cell_markers <- read_tsv('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') 
##cell_markers <- read_tsv('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Single_cell_markers.txt')

cm <- cell_markers %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
cm

m3 <- left_join(m2,bitr(m2$symbol,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% dplyr::select(symbol=SYMBOL,ENTREZID)) %>%
    filter(!is.na(ENTREZID))

## cluster profiler plot?
cres <- compareCluster(ENTREZID ~ cluster, data=m3, fun="enricher", TERM2GENE=cm, minGSSize=1)

fname=paste0(outFolder,"dotplot.celltype.pdf");
pdf(fname,width=15,height=15)
dotplot(cres)
dev.off()


### Top 20
m3 <- left_join(top20,bitr(m2$symbol,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% dplyr::select(symbol=SYMBOL,ENTREZID)) %>%
    filter(!is.na(ENTREZID))

## cluster profiler plot?
cres <- compareCluster(ENTREZID ~ cluster, data=m3, fun="enricher", TERM2GENE=cm, minGSSize=1)

fname=paste0(outFolder,"dotplot.top20.celltype.pdf");
pdf(fname,width=15,height=15)
dotplot(cres)
dev.off()

## should try other types of enrichment...
## or thin the bigger cell-type definitions as it is not helpful... 


##pbmc <- RenameIdents(pbmc, new.cluster.ids)
fname=paste0(outFolder,"UMAP_Harmony.pdf");
pdf(fname,width=5,height=5)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


##table(sc@meta.data$seurat_clusters,sc@meta.data$SNG.BEST.GUESS)



### END- HERE ###
########################################################


