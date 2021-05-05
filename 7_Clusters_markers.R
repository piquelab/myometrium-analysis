
library(tidyverse)
##library(knitr)
library(DESeq2)
##library(annotables)
library(qqman)

library(Seurat)


library("BiocParallel")
register(MulticoreParam(12))


# outFolder <- paste0("7_outputs_DESeq_ConditionsByCluster/")
# system(paste0("mkdir -p ",outFolder))

## Load cells
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")

md <- sc@meta.data

# filter sample "HPL20874"


# md <- filter(md,!Pregnancy_ID %in% filter_sample)
# #md <- filter(md,!Pregnancy_ID %in% c("HPL20874","HPL20875"))
# 
# # 
# length(md$Pregnancy_ID=="HPL20874")
# 
# length(md$Pregnancy_ID %in% c("HPL20874","HPL20875"))

## Load gene annotations. 
anno <- read_rds("3_MergeDemux_Output/anno.rds")

stopifnot(identical(rownames(sc),anno$kbid))

gene_symbol <- anno$gene_name
names(gene_symbol) <- anno$kbid

cluster_markers<-lapply(c(0:23),function(ii)
    {
    cluster.markers <- FindMarkers(sc, ident.1 = ii, min.pct = 0.25)
    myres <- as.data.frame(cluster.markers) %>%
        rownames_to_column("kbid") %>%
        left_join(anno)
    
    cat("# cluster" , ii, "\n")
    return(myres)
    }
       )
names(cluster_markers)<-c(0:23)


outFolder <- paste0("5_harmonyClustersDGE/")
system(paste0("mkdir -p ",outFolder))
#summary_cluster_markers<-()
write_rds(cluster_markers,file=paste0(outFolder,"cluster_markers.rds"))

cluster_markers<-lapply(c(1:24),function(ii)
    {
    res<-cluster_markers[[ii]]  
    res<-res  %>% filter(p_val_adj<0.1)
    res<-res  %>% arrange(desc(-log10(p_val_adj)),desc(avg_log2FC) )
    res$cluster<-ii-1
    return(res)
    })

summary_cluster_markers_all <- do.call(rbind,cluster_markers)
summary_cluster_markers_all<-summary_cluster_markers_all[,rev(colnames(summary_cluster_markers_all))]

write.csv(summary_cluster_markers_all,file=paste0(outFolder,"summary_cluster_markers_all.csv"))


cluster_markers_top20<-lapply(c(1:24),function(ii)
{
    res<-cluster_markers[[ii]]  
    res<-res  %>% filter(p_val_adj<0.1)
    res<-res  %>% arrange(desc(-log10(p_val_adj)),desc(avg_log2FC) )
    res$cluster<-ii-1
    res<-res[1:20, c("cluster","gene_name","avg_log2FC","p_val_adj","pct.2","pct.1","rs")]
    return(res)
})

summary_cluster_markers_top20 <- do.call(rbind,cluster_markers_top20)
write.csv(summary_cluster_markers_top20,file=paste0(outFolder,"summary_cluster_markers_top20.csv"))
#write_tsv(summary_cluster_markers_top10,file=paste0(outFolder,"summary_cluster_markers_top10.tsv"))

################################################################################ 
#              visualization
################################################################################

sc.markers.all <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

sc.markers.all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- sc.markers.all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#heatmap plot
DoHeatmap(sc, features = top10$gene) + NoLegend()


#known markers across classters 

anno$kbid[which(anno$gene_name =="EMILIN2")]
fname=paste0(outFolder,"EMILIN2.pdf");
pdf(fname,width=10,height=6)
#volcano plot
VlnPlot(sc, features = c("ENSG00000132205.11"))


#feature plot
FeaturePlot(sc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))


#############################
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")
markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

outFolder <- paste0("5_harmonyClustersDGE/")
anno <- read_rds("3_MergeDemux_Output/anno.rds")
#markers<-read_rds(paste0(outFolder,"cluster_markers.rds"))


m2 <- markers %>% left_join(select(anno,gene=kbid,symbol=gene_name,rs)) 

m2 <- m2 %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)

#m2<-read_tsv(paste0(outFolder,"ClusterDEG.tsv"))

top20 <- m2 %>% top_n(n = 20, wt = avg_log2FC)

fname=paste0(outFolder,"ClusterDEG.tsv");
write_tsv(m2,fname)


fname=paste0(outFolder,"ClusterDEGtop20.tsv");
write_tsv(top20,fname)


## Try to do the cluster profiler plot. It would be good to do it for cell-type thing to. 
## http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Single_cell_markers.txt

cell_markers <- read_tsv('Human_cell_markers.txt')

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


####################################################################


library(Seurat)

clusters <- unique(as.character(Idents(sc)))
pairwise <- combn(clusters, 2)

results <- list()
top20<-list()
for(i in 1:ncol(pairwise)) {
    cat(i)
    markers <- FindMarkers(sc, ident.1 = pairwise[1, i], ident.2 = pairwise[2, i],only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    #comparisons <- pairwise[, i]
    #markers$comparison <- paste(comparisons[1], comparisons[2], sep = '_')
    markers$clustr1<-pairwise[, i][1]
    markers$clustr2<-pairwise[, i][2]
    markers$gene<-rownames(markers)
    
    m2 <- markers %>% left_join(select(anno,gene=kbid,symbol=gene_name,rs)) 
    
    m2 <- m2 %>% arrange(-avg_log2FC) 
    
    top20[[i]]<- m2 %>% top_n(n = 20, wt = avg_log2FC)
    results[[i]] <- m2
}

# write_rds(results,file=paste0(outFolder,"pairwise_markers_results.rds"))
# write_rds(top20,file=paste0(outFolder,"pairwise_markers_top20.rds"))

results.df <- do.call(rbind, results)
top20.df <- do.call(rbind, top20)

write_rds(results.df,file=paste0(outFolder,"pairwise_markers_results_df.rds"))
write_rds(top20.df,file=paste0(outFolder,"pairwise_markers_top20_df.rds"))

