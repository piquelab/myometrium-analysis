
#############################################################
###  TIL vs TNL across cell types

#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)

#################
##library(SingleR)


clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast", "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Name)<-c(0:23)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")


md <- read_rds("./4_harmony_cellClass_PBMC/sc.NormByLocation.ref.Anchors.rds") %>%
    as.data.frame %>%
    rownames_to_column("BARCODES") %>%
    select(BARCODES,scLabor_ID=predicted.celltype.l2,scLabor_Score=predicted.celltype.l2.score)

md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
    left_join(md) 

identical(md$BARCODES,rownames(sc@meta.data))


outFolder="5b_cellTypeAnalysisByCondition"
system(paste0("mkdir -p ", outFolder))


dim(sc)

table(sc$Library)

table(sc$Location) 


aa <- FetchData(sc,c("seurat_clusters","Location","Condition","Origin","status","Pregnancy_ID")) 
aa$seurat_clusters <- clust2Name[aa$seurat_clusters]


cc <- aa %>% group_by(seurat_clusters,Condition,Pregnancy_ID) %>%
    summarize(n=n()) %>%
    group_by(seurat_clusters) %>% mutate(p0=sum(n)/nrow(aa)) %>%
    group_by(Pregnancy_ID) %>%
    mutate(nt=sum(n),p=n/nt,z=(p-p0)/sqrt(p0*(1-p0)*(1/nt+1/nrow(aa))))


ccdiff <-  cc %>% group_by(seurat_clusters) %>% summarize(wilcox.pval = wilcox.test(p ~ Condition)$p.value,
                                               test.t = t.test(z ~ Condition)$statistic) %>% ungroup() %>%
    mutate(wilcox.padj = p.adjust(wilcox.pval))  

ccdiff_filtered<-ccdiff %>% filter(wilcox.padj<0.1)
write.csv(ccdiff,file=paste0(outFolder,"/wilcox_result.csv"))




aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status","SNG.BEST.GUESS")) 
aa$cluster_name <- clust2Name[aa$seurat_clusters]
aa$seurat_clusters<-as.numeric(aa$seurat_clusters)
aa<-aa[order(aa$seurat_clusters,decreasing =FALSE),]


ccdiff_filtered$seurat_clusters
aa$cluster_name[which(aa$cluster_name %in% ccdiff_filtered$seurat_clusters )]<-paste(aa$cluster_name," *",sep="")[which(aa$cluster_name %in% ccdiff_filtered$seurat_clusters )]


fname=paste0(outFolder,"/UMAP_LocationCondition.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=reorder(cluster_name,-seurat_clusters),fill=Condition)) +
    geom_bar(position = position_dodge()) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Location) + coord_flip() +
    scale_fill_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
    xlab("")+
    theme_bw()
p2
##    theme_black()
dev.off()




