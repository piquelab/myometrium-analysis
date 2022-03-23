## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
##library(SingleR)

source("theme_black.R")

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
               "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
               "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

names(clust2Names)<-as.character(c(0:23))
cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


clusters_annotation_color<-cbind(names(clust2Names),as.character(clust2Names),as.character(cluster.Colors))
colnames(clusters_annotation_color)<-c("Cluster-ID","Cluster-label","Color")
#write.csv(clusters_annotation_color,file="clusters_annotation_color.csv")
###########################################
sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")


md <- read_rds("./4_harmony_cellClass_parturition_elife/sc.NormByLocation.ref.Harmony.singler.res0.8.rds") %>%
    as.data.frame %>%
    rownames_to_column("BARCODES") %>%
    select(BARCODES,pruned.labels,scLabor_ID=pruned.labels,scLabor_Score=tuning.scores.second)


md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
    left_join(md) 

identical(md$BARCODES,rownames(sc@meta.data))



outFolder="5_harmony_cellClass_parturition_elife_res0.8_plots/"
#outFolder="5_harmony_cellClass_parturition_elife_res0.6_plots/"
system(paste0("mkdir -p ", outFolder))

dim(sc)

table(sc$Library)

table(sc$Location) 

table(md$scLabor_ID, md$Location)


table(md$scLabor_ID, md$Origin)

table(md$scLabor_ID, md$Location)

table(md$SNG.BEST.GUESS, md$scLabor_ID)

pregtable <- table(md$Pregnancy_ID,sc$Library)



pregtable %>% as.data.frame %>% filter(Freq>0) %>% arrange(Var2) %>% write_tsv("pregtable.tsv")


cc <- md %>% select(seurat_clusters,scLabor_ID) %>%
    group_by(seurat_clusters,scLabor_ID) %>%
    summarize(n=n()) %>%
    group_by(seurat_clusters) %>%
    mutate(prop=n/sum(n)) %>%
    arrange(-n)

cc2 <- cc %>% top_n(n=1) %>% mutate(cname=paste0(seurat_clusters,":",scLabor_ID)) 
tt <- table(md$seurat_clusters,md$scLabor_ID)
tt <- tt/rowSums(tt)

library(pheatmap)

fname=paste0(outFolder,"Singler.HeatMap.pdf");
pdf(fname,width=7,height=4)
pheatmap(t(tt),cluster_rows=TRUE,cluster_cols=FALSE,scale="none")
dev.off()

##pbmc <- RenameIdents(pbmc, new.cluster.ids)

#UMAP scatter plot 

mycol=cluster.Colors
names(mycol)<-as.character(c(0:23))

fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
DimPlot(sc, cols=mycol,reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.png");
png(fname,width=1600,height=1200)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status")) 
aa$seurat_clusters <- clust2Names[aa$seurat_clusters]
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Location)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=8),title="Location")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()



## Condition
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Condition")) +
     theme(legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))+
    scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
##    facet_wrap(~LocTime) +
 theme_black() #theme_bw()
p1
##    theme_black()
dev.off()


## Condition
fname=paste0(outFolder,"UMAP_ConditionHarmony2.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Condition")) +
    theme(legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))+
    scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
    ##    facet_wrap(~LocTime) +
    
    theme_black() + #theme_bw()++
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p1
##    theme_black()
dev.off()



fname=paste0(outFolder,"UMAP_OriginHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Origin")) +
    theme(legend.text=element_text(size=20,face="bold"), legend.title =element_text(size=20,face="bold") , axis.text=element_text(size=20,face="bold"))+
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_StatusSoC_Harmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=status)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="status")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()



## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
    facet_wrap(Condition ~ Location) +
    theme_bw()
p1
##    theme_black()
dev.off()



fname=paste0(outFolder,"UMAP_Location.Barplot.pdf");
pdf(fname,width=10,height=6)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status","SNG.BEST.GUESS")) 
aa$seurat_clusters <- clust2Names[aa$seurat_clusters]

p2 <- ggplot(aa,aes(x=seurat_clusters,fill=SNG.BEST.GUESS)) +
    geom_bar(position="stack") +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
    facet_grid(.~Location) + coord_flip() +
    theme_bw()
p2
##    theme_black()
dev.off()

fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=SNG.BEST.GUESS)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
    facet_wrap(~Location) +
    theme_bw()
p2
##    theme_black()
dev.off()





### END- HERE ###
########################################################


