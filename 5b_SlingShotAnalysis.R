
library(Seurat)
library(Matrix)
library(tidyverse)

library(future)


#################
##library(SingleR)

clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

#clust2Name<-paste0(c(0:23),"_",clust2Name)
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


outFolder="5b_SlingShotAnalysis/"
system(paste0("mkdir -p ", outFolder))


sc@meta.data$cluster_name <- clust2Name[sc@meta.data$seurat_clusters]


dim(sc)

## table(sc$Library)

## table(sc$Location) 

## table(md$scLabor_ID, md$Location)


## table(md$scLabor_ID, md$Origin)

## table(md$scLabor_ID, md$Location)

## table(md$SNG.BEST.GUESS, md$scLabor_ID)


## cc <- md %>% select(seurat_clusters,scLabor_ID) %>%
##     group_by(seurat_clusters,scLabor_ID) %>%
##     summarize(n=n()) %>%
##     group_by(seurat_clusters) %>%
##     mutate(prop=n/sum(n)) %>%
##     arrange(-n)



aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status","SNG.BEST.GUESS")) 
aa$cluster_name <- clust2Name[as.character(aa$seurat_clusters)]
# aa<-aa[order(as.numeric(aa$seurat_clusters),decreasing =TRUE),]
# aa$seurat_clusters <- factor(aa$seurat_clusters,levels=unique(aa$seurat_clusters))
# aa$cluster_name <- factor(aa$cluster_name,levels=unique(aa$cluster_name))



#UMAP scatter plot 

fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
dev.off()



cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")


# UMAP with colors
fname=paste0(outFolder,"UMAP_LocationHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.1) +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    ##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()

##################
##################
##library(scater)

library(slingshot)


selCells <- c("Stromal-1","Endothelial-1","Decidual","LED","Stromal-2","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","Smooth muscle cells-2","Smooth muscle cells-3","Macrophage-4","B-cell")

ssc = subset(scilab, idents = selCells)


## ## Maybe not needed or embed from PC to UMAP.
## sce <- Seurat::as.SingleCellExperiment(sc)

## sce <- slingshot(sce, clusterLabels = 'cluster_name', reducedDim = 'PCA')

## ssds <- SlingshotDataSet(sce)

## ##c1 <- getCurves(ssds)


## ### Maybe only need this part. ad not the top. 


sds <- slingshot(Embeddings(sc,"umap") , clusterLabels = sc$cluster_name,
                                  start.clus = "Stromal-1", stretch = 0)



## Make a simple plot here:
##fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
##png(fname,width=1600,height=1200)


####pdf("./SlingshotTrajectory.pdf")
plot(reducedDim(sds), col = cluster.Colors[sc$cluster_name], pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
##lines(sds, lwd = 2, col = 'black',)
dev.off()


##res.0.7 , start 3 blue, or 15 yellow


sds2 <- slingshot(Embeddings(troph,"umap") , clusterLabels = troph$res.0.7,
                                  start.clus = 3, stretch = 0)





mycol = RColorBrewer::brewer.pal(8, "Set1")
names(mycol)= unique(troph$res.0.7)
mycol

table(troph$res.0.7,troph$FinalNames)



pdf("./SlingshotTrajectory.v2.pdf")
plot(reducedDim(sds2), col = mycol[troph$res.0.7], pch = 16, cex = 0.5)
lines(sds2, lwd = 2, type = 'lineages', col = 'black')
##lines(sds, lwd = 2, col = 'black',)
dev.off()

pdf("./SlingshotTrajectory.v3.pdf")
plot(reducedDim(sds2), col = group.colors[troph$FinalNames], pch = 16, cex = 0.5, xlim=c(-10,6),ylim=c(3,15))
lines(sds2, lwd = 2, type = 'lineages', col = '#000000AA')
##lines(sds2, lwd = 2, col = 'black',)
dev.off()


####




## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=10),title="Condition")) +
     theme(legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))+
    scale_color_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()




## Make a simple plot here:
fname=paste0(outFolder,"UMAP_OriginHarmony_v2.png");
#fname=paste0(outFolder,"UMAP_OriginHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    theme_bw()+
    guides(colour = guide_legend(override.aes = list(size=5),title="Origin")) +
    #theme(legend.text=element_text(size=50,face="bold"), legend.title =element_text(size=50,face="bold") , axis.text=element_text(size=50,face="bold"))+
    theme(text = element_text(size=30,face = "bold"),
          plot.title = element_text(size = 25, face = "bold"),
          legend.title=element_text(size=25,face="bold"), 
          legend.text=element_text(size=25,face="bold"))+
    
   # scale_color_manual(values=c("M"="#DE70EA","F"="#A61BB5"))
scale_color_manual(values=c("M"="#D1D1D1","F"="#A61BB5"))

##    facet_wrap(~LocTime) +
    
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





#changed 

clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

names(clust2Name)<-c(0:23)


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status","SNG.BEST.GUESS")) 
aa$cluster_name <- clust2Name[as.character(aa$seurat_clusters)]
aa<-aa[order(as.numeric(aa$seurat_clusters),decreasing =TRUE),]
aa$seurat_clusters <- factor(aa$seurat_clusters,levels=unique(aa$seurat_clusters))
aa$cluster_name <- factor(aa$cluster_name,levels=unique(aa$cluster_name))


fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=reorder(cluster_name,-seurat_clusters),fill=Condition)) +
    geom_bar(position = position_stack(reverse = TRUE)) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Location) + coord_flip() +
    scale_fill_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
    xlab("")+
    theme_bw()
p2
##    theme_black()
dev.off()

# fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
# pdf(fname,width=10,height=6)
# #p2 <- ggplot(aa,aes(x=reorder(cluster_name,-seurat_clusters),fill=Condition)) +
# p2 <- ggplot(aa,aes(x=cluster_name ,y=Condition,fill=Condition)) +  
#     #geom_bar(position = position_stack(reverse = TRUE)) +
#     geom_bar (stat="identity", width = 0.9,position =position_dodge(width = 0.8))+
#     ##    scale_color_manual(values=group.colors) +
#     #guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
#     facet_grid(.~Location) + coord_flip() +
#     scale_fill_manual(values=c("TNL"="#333399","TIL"="#A50021"))+
#     xlab("")+
#     theme_bw()
# p2
# ##    theme_black()
# dev.off()


### END- HERE ###
########################################################


