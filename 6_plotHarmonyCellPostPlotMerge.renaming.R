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

##sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")

##sc <- read_rds("./4_harmony/sc.NormByLibrary.Harmony,StringentFiltering.rds")
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")

## ## To Fix in previous script. 
## md <- sc@meta.data
## md <- md %>% column_to_rownames("BARCODES")
## identical(colnames(sc),rownames(md))
## ## md <- md[colnames(sc),]
## sc@meta.data <- md
## identical(colnames(sc),rownames(sc@meta.data))
## write_rds(sc,"./5_harmony_cellClass_plots_res0.8/SeuratObject.rds")

identical(colnames(sc),rownames(sc@meta.data))



##################################
# rename the clusters 
##################################

# new_names <- read_tsv("./6_harmony_cellClass_plots_res0.8_final/ClusterAssignment.tsv")
# clust2Names <- new_names$scLabor_ID
# names(clust2Names) <- new_names$seurat_clusters
# cc <- new_names %>% select(scLabor_ID,color) %>% unique 
# cluster.Colors <- cc$color
# names(cluster.Colors) <- cc$scLabor_ID
# 
# tempcol<-cluster.Colors["B-cell"]
# cluster.Colors["B-cell"]<-cluster.Colors["T-cell"]
# cluster.Colors["T-cell"]<-tempcol

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





sc@meta.data$cluster_name <- clust2Names[sc@meta.data$seurat_clusters]

table(sc$Library)

table(sc$Location) 

table(sc$cluster_name, sc$Location)

table(sc$cluster_name, sc$Origin)

table(sc$SNG.BEST.GUESS, sc$cluster_name)


outFolder="./6_harmony_rename_res0.8_plots/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

write_rds(sc,paste0(outFolder,"SeuratObject.rds"))

Idents(sc) <- "cluster_name"

##pbmc <- RenameIdents(pbmc, new.cluster.ids)
fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

aa <- FetchData(sc,c("UMAP_1","UMAP_2","Location","Condition","Origin","status","FetalSex","seurat_cluster","cluster_name")) 
head(aa)


## Make a simple plot here:
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


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    ##    facet_wrap(~LocTime) +
    scale_color_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))

p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_OriginHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Origin")) +
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
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.1) +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
    facet_grid(Condition ~ Location) +
    theme_bw()
p1
##    theme_black()
dev.off()




fname=paste0(outFolder,"UMAP_Location.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=cluster_name,fill=Origin)) +
    geom_bar(position="stack") +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
    facet_grid(.~Location) + coord_flip() +
    theme_bw()
p2
##    theme_black()
dev.off()


# fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
# pdf(fname,width=10,height=6)
# p2 <- ggplot(aa,aes(x=cluster_name,fill=Condition)) +
#     geom_bar(position="stack") +
# ##    scale_color_manual(values=group.colors) +
#     guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
#     facet_grid(.~Location) + coord_flip() +
#     theme_bw()
# p2
# ##    theme_black()
# dev.off()

#changed 
fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=cluster_name,fill=Condition)) +
    geom_bar(position = position_stack(reverse = TRUE)) +
    ##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Location) + coord_flip() +
    scale_fill_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))
p2
##    theme_black()
dev.off()

### END- HERE ###
########################################################


