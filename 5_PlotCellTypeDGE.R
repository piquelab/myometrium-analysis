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

anno <- read_rds("3_MergeDemux_Output/anno.rds")
#sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")


outFolder="./5_PlotClustersDGE/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 


m2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")

Hmax=log2(max(m2$cluster)+1)

m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
  group_by(gene) %>%
  mutate(H=log2(length(cluster))) %>%
  filter(H<=1) %>%
  ungroup()

table(m3$cluster)


top20 <- m3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% ungroup()
table(top20$cluster)


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Rep","SNG.BEST.GUESS"))

myscale = 1/colSums(sc@assays$RNA@counts)*1000000

for(i in 0:max(top20$cluster)){
  ##    
  genesel <- filter(top20,cluster==i)
  ##        FetchData(sc,genesel$gene)
  genexpr <- map_dfr(1:nrow(genesel),function(i){
    aux <- FetchData(sc,genesel$gene[i])
    colnames(aux) <- "expr"
    aux <- cbind(aa,aux*myscale)
    aux$gene=genesel$gene[i]
    aux$symbol=genesel$symbol[i]
    aux
  })
  ##    
  cat(i,dim(genexpr),"\n")    
  ##
  fname=paste0(outFolder,"cluster_",i,"_umap.png");
  ##png(fname,width=1000,height=800)
  p1 <- genexpr %>% arrange(symbol,expr) %>%
    ggplot(aes(UMAP_1,UMAP_2,color=log10(0.1+expr))) +
    geom_point(size=0.1) +
    ##        scale_fill_viridis_c(option = "plasma") +
    scale_color_gradient(low = "lightblue", high = "darkred") +
    facet_wrap(~symbol) +
    theme_bw()
  ggsave(fname,p1,width=10,height=8)
  ##dev.off()
  cat(fname,"\n")
}

top20 <- m3 %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC) %>% ungroup()
table(top20$cluster)


for(i in 0:max(top20$cluster)){
  ##    
  genesel <- filter(top20,cluster==i)
  ##        FetchData(sc,genesel$gene)
  genexpr <- map_dfr(1:nrow(genesel),function(i){
    aux <- FetchData(sc,genesel$gene[i])
    colnames(aux) <- "expr"
    aux <- cbind(aa,aux*myscale)
    aux$gene=genesel$gene[i]
    aux$symbol=genesel$symbol[i]
    aux
  })
  ##    
  cat(i,dim(genexpr),"\n")
  ##
  ##    sum_rec <- genexpr %>% dplyr::group_by(symbol,Location,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))
  sum_rec <- genexpr %>% dplyr::group_by(symbol,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))
  ##
  fname=paste0(outFolder,"cluster_",i,"_dotplot2.png");
  ##    png(fname,width=1600,height=1000,pointsize=48)
  ##    p0 <- ggplot(sum_rec,aes(x=Location,y=seurat_clusters,color=Expr,size=Prop)) +
  p0 <- ggplot(sum_rec,aes(x=seurat_clusters,y=symbol,color=Expr,size=Prop)) +
    geom_point() +
    scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
    scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
    ##facet_grid( .~ symbol) +
    theme_classic() +
    ##        labs(x="Tissue",y="Cell type",size="Prop.",color="Expr.(TPM)") +
    labs(y="Gene Symbol",x="Cell type",size="Prop.",color="Expr.(TPM)") + 
    ## theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45,hjust=1)) 
  ##dev.off()
  ##    ggsave(fname,p0,width=16,height=10)
  ggsave(fname,p0,width=8,height=10)
  cat(fname,"\n")
}




clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")

names(clust2Name)<-c(0:23)
clust2Name<-paste0(names(clust2Name),"_",clust2Name)

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.V2.pdf");
pdf(fname,width=10,height=6)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Rep","SNG.BEST.GUESS")) 
aa$seurat_clusters <- clust2Name[aa$seurat_clusters]
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
  geom_point(size=0.1) +
  ##    scale_color_manual(values=group.colors) +
  guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
  facet_wrap(~Location) +
  theme_bw()
p1
##    theme_black()
dev.off()


##pbmc <- RenameIdents(pbmc, new.cluster.ids)
fname=paste0(outFolder,"UMAP_Harmony.pdf");
pdf(fname,width=5,height=5)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


##table(sc@meta.data$seurat_clusters,sc@meta.data$SNG.BEST.GUESS)



################################
# heatmap


clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-paste0(c(0:23),"_",clust2Names)
#sc$seurat_clusters<-clust2Names[sc$seurat_clusters]


scale.data<-sc@assays$RNA@scale.data
rw<-anno$gene_name[which(anno$kbid %in% rownames(scale.data))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(scale.data))]
scale.data<-scale.data[names(rw),]
rownames(scale.data)<-as.character(rw)
sc@assays$RNA@scale.data<-scale.data
#colnames(scale.data)<-newsc$seurat_clusters[colnames(scale.data)]

m2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
  group_by(gene) %>%
  mutate(H=log2(length(cluster))) %>%
  filter(H<=1) %>%
  ungroup()

top10 <- m3 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

fname=paste0(outFolder,"Heatmap-top10.pdf");
pdf(fname,width=35,height=25)
#DoHeatmap(newsc, features = top10$symbol) + NoLegend()
DoHeatmap(sc, features = top10$symbol)+ #+ NoLegend()
theme(text = element_text(size = 10),axis.text.y = element_text(size = 10)) + NoLegend()#
dev.off()



outFolder="./feature_importance/"

scale.data<-sc@assays$RNA@scale.data
colnames(scale.data)<-sc$seurat_clusters[colnames(scale.data)]
rw<-anno$gene_name[which(anno$kbid %in% rownames(scale.data))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(scale.data))]
scale.data<-scale.data[names(rw),]
rownames(scale.data)<-as.character(rw)
scale.data2<-t(scale.data)
scale.data2<-as.data.frame(scale.data2)
scale.data2$target<-colnames(scale.data)
write.csv(scale.data2,file=paste0(outFolder,"scale.data2.csv"))

data<-sc@assays$RNA@data
colnames(data)<-sc$seurat_clusters[colnames(data)]
rw<-anno$gene_name[which(anno$kbid %in% rownames(data))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(data))]
data<-data[names(rw),]
rownames(data)<-as.character(rw)

data2<-t(as.matrix(data))
data2<-as.data.frame(data2)
data2$target<-colnames(data)
write.csv(data2,file=paste0(outFolder,"data2.csv"))
