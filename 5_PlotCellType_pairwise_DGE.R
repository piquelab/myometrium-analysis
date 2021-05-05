########################################################
# pairwise cluster marker
########################################################

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


clusters <- unique(as.character(Idents(sc)))
pairwise <- combn(clusters, 2)


#0,6,9
#1,2, 14,21
#12,17,20
#5,7,11
#3,15

pairwise_sel<-cbind( c(0,6),c(0,9),c(6,9),c(1,2),c(2,14),c(14,21),c(1,21),c(2,21),c(12,17),c(17,20),c(12,20),c(5,7),c(5,11),c(7,11),c(3,15),c(4,16),c(4,2),c(6,13),c(23,18),c(22,19),c(8,23))
pairwise_sel<-apply(pairwise_sel,c(1,2),as.character)



which_pair<-function(pairwise_sel,pairwise)
{
  index_sel<-c()
  index_sel_sel<-c()
  
  for (j in 1: ncol(pairwise_sel))
  {
    for (i in 1:ncol(pairwise)){
      #if (length(intersect(pairwise[,i], pairwise_sel[,j])) ==2) 
      if (  all(pairwise[,i] %in%  pairwise_sel[,j])  )
      {
        index_sel<-c(index_sel,i)
        index_sel_sel<-c(index_sel_sel,j)
      }
        
    }
  }
  
  return(cbind(index_sel,index_sel_sel))
}

res<-which_pair(pairwise_sel,pairwise)
index_sel<-res[,1]
index_sel_sel<-res[,2]
outFolder="5_PlotPairwiseClustersDGE/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 


m2 = read_rds("5_harmonyClustersDGE/pairwise_markers_results_df.rds")


m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
  group_by(symbol) %>%
  # mutate(H1=log2(length(clustr1)), H2=log2(length(clustr2))) %>%
  # filter(H1<=1 & H2<=1) %>%
  ungroup()




top20 <- m3 %>% group_by(clustr1,clustr2) %>% top_n(n = 20, wt = avg_log2FC) %>% ungroup()
#table(top20$clustr1, top20$clustr2)


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Rep","SNG.BEST.GUESS"))

myscale = 1/colSums(sc@assays$RNA@counts)*1000000



for(ii in 1:length(index_sel)){
  ##   
  
    
  clusters<-pairwise[,index_sel[ii]]
  
  genesel <- filter(top20,clustr1==clusters[1], clustr2==clusters[2])
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
  
  clusterss<-pairwise_sel[,index_sel_sel[ii]]
  
  fname=paste0(outFolder,"clusters_",clusterss[1],"_",clusterss[2],"_umap.png");
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

#top20 <- m3 %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC) %>% ungroup()
table(top20$clustr1,top20$clustr2)


for(ii in 1:length(index_sel)){
  ##   
  
  
  clusters<-pairwise[,index_sel[ii]]
  ##    
  genesel <- filter(top20,clustr1==clusters[1], clustr2==clusters[2])
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
  #fname=paste0(outFolder,"cluster_",i,"_dotplot2.png");
  
  clusterss<-pairwise_sel[,index_sel_sel[ii]]
  
  fname=paste0(outFolder,"clusters_",clusterss[1],"_",clusterss[2],"_dotplot2.png");
  
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






##table(sc@meta.data$seurat_clusters,sc@meta.data$SNG.BEST.GUESS)

################################
# heatmap
#https://github.com/satijalab/seurat/issues/287

scale.data<-sc@assays$RNA@scale.data
rw<-anno$gene_name[which(anno$kbid %in% rownames(scale.data))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(scale.data))]
scale.data<-scale.data[names(rw),]
rownames(scale.data)<-rw
sc@assays$RNA@scale.data<-scale.data
#colnames(scale.data)<-newsc$seurat_clusters[colnames(scale.data)]

m2 = read_rds("5_harmonyClustersDGE/pairwise_markers_results_df.rds")


m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>4) %>%
  group_by(symbol) %>%
  # mutate(H1=log2(length(clustr1)), H2=log2(length(clustr2))) %>%
  # filter(H1<=1 & H2<=1) %>%
  ungroup()


top10 <- m3 %>% group_by(clustr1,clustr2) %>% top_n(n = 5, wt = avg_log2FC) %>% ungroup()

fname=paste0(outFolder,"Heatmap-top10.pdf");
pdf(fname,width=40,height=20)
#DoHeatmap(newsc, features = top10$symbol) + NoLegend()
DoHeatmap(sc, features = unique(top10$symbol))+ #+ NoLegend()
theme(text = element_text(size = 10),axis.text.y = element_text(size = 10)) + NoLegend()#
dev.off()

