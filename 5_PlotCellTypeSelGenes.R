## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
##library(SingleR)



anno <- read_rds("3_MergeDemux_Output/anno.rds")
#sc <- read_rds("4_harmony/sc.NormByLocationRep.Harmony.rds")
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

#load("3_MergeDemux_Output/annoOnly.Rdata")




outFolder="./5_PlotClustersSelGenes/"
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



##for(i in 0:max(top20$cluster)){
    ##    

##i <- "Endothelial" 
##genesym=c("CD74","SELP","THBD","PECAM1","CD93","FLT1","CCL21")


i <- "SMC" 
genesym=c("ACTG2","MYH11","ACTA2","MYO1B","MYLK")


genesel <- filter(anno,gene_name %in% genesym)

##        FetchData(sc,genesel$gene)
    genexpr <- map_dfr(1:nrow(genesel),function(i){
        aux <- FetchData(sc,genesel$kbid[i])
        colnames(aux) <- "expr"
        aux <- cbind(aa,aux*myscale)
        aux$gene=genesel$kbid[i]
        aux$symbol=genesel$gene_name[i]
        aux
    })
##    
    cat(i,dim(genexpr),"\n")
    
  ##
    fname=paste0(outFolder,"Ref_",i,"_umap.png");
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


##    sum_rec <- genexpr %>% dplyr::group_by(symbol,Location,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))
    sum_rec <- genexpr %>% dplyr::group_by(symbol,seurat_clusters) %>% dplyr::summarize(Prop=mean(expr>0),Expr=mean(expr[expr>0]))

##
    fname=paste0(outFolder,"Ref_",i,"_dotplot2.png");
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
    ggsave(fname,p0,width=8,height=4.1)
    cat(fname,"\n")



## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.V2.pdf");
pdf(fname,width=10,height=2)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Rep","SNG.BEST.GUESS")) 
aa$seurat_clusters <- clust2name[aa$seurat_clusters]
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



### END- HERE ###
########################################################


