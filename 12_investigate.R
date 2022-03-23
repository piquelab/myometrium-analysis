library(Seurat)
library(Matrix)
library(tidyverse)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



outFolder="./12_investigate/"
system(paste0("mkdir -p ", outFolder))



anno <- read_rds("3_MergeDemux_Output/anno.rds")
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")


aa <- FetchData(sc,c("UMAP_1","UMAP_2","Location","Condition","Origin","status","seurat_clusters","cluster_name","Library","Pregnancy_ID")) 
head(aa)



myscale = 1/colSums(sc@assays$RNA@counts)*1000000
genesym=c("AGTR2")

genesel <- filter(anno,gene_name== "AGTR2")

genexpr <- map_dfr(1:nrow(genesel),function(i){
  aux <- FetchData(sc,genesel$kbid[i])
  colnames(aux) <- "expr"
  aux <- cbind(aa,aux*myscale)
  aux$gene=genesel$kbid[i]
  aux$symbol=genesel$gene_name[i]
  aux
})

fname=paste0(outFolder,"AGTR2_umap.png");
##png(fname,width=1000,height=800)
p1 <- genexpr %>%  arrange(symbol,expr) %>%
  ggplot(aes(UMAP_1,UMAP_2,color=log10(0.1+expr))) +
  geom_point(size=0.1) +
  ##        scale_fill_viridis_c(option = "plasma") +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  facet_wrap(~symbol) +
  theme_bw()
ggsave(fname,p1,width=10,height=8)