
#############################################################
###  Heatmap plot showing the comparison between the clusters based of log2fc
### 
#############################################################

library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)

outFolder="6_celltype_cor_heatmap_plot/"
system(paste0("mkdir -p ", outFolder))
#########################


res_de <-read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
de_gene<-unique(res_de$gene_name)

res <-read_tsv("7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
res$cluster<-unlist(strsplit(res$cname,"_"))[seq(1,2*length(res$cname),by=2)]
res<-res  %>% filter ( gene_name %in% de_gene)

param<-names(table(res$cluster))[table(res$cluster)>20]

clusters<-unique(param)
cluster2name<-read.csv("clust2Name.csv",stringsAsFactors = FALSE)
colnames(cluster2name)<-c("id","name")
clusternames<-paste0(cluster2name$id,"_",cluster2name$name)
names(clusternames)<-as.character(cluster2name$id)



cor_matrix<-matrix(NA,nrow=length(clusters),ncol=length(clusters))


for (i in 1:length(clusters))
{
  for (j in 1: length(clusters))
  {
      resi<-res %>% dplyr::filter(cluster == clusters[i]) %>% dplyr::select(gene_name,log2FoldChange,lfcSE)
      resj<-res %>% dplyr::filter(cluster == clusters[j])%>% dplyr::select(gene_name,log2FoldChange,lfcSE)
      colnames(resj)<-c("gene_name","log2FoldChange2","lfcSE2")
      res_intersect<-resi%>% inner_join(resj)
      res_intersect<-res_intersect%>% dplyr::select(log2FoldChange,log2FoldChange2)
      if(nrow(res_intersect)>5)
      {
        cr<-cor(res_intersect)[1,2]
        cor_matrix[i,j]<-cr
        cor_matrix[j,i]<-cr
      }
      }
  }
 
rownames(cor_matrix)<-clusternames[clusters]
colnames(cor_matrix)<-clusternames[clusters]

melted_cormat <- melt(cor_matrix)
head(melted_cormat)
colnames(melted_cormat)[c(1,2)]<-c("x","y")




###############
fname=paste0(outFolder,"heatmap_celltype_cor.png");
png(fname,width=1600,height=1200)
ggplot(data = melted_cormat, aes(x=x, y=y, fill=value),xlab="") + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=25),legend.text=element_text(size=20))+
  scale_fill_gradient2(low = "darkblue", high = "red",mid="white",na.value="black", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")
dev.off()
