#"/wsu/home/groups/prbgenomics/firs"
library(igraph)
library(readxl)
library(tidyverse)
library(factoextra)
library(ggfortify)
library(Rtsne)
require(data.table)
library(umap)
library(mclust)
library(clusterProfiler)
library(annotables)
library(Matrix)
outFolder<-"./14_PCA_plots/"

outFolder<-"./14_PCA_plots_signaturesCombined/"
system(paste0("mkdir -p ",outFolder))


##### single cell signature 
m2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
  group_by(gene) %>%
  mutate(H=log2(length(cluster))) %>%
  filter(H<=1) %>%
  ungroup()

table(m3$cluster)
dat <- m3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% ungroup()
table(dat$cluster)


clust2Names<-c("Stromal","Macrophage","Macrophage","Endothelial","Monocyte","T-cell","Decidual","T-cell","LED","Stromal","ILC","NK-cell","Smooth muscle cells","Myofibroblast","Macrophage","Endothelial","DC","Smooth muscle cells","EVT","Plasmablast","Smooth muscle cells","Macrophage","B-cell","Unciliated Epithelial")

#clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
#clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)
dat$cluster<-clust2Names[as.character(dat$cluster)]

dat<-dat %>% select(cluster,gene,symbol,pct.1,pct.2,avg_log2FC,p_val,p_val_adj)
colnames(dat)<-c("cluster" ,  "ENSG" ,     "gene"   ,   "pct.1"   ,  "pct.2"    , "avg_logFC" ,"p_val"  ,   "p_val_adj")


generate_gene_cell_matrix<-function(dat) 
{
      dat<-as.data.frame(dat)
      edge_list<-dat[,c("gene","cluster","avg_logFC")]
      
      edge_list<-as.data.frame(edge_list)
      colnames(edge_list)[3]<-"weight"
      edge_list<-edge_list[!duplicated(edge_list), ]
      edge_list$weight<-1
      cl<-unique(as.character(edge_list[,"cluster"]))
      rw<-unique(as.character(edge_list[,"gene"]))
      edge_list$cluster<-as.character(edge_list$cluster)
      gene_cell_matrix<-matrix(0,nrow=length(rw),ncol = length(cl))
      rownames(gene_cell_matrix)<-rw
      colnames(gene_cell_matrix)<-cl
      gene_cell_matrix[as.matrix(edge_list[,1:2])] <- edge_list[,3]
      return(gene_cell_matrix)

}



### Bulk data 

load("TLTNLmyoToCaseWest.rdata")
eset_bulk<-eset
colnames(eset_bulk)<-unlist(strsplit(colnames(eset),"_"))[seq(2,2*length(colnames(eset)),by=2)]
sample_names<-colnames(eset_bulk)
rownames(eset_bulk)<-unlist(strsplit(rownames(eset),"_"))[seq(2,2*length(rownames(eset)),by=2)]

tl<-table(rownames(eset_bulk))
nonrepeat<-eset_bulk[rownames(eset_bulk) %in% names(tl)[tl==1],]
repeats<-sapply(names(tl)[tl>1], function(x) {
  print(x)
  mat<-eset_bulk[which(rownames(eset_bulk)==x),]
  mat[which(rowSums(mat)==max(rowSums(mat))),]
  
})
repeats<-t(repeats)
rownames(repeats)<-names(tl)[tl>1]
bulk_matrix<-rbind(nonrepeat,repeats)
Rcount<-bulk_matrix



gene_cell_matrix<-generate_gene_cell_matrix(dat) #gene = "ENSG",clustertype = "cluster",

rw1<-rownames(Rcount)
rw2<-rownames(gene_cell_matrix)
rw<-intersect(rw1,rw2)
gene_cell_matrix<-gene_cell_matrix[rw,]
Rcount<-Rcount[rw,]
gene_cell_matrix<-t(gene_cell_matrix)  # cell types * gene



# Scaling/Normalization
#log normalized
# Generating metagene matrix 

# // non- corrected data
#myscale = 1/colSums(Rcount)
#Rcount=Rcount*myscale



Metgene<-gene_cell_matrix  %*% Rcount
Metgene<-t(Metgene)



#log normalized // non- corrected data

#Metgene <- log10(as.matrix(Metgene)+1)


lg<-c("TNL" ,"TL")
sg<-colnames(Rcount)
names(sg)<-sample_names
sgCol<-rep("red",length(sg))
names(sgCol)<-rownames(Metgene)
sgCol[which(sg[names(sgCol)]=="TNL")]<-"#333399"
sgCol[which(sg[names(sgCol)]=="TL")]<-"#A50021"
sg<-sg[names(sgCol)]

res.sample.pca <- prcomp(t(Rcount))
pdf(paste0(outFolder,"Rcount_pca_1_2_plot.pdf"))
plot(res.sample.pca$x[,1], res.sample.pca$x[,2], col=sgCol, pch=16, main='PCA',cex=1.5,cex.axis=1,font=2)
legend("bottomright", pch = 20, col=c("#333399" , "#A50021") , legend=lg, bty='n', cex=.9)
dev.off()


pdf(paste0(outFolder,"Rcount_pca_plot_2_3.pdf"))
plot(res.sample.pca$x[,2], res.sample.pca$x[,3], col=sgCol, pch=16, main='PCA',,cex=1.5,cex.axis=1,font=2)
legend("bottomleft", pch = 20, col=col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=.9)
dev.off()


pdf(paste0(outFolder,"Rcount_pca_plot_3_4.pdf"))
plot(res.sample.pca$x[,3], res.sample.pca$x[,4], col=sgCol, pch=16, main='PCA',cex=1.5,cex.axis=1,font=2)
legend("bottomleft", pch = 20, col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=.9)
dev.off()
#Metgene<-scale(Metgene,center = TRUE)

pdf(paste0(outFolder,"Rcount_pca_plot_2_4.pdf"))
plot(res.sample.pca$x[,2], res.sample.pca$x[,4], col=sgCol, pch=16, main='PCA',cex=1.5,cex.axis=1,font=2)
legend("bottomleft", pch = 20, col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=.9)
dev.off()


library(rgl)



#Compute PCA
res.pca <- prcomp(Metgene)

# eig.val <- get_eigenvalue(res.pca)
# eig.val

# for coloring samples

var<-get_pca_var(res.pca)
head(var$contrib)
which(var$contrib[,3]==max(var$contrib[,3]))
which(var$contrib[,4]==max(var$contrib[,4]))


pdf(paste0(outFolder,"pca1_2_plot.pdf"))
par(mar=c(5,6,4,1)+.1)
plot(res.pca$x[,1], res.pca$x[,2], col=sgCol, pch=16, main='Principal component analysis',cex=2,cex.axis=1,font=2,cex.lab = 1.5,xlab="PC1", ylab="PC2")
legend("topleft", pch = 20, col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=1.5)
dev.off()



pdf(paste0(outFolder,"pca2_3_plot.pdf"))
plot(res.pca$x[,2], res.pca$x[,3], col=sgCol, pch=16,main='Principal component analysis',cex=2,cex.axis=1,font=2,cex.lab = 1.5,xlab="PC2", ylab="PC3")
legend("bottomleft", pch = 20, col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=.75)
dev.off()

pdf(paste0(outFolder,"pca3_4_plot.pdf"))
plot(res.pca$x[,3], res.pca$x[,4], col=sgCol, pch=16,main='Principal component analysis',cex=2,cex.axis=1,font=2,cex.lab = 1.5,xlab="PC3", ylab="PC4")
legend("bottomleft", pch = 20, col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=.75)
dev.off()


pdf(paste0(outFolder,"pca2_4_plot.pdf"))
plot(res.pca$x[,2], res.pca$x[,4], col=sgCol, pch=16,main='Principal component analysis',cex=2,cex.axis=1,font=2,cex.lab = 1.5,xlab="PC2", ylab="PC4")
legend("bottomleft", pch = 20, col=c("#333399" , "#A50021"), legend=lg, bty='n', cex=.75)
dev.off()





#Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.

pdf(paste0(outFolder,"fviz_pca_dim1_2_variables.pdf"))
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()


pdf(paste0(outFolder,"fviz_pca_dim2_3_variables.pdf"))
fviz_pca_var(res.pca,axes = c(2, 3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf(paste0(outFolder,"fviz_pca_dim3_4_variables.pdf"))
fviz_pca_var(res.pca,axes = c(3, 4),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()


pdf(paste0(outFolder,"fviz_pca_dim1-3_variables.pdf"))
fviz_pca_var(res.pca,axes = c(1, 3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf(paste0(outFolder,"fviz_pca_dim2-4_variables.pdf"))
fviz_pca_var(res.pca,axes = c(2, 4),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()


pdf(paste0(outFolder,"umap.pdf"),onefile=FALSE)
umap_out <- umap(Metgene,n_neighbors=5)
plot(umap_out$layout, col=sgCol, pch=16, main="UMAP")
dev.off()























##### [PC1:PC4, samples]######

##################### t-test #####################



comparing_samples<-function(samples,res.pca,select1,select2,PC=1)
{
  data_pca<-t(as.matrix(res.pca$x[,PC])) # PC1, PC2, PC3
  x<-data_pca[,which(colnames(data_pca)==select1)]
  y<-data_pca[,which(colnames(data_pca)==select2)]
  res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE,conf.level = 0.95)
  return(res)
}

result_PCs<-array() #matrix(nrow=3,ncol=3)
for (PC in 1:4)
{
  comparison_TL_TNL<-comparing_samples(samples,res.pca,select1="TL" ,select2="TNL",PC)
  pvalue<-comparison_TL_TNL$p.value
  
  result_PCs<-cbind(result_PCs,pvalue)
}

result_PCs<-result_PCs[,-1]
#rownames(result_PCs)<-c("TIL vs TNL")
names(result_PCs)<-c("PC1","PC2","PC3","PC4")

result_PCs_adjusted<-matrix(p.adjust(result_PCs,"bonferroni"),nrow=1, ncol=4)
rownames(result_PCs_adjusted)<-c("TIL vs TNL")
colnames(result_PCs_adjusted)<-c("PC1","PC2","PC3","PC4")

data<-rbind(result_PCs_adjusted,result_PCs)

write.csv(data,file=paste0(outFolder,"result_PCs_TIL_TNL_comparison.csv"))


###################

resJoin<-cbind(res.pca$x[,2], res.pca$x[,4])
colnames(resJoin)<-c("PC2","PC4")
resJoin<-as.data.frame(resJoin)
p2 <- resJoin %>% ggplot(aes(PC2,PC4)) +
  geom_point(color="black")+ #aes(colour = ref_color)) +
  geom_smooth(method=lm, se=FALSE,linetype = "dashed", color="black")+
  theme_bw()
fname=paste0(outFolder,paste0("scatterplotPC2_4.png"))
ggsave(fname,p2,width=6,height=4.5)

