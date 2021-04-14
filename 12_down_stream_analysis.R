library(UpSetR)
library(tidyverse)
library(ComplexHeatmap)
library(mashr)
library(rmeta)
library(ashr)
library(mashr)
library(pheatmap)


outFolder <- paste0("12_downstream_analysis/")
system(paste0("mkdir -p ",outFolder))


res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-clust2Names

res$cluster_colors<-cluster.Colors[res$Cell_type]

res$Cell_type<-clust2Names[res$Cell_type]



################################################ 
### upset plot
################################################
res<-res %>% filter(abs(log2FoldChange)>0 & padj< 0.1) #%>% select(gene_name,Cell_type)

res <- res%>% mutate (DE=ifelse(abs(log2FoldChange)>0 & padj< 0.1 ,1,0))

res_df<-matrix(0,nrow=length(unique(res$gene_name)),ncol=length(unique(res$Cell_type)))
colnames(res_df)<-unique(res$Cell_type)
rownames(res_df)<-unique(res$gene_name)
cell_types<-unique(res$Cell_type)
for (i in 1:length(cell_types))
{
 res_celltype<-res %>% filter(Cell_type %in% cell_types[i]) %>%select(gene_name)
 gn<-res_celltype$gene_name
 res_df[gn,cell_types[i]]<-1
}
res_df<-as.data.frame(res_df)
mt = make_comb_mat(res_df)#, top_n_sets = 25)
m2 <- mt[,comb_size(mt)>20]
#fig0 <- UpSet(m2, set_order=colnames(res_df),comb_order=order(-comb_size(m2)), row_names_side="left",right_annotation=NULL)#sets.bar.color = cluster.Colors[colnames(res_df)]
fig0 <- UpSet(m2, set_order=colnames(res_df),comb_order=order(comb_size(m2)))#sets.bar.color = cluster.Colors[colnames(res_df)]
pdf(file=paste0(outFolder,"upsetplot_comb_size25.pdf"))#,width=45,height=20)# or other device
fig0
dev.off()



##############

set.seed(1)



########################################################################
#load data 
########################################################################
#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")

res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-clust2Names

res$cluster_colors<-cluster.Colors[res$Cell_type]

res$Cell_type<-clust2Names[res$Cell_type]


########################################################################
#### mash analysis 
########################################################################

hist(table(res$kbid))
max(table(res$kbid))
sum(table(res$kbid)==16)
selgenes <- names(which(table(res$kbid)==16))

Bhat<-res %>% filter(kbid %in% selgenes) %>% select(kbid,Cell_type,log2FoldChange) %>%
  pivot_wider(names_from = Cell_type, values_from = log2FoldChange)
rw<-Bhat$kbid
Bhat<-as.matrix(Bhat)
rownames(Bhat)<-Bhat[,"kbid" ]
Bhat<-Bhat[,-1]
Bhat<-apply(Bhat,c(1,2),as.numeric)

Shat<-res %>% filter(kbid %in% selgenes) %>% select(kbid,Cell_type,lfcSE) %>%
  pivot_wider(names_from = Cell_type, values_from = lfcSE)

rw<-Shat$kbid
Shat<-as.matrix(Shat)
rownames(Shat)<-Shat[,"kbid" ]
Shat<-Shat[,-1]
Shat<-apply(Shat,c(1,2),as.numeric)

########################################################################
# data
########################################################################
data = mash_set_data(Bhat, Shat)


########################################################################
#covariance matrices
########################################################################
U.c = cov_canonical(data)
print(names(U.c))


########################################################################
#fit the model
########################################################################
m.c = mash(data, U.c)


########################################################################
#Extract Posterior Summaries
########################################################################

#local false sign rate

lfsr_m.c<-get_lfsr(m.c)

head(get_lfsr(m.c))

#######################


#############################################
#meta plot 
#############################################

# selecting genes for meta plot 

numsig.rel <- rowSums(lfsr_m.c<0.9)
mytop <- (numsig.str==1 & numsig.rel==1)
sum(mytop)


nabstr <- rowSums(lfsr_m.c<0.01)
nberel <- rowSums(lfsr_m.c>0.2)
mytop <- (nabstr<3 & nabstr>1 & (nberel+nabstr>10))
sum(mytop)

mytop <- (nabstr<3 & nabstr>1 & (nberel+nabstr>8))
sum(mytop)

#mash_plot_meta(m.c,get_significant_results(m.c)[1])

nabstr <- rowSums(lfsr_m.c<0.01)
nberel <- rowSums(lfsr_m.c>0.2)
mytop <- (nabstr<4 & nabstr>1 & (nberel+nabstr>12))
sum(mytop)



outFolder2<-"./8_outputs_DESeq_Plots/RNASeq/"
outFolder2<-"./8_outputs_DESeq_Plots/DriverMap/"
intersected_genes<-read.csv(paste0(outFolder2,"intersected_genes.csv"),stringsAsFactors = FALSE)

#outFolder<-"12_downstream_analysis/DriverMap/"
outFolder<-"12_downstream_analysis/RNASeq/"

mytop<-rep(FALSE,length(rownames(lfsr_m.c)))
names(mytop)<-rownames(lfsr_m.c)

#some specific genes 
sample_genes<-res$kbid[which(res$gene_name %in% hub_scores$gene_name[1:50])]
#sample_genes<-res$kbid[which(res$gene_name %in% c("OXTR","COX2","CAM","Cx43","GJA1","PTGES2")) ]
sample_genes<-res$kbid[which(res$gene_name %in% c("TCIRG1","ADAM15","MAP1B","PNPLA2","PLEKHO1","IFI16","HTRA1","FAM129A","CD52","EFHD2","CD82","ACOT7","WDR1","RRP9","OXTR","COX2","CAM","Cx43","GJA1","PTGES2","SCD")) ]

sample_genes<-res$kbid[which(res$gene_name %in% c("MMP9","STOM","METTL9","MSR1","MCEMP1","TFDP1","TSTA3","GMPR","HK1","FOXO3","P2RX7","MYOF","ABCG1","OLR1")) ]

sample_genes<-res$kbid[which(res$gene_name %in% unique(intersected_genes$gene_name)) ]
mytop[which(names(mytop) %in%sample_genes)]<-TRUE

# meta plot
plot_func<-function(m.c,mytop,k=1)
{
  for ( k in 1:length(which(mytop)))
  {
    par(mar=c(2, 1 ,4 ,3))
    i<-which(mytop)[k]
    plot.title<-res$gene_name[which(res$kbid==names(i))[1]]
    system(paste0("mkdir -p ",outFolder))
    fname=paste0(outFolder,plot.title,".pdf");
    pdf(fname,width=7,height=7)
    
    print(plot.title)
    
    metaplot(get_pm(m.c)[i,],get_psd(m.c)[i,],colors = meta.colors(box = as.character(cluster.Colors))
             ,xlim = c(-1,1),xlab = "",ylab = "")
    #legend("topleft",legend=as.character(unique(res$Cell_type)), fill=cluster.Colors[unique(res$Cell_type)],bty = "n",cex=0.8) #
    title(plot.title)
    dev.off()
  }
  
  
}


# head(get_pm(m.c))
# 
# #posteriore standard deviation
# head(get_psd(m.c))  
# 
# 
# # effects that are “significant”,
# head(get_significant_results(m.c))
# print(length(get_significant_results(m.c)))
# 
# 
# print(head(get_significant_results(m.c, conditions=1)))
# 
# #assess sharing of significant signals among each pair of conditions
# print(get_pairwise_sharing(m.c)) # heatmap
# 
# #For example, here by setting the factor to be 0 you assess only if they are the same sign:
# print(get_pairwise_sharing(m.c, factor=0))
# 
# #Measure of fit (log-likelihood)
# print(get_loglik(m.c))
# 
# print(get_estimated_pi(m.c))
# 
# 
# #barplot
# par(mar=c(14, 4 ,4 ,3))
# barplot(get_estimated_pi(m.c),las = 2)
# 
# #Metaplot
# mash_plot_meta(m.c,get_significant_results(m.c)[1])
# 
# 
# metaplot(get_pm(m.c)[4217,],get_psd(m.c)[4217,],colors = meta.colors(box = as.character(cluster.Colors))
#          ,xlim = c(-1,1),xlab = "",ylab = "")
#                      
# mash_plot_meta(m.c,which(mytop)[1])
####################################

# Bhat<-res %>% filter(kbid %in% selgenes) %>% select(kbid,Cell_type,log2FoldChange) %>%
#   pivot_wider(names_from = Cell_type, values_from = log2FoldChange)
# rw<-Bhat$kbid
# Bhat<-as.matrix(Bhat)
# rownames(Bhat)<-Bhat[,"kbid" ]
# Bhat<-Bhat[,-1]
# Bhat<-apply(Bhat,c(1,2),as.numeric)


# wider matrix generation

selgenes <- names(which(table(res$kbid)==16))
padj<-res %>% filter(!is.na(pvalue) & kbid %in% selgenes) %>% select(kbid,Cell_type,padj) %>%
  pivot_wider(names_from = Cell_type, values_from = padj)
rw<-padj$kbid
padj<-as.matrix(padj)
rownames(padj)<-padj[,"kbid" ]
padj<-padj[,-1]
padj<-apply(padj,c(1,2),as.numeric)
padj[padj>=0.1]<--1
padj[padj<=0.1 & padj>=0]<-1
padj[padj==-1]<-0
padj<-na.omit(padj)
DEgene_cluster<-padj
hist(rowSums(DEgene_cluster))


#############################################
## plots from result of mash analysis
#############################################

## heatmap 



pairwise_sharing_data<-get_pairwise_sharing(m.c)

fname=paste0(outFolder,"heatmap_pairwise_sharing.pdf");
pdf(fname,width=7,height=7)
pheatmap(pairwise_sharing_data,cluster_rows=TRUE,cluster_cols=TRUE,scale="none")
dev.off()


#############################################
### upset plot 
#############################################

lfsr<- get_lfsr(m.c)
lfsr[which(lfsr<0.1)]<--1
lfsr[which(lfsr>=0.1)]<-0
lfsr[which(lfsr==-1)]<-1

lfsr<-as.data.frame(lfsr)
mt = make_comb_mat(lfsr,complement_size = 0)#, top_n_sets = 25)

m2 <- mt[,comb_size(mt)>10]
#m2 <- mt[comb_degree(mt)>0]
#fig0 <- UpSet(m2, set_order=colnames(lfsr),comb_order=order(comb_size(m2)))#sets.bar.color = cluster.Colors[colnames(res_df)]
fig0 <- UpSet(m2,set_order=colnames(lfsr),comb_order=order(comb_size(m2)))
pdf(file=paste0(outFolder,"upsetplot_lfsr_comb_size25.pdf"))#,width=45,height=20)# or other device

fig0
dev.off()



