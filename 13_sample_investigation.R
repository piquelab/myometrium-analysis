library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)


library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ReactomePA)
#################
##library(SingleR)




future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)



outFolder="13_sample_investigation_plots/"
system(paste0("mkdir -p ", outFolder))


###########################################
clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
              "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast",
              "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

#clust2Name<-paste0(c(0:23),"_",clust2Name)
names(clust2Name)<-c(0:23)

sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")


md <- read_rds("./4_harmony_cellClass_PBMC/sc.NormByLocation.ref.Anchors.rds") %>%
  as.data.frame %>%
  rownames_to_column("BARCODES") %>%
  select(BARCODES,scLabor_ID=predicted.celltype.l2,scLabor_Score=predicted.celltype.l2.score)

md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
  left_join(md) 

identical(md$BARCODES,rownames(sc@meta.data))



sc@meta.data$cluster_name <- clust2Name[sc@meta.data$seurat_clusters]


aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status","SNG.BEST.GUESS","Pregnancy_ID","cluster_name")) 


cc <- read.csv("ParturitionAndMyoSamples.csv",stringsAsFactors = FALSE) 


cc %>% filter (Pregnancy_ID %in% c("MRN2762","MRN2215","MRN8739","MRN6897","HPL20914") ) %>% select(Pregnancy_ID,Sample_ID,Condition) 

#cc %>% filter (extract_numeric(Pregnancy_ID) %in% c("2762","2215","8739","6897","20914") ) %>% select(Pregnancy_ID,Sample_ID,Condition) 


head(cc)

cc<- cc %>% mutate(Traverse= ifelse( Pregnancy_ID %in% c("MRN2762","MRN2215","MRN8739","MRN6897","HPL20914"), "High", "Low"))
cc<-cc  %>%  select(Pregnancy_ID,Traverse)
aa <- aa %>% inner_join(cc)
aa$cluster_name <- clust2Name[as.character(aa$seurat_clusters)]




############ bar plot ####################################



data<-aa %>% group_by(Pregnancy_ID,cluster_name,Traverse) %>% summarize(n=n()) %>% 
  group_by(Traverse)%>% mutate(prop=n/sum(n)) 

fname=paste0(outFolder,"Barplot_myo_proportion_v2.pdf");
pdf(fname,width=30,height=10)
p2 <- ggplot(data,aes(x=Pregnancy_ID,y=prop, fill=Traverse)) +
  geom_bar(stat="identity",position = position_dodge()) +
  guides(fill = guide_legend(override.aes = list(size=5),title="Transverse incision")) +
  facet_grid(.~cluster_name) +
  coord_flip() +
  scale_fill_manual(values=c("Low"="#333399","High"="#A50021"))+
  xlab("")+
  #theme(legend.text=element_text(size=30), axis.text=element_text(size=30,angle = 45), axis.title=element_text(size=30,face="bold"))+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  theme_bw()
p2
##    theme_black()
dev.off()



  
  
  
data<-aa%>% group_by(cluster_name)%>% summarise(total_count = n(),Low=count(Traverse=="Low") ,High=count(Traverse=="High")) 
data$Low<-data$Low/data$total_count
data$High<-data$High/data$total_count
data<-data %>% select(cluster_name,Low,High)

library(reshape2)
library(ggplot2)
dat2 <- melt(data)
colnames(dat2)<-c("cluster","proportion","value")
dat2$clusternumber<-names(clust2Name)[which(clust2Name %in% as.character(dat2$cluster))]
dat2$clusternumber<-as.numeric(dat2$clusternumber)
# dat2<-dat2[order(dat2$clusternumber,decreasing = FALSE),]


fname=paste0(outFolder,"Barplot_myo_proportion.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(dat2,aes(x=reorder(cluster,-clusternumber),y=value, fill=proportion)) +
  geom_bar(stat="identity",position = position_dodge()) +
  guides(fill = guide_legend(override.aes = list(size=5),title="Transverse incision")) +
  #facet_grid(.~Location) +
  coord_flip() +
  scale_fill_manual(values=c("Low"="#333399","High"="#A50021"))+
  xlab("")+
  theme(legend.text=element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=20,face="bold"))+
theme_bw()
p2
##    theme_black()
dev.off()



aa<-aa[order(as.numeric(aa$seurat_clusters),decreasing =TRUE),]
aa$seurat_clusters <- factor(aa$seurat_clusters,levels=unique(aa$seurat_clusters))
aa$cluster_name <- factor(aa$cluster_name,levels=unique(aa$cluster_name))
fname=paste0(outFolder,"Barplot_myo.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=reorder(cluster_name,-seurat_clusters),fill=Traverse)) +
  geom_bar(position = position_dodge()) +
  guides(fill = guide_legend(override.aes = list(size=5),title="Transverse incision")) +
  facet_grid(.~Location) + coord_flip() +
  scale_fill_manual(values=c("Low"="#333399","High"="#A50021"))+
  xlab("")+
  theme(legend.text=element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=20,face="bold"))
  theme_bw()
p2
##    theme_black()
dev.off()



######### UMAP plot
fname=paste0(outFolder,"UMAP_myo.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Traverse)) +
  geom_point(size=0.1) +
  ##    scale_color_manual(values=group.colors) +
  guides(colour = guide_legend(override.aes = list(size=10),title="Transverse incision")) +
  scale_color_manual(values=c("Low"="#333399","High"="#A50021"))+
  theme(legend.text=element_text(size=20), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
##    facet_wrap(~LocTime) +
theme_bw()
p1
##    theme_black()
dev.off()

################################################################################
# cell count comparison // low vs high 

################################################################################
aa$seurat_clusters <- clust2Name[aa$seurat_clusters]

cc <- aa %>% group_by(seurat_clusters,Traverse,Pregnancy_ID) %>%
  summarize(n=n()) %>%
  group_by(seurat_clusters) %>% mutate(p0=sum(n)/nrow(aa)) %>%
  group_by(Pregnancy_ID) %>%
  mutate(nt=sum(n),p=n/nt,z=(p-p0)/sqrt(p0*(1-p0)*(1/nt+1/nrow(aa))))


ccdiff <-  cc %>% group_by(seurat_clusters) %>% summarize(wilcox.pval = wilcox.test(p ~ Traverse)$p.value,
                                                          test.t = t.test(z ~ Traverse)$statistic) %>% ungroup() %>%
  mutate(wilcox.padj = p.adjust(wilcox.pval))  

#ccdiff_filtered<-ccdiff %>% filter(wilcox.padj<0.1)

write_tsv(ccdiff,file=paste0(outFolder,"difference_wilcox_cell_count_traverse.tsv"))



################################################################################################################################################################################################################################################

#                                                        Enrichment analysis
################################################################################################################################################################


##############################################################################
# SMC-1
################################################################################

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
## cluster colors 

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
res$Cell_type<-clust2Names[res$Cell_type]

#ENTREZID
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

res<-res %>% filter(Cell_type=="Smooth muscle cells-1")


geneList <- -log10(res$pvalue)
names(geneList) <- res$ENTREZID
geneList = sort(geneList, decreasing = TRUE)


## Go term
gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
#load(paste0(outFolder,"gseGO.res_SMC-1.RData"))
save(gseGO.res,file=paste0(outFolder,"gseGO.res_SMC-1.RData"))

res_df_gseGO<-gseGO.res@result%>% filter(qvalues<=0.1)

#pdf(paste0(outFolder,"gseGO_SMC_DotPlot.pdf"),width=10,height=10)

fname<-paste0(outFolder,"gseGO_SMC_DotPlot.png")
p2<-ggplot(res_df_gseGO, 
       aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  theme_bw()
  #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
#dev.off()
ggsave(fname,p2,width=6,height=5)


fname<-paste0(outFolder,"gsea-enrichmentscore-plot.png")
p2 <- gseaplot2(gseGO.res, geneSetID = 1:6, subplots = 1)
ggsave(fname,p2,width=12,height=10)


#http://www.informatics.jax.org/go/report.txt?goID=GO:0060047
GO0060047_heart_contraction<-read_tsv("13_sample_investigation_plots/GO_term_summary_20210510_172542.txt")
colnames(GO0060047_heart_contraction)[2]<-"gene_name"
GO0060047_heart_contraction <- GO0060047_heart_contraction %>% left_join(eg) %>% filter(!is.na(ENTREZID))
GO0060047_heart_contraction$ENTREZID[which( GO0060047_heart_contraction$ENTREZID %in% names(geneList))]


core_genes<-gseGO.res@result %>% filter(ID=="GO:0060047") %>% select(core_enrichment) %>% unlist
core_genes<-unlist(strsplit(core_genes,"/"))

length(res %>% filter ( ENTREZID %in% core_genes & padj<0.1) %>% select(ENTREZID) %>% unlist)

res %>% filter ( ENTREZID %in% core_genes & padj<0.1) %>% select(gene_name) %>% unlist




###### Reactome
message("gsePathway")
gseRPath.res <- gsePathway(geneList,pvalueCutoff = 1)
print(head(gseRPath.res))
save(gseRPath.res,file=paste0(outFolder,"gsePathway_SMC-1.RData"))

res_df<-gseRPath.res@result
which(res_df$ID =="R-HSA-445355")

core_genes<-res_df %>% filter(ID=="R-HSA-445355") %>% select(core_enrichment) %>% unlist
core_genes<-unlist(strsplit(core_genes,"/"))
length(res %>% filter ( ENTREZID %in% core_genes & padj<0.1) %>% select(ENTREZID) %>% unlist)






res_df<-res_df %>% filter(qvalues<=0.1) #[1:10,]

fname<-paste0(outFolder,"/gsePathway_SMC_DotPlot.png")
#pdf(paste0(outFolder_cname_plots,"gsePathway_SMC_DotPlot.pdf"),width=10,height=10)
p3<-ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  #theme_bw(base_size = 11) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
  ylab(NULL) +
  theme_bw()
  #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
#dev.off()
ggsave(fname,p3,width=8,height=5)



gseKEGG.res <-gseKEGG( geneList)
print(head(gseKEGG.res))
save(gseKEGG.res,file=paste0(outFolder,"gseKEGG.res_SMC-1.RData"))

res_df<-gseKEGG.res@result %>% filter(qvalues<=0.1) #[1:10,]

fname<-paste0(outFolder,"gseKEGG.res_SMC_DotPlot.png")
#pdf(paste0(outFolder_cname_plots,"gseKEGG.res_SMC_DotPlot.pdf"),width=10,height=10)
p1<-ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  #theme_bw(base_size = 11) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  theme_bw()
  #theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
#dev.off()
ggsave(fname,p1,width=8,height=5)



##########################################################
# ORA  
##########################################################
res<-res %>% filter(!is.na(padj) & !is.na(log2FoldChange)) 
#genes <- filter(res,padj<0.1,abs(log2FoldChange)>=0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique

genes <- filter(res,padj<0.1) %>% unlist %>% unique


geneUniv <- res %>% dplyr::select(ENTREZID) %>% unlist %>% unique



ego<- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
res_df_enrichGO<-ego@result
res_df_enrichGO<-res_df_enrichGO%>% filter(p.adjust<0.1)
res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichGO_SMC-1_DotPlot.pdf"),width=10,height=4)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()



enrichKEGG.res <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")

res_df_enrichKEGG<-enrichKEGG.res@result

res_df_enrichKEGG<-res_df_enrichKEGG%>% filter(qvalue<0.1)


res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})
# 
# pdf(paste0(outFolder,"enrichKEGG_SMC-1_DotPlot.pdf"),width=10,height=15)
# ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = cname, y = Description)) + 
#   geom_point(aes(size = GeneRatio, color = p.adjust)) +
#   theme_bw(base_size = 11) +
#   #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#   scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.1))+
#   theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#   labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#   ylab(NULL)+ 
#   xlab("GeneRatio") 
# dev.off()



enrichRPath.res <- enrichPathway(gene=genes,universe=geneUniv)

res_dfRPath<-enrichRPath.res@result
res_dfRPath<-res_dfRPath%>% filter(qvalue<0.1)
res_dfRPath$GeneRatio<-sapply(res_dfRPath$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})



pdf(paste0(outFolder,"enrichRPath_SMC-1_DotPlot.pdf"),width=10,height=4)
ggplot(res_dfRPath, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+ 
  xlab("GeneRatio") 
dev.off()









################################################################################
# Wiki pathways- SMC-1
################################################################################
#library(rWikiPathways)
#wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")

library(magrittr)
library(clusterProfiler)


#gene <- names(geneList)[abs(geneList) > 2]
gene<-genes
#gene <- filter(res,padj<0.1,abs(log2FoldChange)>=0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique

wp2gene <- read.gmt("13_sample_investigation_plots/wikipathways-20191210-gmt-Homo_sapiens.gmt")
#wp2gene <-read.gmt(wpgmtfile)

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME



genes_WP289<-as.character(unique(wp2gene %>% filter(wpid=="WP289") %>% select(gene) %>% unlist) )



exist_genelist<-genes_WP289[which(genes_WP289 %in% names(geneList))]
res %>% filter(ENTREZID %in% exist_genelist &  padj<0.1)%>%select(gene_name)%>% unlist
#"IGFBP5"  "GUCY1A1"      "JUN" 



ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)



ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE,pvalueCutoff = 1 )#,minGSSize=5, maxGSSize=3000,eps =0,pvalueCutoff = 1)
head(ewp2)



core_genes<-ewp2@result %>% filter(ID=="WP289") %>% select(core_enrichment) %>% unlist
core_genes<-unlist(strsplit(core_genes,"/"))


length(res %>% filter ( ENTREZID %in% genes_WP289 & padj<0.1) %>% select(ENTREZID) %>% unlist)



################################################################################
# Wiki pathways- all cell types
################################################################################



res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
## cluster colors 

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
res$Cell_type<-clust2Names[res$Cell_type]

#ENTREZID
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

geneList <- -log10(res$pvalue)
names(geneList) <- res$ENTREZID
geneList = sort(geneList, decreasing = TRUE)





#qvalueCutoff  = 0.05

pathway_enrich<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)
  print(cname_select)
 
  aux <- res %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    print(length(gene))
    
    if (length(gene)>0)
    {
      
      ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
      
      if (!is.null(ewp))
      {
        res_en<-ewp@result
        res_en<-res_en%>%filter(p.adjust<0.1)
        dim1<-nrow(res_en)
        #print(dim1)
        if(min(dim1,10)>0)
        {
          
          
          res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
          res_en$Cell_type<-rep(cname_select,min(dim1,10))
        }
        return (res_en)
      }
      
    }
    
  }
 
}

res_enrichwikilist<-lapply(unique(res$Cell_type), function(x) {pathway_enrich(cname_select=x)})  
res_df_enrichwiki <- do.call(rbind,res_enrichwikilist)

res_df_enrichwiki$GeneRatio<-sapply(res_df_enrichwiki$GeneRatio, function(x){
  numden<-unlist(strsplit(x,"/"))
  return (as.numeric(numden[1])/as.numeric(numden[2]))
})





res_df<-res_df_enrichwiki[1:20,]

pdf(paste0(outFolder,"enrich_wikipathways_cname_DotPlot.pdf"),width=15,height=15)
ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = Cell_type, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  #theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
  ylab(NULL)+
  xlab(NULL)+
  coord_fixed(ratio = 1)+
  #theme_black()+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.y = element_text(hjust = 1))+
  theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 

#ggtitle("GO pathway enrichment")
dev.off()






pathway_GSEA<-function(res_gene=res,cname_select ,padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
  
  print(cname_select)

  
  aux <- res %>% filter(Cell_type==cname_select)
  
  if(nrow(aux)>0)
  {
    gene<-genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    #length(genes)
    
    print(length(gene))
    
    ewp <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE,pvalueCutoff = 1 )#,minGSSize=5, maxGSSize=3000,eps =0,pvalueCutoff = 1)
      
      if (!is.null(ewp))
      {
        res_en<-ewp@result
        res_en<-res_en%>%filter(p.adjust<0.1)
        dim1<-nrow(res_en)
        #print(dim1)
        if(min(dim1,10)>0)
        {
          
          
          res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"enrichmentScore","p.adjust")]
          res_en$Cell_type<-rep(cname_select,min(dim1,10))
        }
        return (res_en)
      }

  }
  
}

res_gseawikilist<-lapply(unique(res$Cell_type), function(x) {pathway_GSEA(cname_select=x)})  
res_df_gseawiki <- do.call(rbind,res_gseawikilist)



res_df_gseawiki<-res_df_gseawiki[1:15,]
pdf(paste0(outFolder,"gse_wikipathways_cname_DotPlot.pdf"),width=15,height=15)
ggplot(res_df_gseawiki, 
       aes(x = Cell_type, y = Description)) + 
  geom_point(aes(size = enrichmentScore, color = p.adjust)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
  labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
  ylab(NULL)+ 
  #theme_black()+
  theme_bw()+
  theme(text = element_text(size=30)) +
  theme(axis.text.x = element_text(angle = 45))+
  #xlab(NULL) +
  theme(axis.text.y = element_text(hjust = 1))+
  theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=10)) 

#+
#theme_black()
#coord_fixed(ratio = .8)
#ggtitle("GO pathway enrichment")
dev.off()









#######################################################################
# Forest plot
#######################################################################

outFolder="13_sample_investigation_plots/"

res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")


# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))

clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]
res<-res %>% filter(Cell_type=="Smooth muscle cells-1")
res<-res[order(abs(res$log2FoldChange), abs(res$baseMean),decreasing = TRUE),]
res<-res %>% filter (padj <=0.1)

res_select<-res %>% group_by(Cell_type) %>% filter (padj <=0.1 & abs(log2FoldChange)>=0.5)


DE_number_per_celltype<-tapply(res_select$padj ,factor(res_select$Cell_type),length )

res_select$DE_number_per_celltype<-DE_number_per_celltype[res_select$Cell_type]

res_select<-res_select[order(res_select$DE_number_per_celltype,decreasing = TRUE),]


res_select$Cell_type <- factor(res_select$Cell_type,levels=unique(res_select$Cell_type))




# forestPlot 
#coloring
res_select$padj[res_select$padj<1E-6]<- 1E-6

p<-res_select %>% 
  ggplot(aes(x=reorder(gene_name,(log2FoldChange),color="black"),y=log2FoldChange),color="black") +
  geom_errorbar(aes(ymax = log2FoldChange + 1.96*lfcSE, ymin = log2FoldChange - 1.96*lfcSE),size=1.2,color="black",alpha=1,width=0,position=position_dodge(width=0.2)) +
  geom_point(aes(size=padj),position=position_dodge(width=0.2),color="#111111") +
  ##scale_color_manual(guide = guide_legend(reverse = TRUE) ) +
  scale_size("q-values", trans="log10", range=c(7, 1),limits=c(1E-10,1), breaks=c(1E-12,1E-6,0.001,0.01,0.1)) +
  ## scale_alpha_manual(values=c(0.3, 1.0),guide=FALSE) +
  geom_hline(yintercept=0,lty=2) + 
  xlab("Gene") + ylab(expression(log[2](Fold~Change))) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0,hjust=0,vjust=0.5),strip.background=element_rect(fill="white",color="white")) +
  coord_flip() 
  #theme(strip.background =element_rect(fill=cluster.Colors))+
  #facet_grid(Cell_type ~ . ,scales = "free_y",space="free") 
ggsave(paste0(outFolder,"forestPlot_DEgenes_SMC-1.pdf"),p,width=10,height=30)



############### Forest plot for specific genes 

# so same code as in the forest plot, but you only show those genes, and for all cell-types available, 
# and you arrange/sort by gene name and color by cell-type, this way we can see something like the meta plot but using the forest plot.
# basically you can use the same code as in the forest plot, but you color by cell-type, and you can sort/arrange/facet by gene name.


cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")



outFolder="13_sample_investigation_plots/"

res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")


# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))

clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]

res_genes<-res[res$gene_name %in% c("PRKG1", "PDE5A",  "PLN"),]


p<-res_genes %>% 
  ggplot(aes(x=gene_name,y=log2FoldChange,color=Cell_type)) +
  geom_errorbar(aes(ymax = log2FoldChange + 1.96*lfcSE, ymin = log2FoldChange - 1.96*lfcSE),size=1.5,color="black",alpha=1,width=0,position=position_dodge(width=0.7)) +
  geom_point(aes(size=padj,color=Cell_type),position=position_dodge(width=0.2)) +
  ##scale_color_manual(guide = guide_legend(reverse = TRUE) ) +
  scale_size("q-values", trans="log10", range=c(7, 1),limits=c(1E-10,1), breaks=c(1E-12,1E-6,0.001,0.01,0.1)) +
  geom_hline(yintercept=0,lty=2) + 
  scale_color_manual(values=cluster.Colors) +
  xlab("Gene") + ylab(expression(log[2](Fold~Change))) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0,hjust=0,vjust=0.5),strip.background=element_rect(fill="white",color="white")) +
  coord_flip()+
  xlab(NULL)+
#facet_grid(gene_name ~ Cell_type) 
facet_grid(~gene_name ) 
#facet_grid(gene_name ~ . ,scales = "free_y",space="free") 
ggsave(paste0(outFolder,"forestPlot_metaplot.pdf"),p,width=10,height=5)



###############


#############################################################
# select genes from smooth muscle cell-1 and Mittal
############################################################
res<-res %>% filter(Cell_type=="Smooth muscle cells-1")


res<-res %>% filter(padj <0.1)
write.csv(res,file=paste0(outFolder,"DEpadj0.1_Smooth muscle cells-1.csv"))


res<-res %>% filter(padj <=0.05)
write.csv(res,file=paste0(outFolder,"DEpadj0.05_Smooth muscle cells-1.csv"))


ref_data<-load_ref_data(fl="myometrium_term_TL-TNL_ALLList")
ref_data<-ref_data %>% filter(Rpadj<=0.01 & abs(R.Log2FC)>=1)
write_csv(ref_data,file="13_network_analysis/Mittal_padj0.01_absl2fc1.0.csv")


ref_data<-load_ref_data(fl="myometrium_term_TL-TNL_ALLList")
ref_data<-ref_data %>% filter(Rpadj<=0.05 & abs(R.Log2FC)>=1)
write_csv(ref_data,file="13_network_analysis/Mittal_padj0.05_absl2fc1.0.csv")
