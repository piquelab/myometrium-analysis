library(tidyverse)
library(dplyr)
##library(knitr)
#library(DESeq2)
##library(annotables)
library(qqman)
library(org.Hs.eg.db)
library(clusterProfiler)



### Comparison with Mittal, Pooja, Roberto Romero, Adi L. Tarca, Juan Gonzalez, Sorin Draghici, Yi Xu, Zhong Dong et al. "Characterization of the myometrial transcriptome and biological pathways of spontaneous human labor at term." Journal of perinatal medicine 38, no. 6 (2010): 617-643. 






######################################################################### 
## load data
######################################################################### 
cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

names(cluster.Colors)<-paste0(c(0:23),"_",names(cluster.Colors))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
outFolder<-"./8_outputs_DESeq_Plots/"
system(paste0("mkdir -p ",outFolder))
# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
#res$Cell_type<-clust2Names[res$Cell_type]
#anno <- read_rds("3_MergeDemux_Output/anno.rds")
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

#ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))



load_ref_data<-function(outFolder="reference_Adi/",fl="CELLECTA.rds")
{
  
  if (fl=="TL-TNL_21vs28")
  {
    ref_data <- read.csv("TL-TNL_21vs28.csv",stringsAsFactors = FALSE)
    ref_data<-ref_data %>% select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
    ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
    return(ref_data)
    
  }
  if (fl=="myometrium_term_TL-TNL_ALLList")
  {
    ref_data <- read.delim("myometrium_term_TL-TNL_ALLList.txt")
    ref_data<-ref_data %>% select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
    ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
    ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
    return(ref_data)
  }
  else
  if(fl=="myometrium_bulk")
  {
    ref_data <- read_tsv("myo_bulk_TIN_TNL.txt")
    ref_data<-ref_data %>% select(SYMBOL,FoldChange,pval.fdr,ENTREZ ,t)
    colnames(ref_data) <- c("R.gene_name","R.Log2FC","Rpadj","ENTREZID","Rt")
    ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
    ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
    return(ref_data)
  }
  
  else if (fl=="PMID31921132")
  {
    ref_data<-read.delim("PMID31921132.txt")
    ref_data<-ref_data %>% select(SYMBOL,FC,P.Value,adj.P.Val,ENTREZ,t )
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
    ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
    ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
    return (ref_data)
    }
  else 
  {
    
    d<-readRDS(paste0(outFolder,fl))
    ref_data<-d[["LaborEffect"]]
    if(fl=="PCR.rds")
      ref_data<-ref_data %>% select(SYMBOL,logFC,P.Value,adj.P.Val,t)
    else
      ref_data<-ref_data %>% select(SYMBOL,log2FoldChange,P.Value,adj.P.Val,t)
    
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","Rt")
    return (ref_data)
  }
  
  }



#########################################################################  
## barplot and scater plot showing correlation with ref 
#########################################################################


## calculate the  correlation function
cor_with_ref<-function(experiment="RNASeq",ref_cutoff="All",plotcolor=TRUE)
{
  
  if (experiment=="TL-TNL_21vs28")
    ref_data<-load_ref_data(fl=experiment) 
  else 
  if (experiment=="myometrium_term_TL-TNL_ALLList")
    ref_data<-load_ref_data(fl=experiment) 
  else 
    if(experiment=="PMID31921132")
      ref_data<-load_ref_data(fl=experiment) 
    else
      if(experiment!="myometrium_bulk")
        ref_data<-load_ref_data(fl=paste0(experiment,".rds")) 
      else
    ref_data <- load_ref_data(fl="myometrium_bulk")
  
    if(experiment=="CELLECTA")
    outFolder<-paste0("./8_outputs_DESeq_Plots/DriverMap/",ref_cutoff,"/")  else outFolder<-paste0("./8_outputs_DESeq_Plots/", experiment,"/",ref_cutoff,"/")
    system(paste0("mkdir -p ",outFolder))

    
   # calculate the correlation 
  correlation_result<-sapply(unique(res$Cell_type),function(x){
    
    result<-c(rep(NA,4))
    
    #res_intersect<-res %>% filter(!gene_name %in% ref_data$gene_name)
    res <- res %>% filter(!is.na(padj) & !is.na(log2FoldChange)& !is.na(ENTREZID))
    res4 <- res %>% filter(padj<0.1)
    ref_data <- ref_data %>% filter(!is.na(Rpadj) & !is.na(R.Log2FC)& !is.na(ENTREZID))
    ref_data2 <- ref_data %>% filter(Rpadj<0.1)
    length(which(ref_data2$ENTREZID %in% unique(res4$ENTREZID)))
    
    res4 <- res4 %>% filter(Cell_type==x)
    
    
    if (ref_cutoff!="All")
    resJoin <- res4 %>% inner_join(ref_data2)  else resJoin <- res4 %>% inner_join(ref_data) 
      
    resJoin <- resJoin %>% filter(!is.na(padj))
    if(nrow(resJoin)>5)
    {
      print(nrow(resJoin))
      print(clust2Names[x])
      spearman_all<-cor.test(resJoin$log2FoldChange/resJoin$lfcSE,resJoin$Rt,method="spearman",na.rm=TRUE)
      spearman_all_pvalue<-as.numeric(spearman_all$p.value)
      spearman_all_cor<-as.numeric(spearman_all$estimate)
      spearman_all<-c(spearman_all_cor,spearman_all_pvalue)
      
      pearson_all<-cor.test(resJoin$log2FoldChange/resJoin$lfcSE,resJoin$Rt,method="pearson",na.rm=TRUE)
      pearson_all_pvalue<-as.numeric(pearson_all$p.value)
      pearson_all_cor<-as.numeric(pearson_all$estimate)
      pearson_all<-c(pearson_all_cor,pearson_all_pvalue)
      
      
      
      ## scater plot showing all genes 
      
      resall_celltype <- res %>% filter(Cell_type==x)
      resJoin <- resall_celltype %>% inner_join(ref_data)

      # color of points for scatter plot
      resJoin$ref_color <- "#CCCCCC" ## Light gray default no sig.
      resJoin$ref_color[resJoin$padj<0.1] <- "#BDD7EE"  #light blue 
      resJoin$ref_color[resJoin$Rpadj<0.1] <- "#800080" #purpule 
      resJoin$ref_color[resJoin$padj<0.1 & resJoin$Rpadj<0.1]="#0000EE" #blue  
      ## resJoin$ref_color<- sapply(resJoin$ENTREZID,function(y)
      ##   {
      ##    res_e<-resJoin %>% filter(ENTREZID ==y)
      ##    if (res_e$padj <0.1 & res_e$Rpadj <0.1) x<-"#0000EE" #blue
          
      ##      if (res_e$padj <0.1 & res_e$Rpadj >=0.1) x<-"#BDD7EE" #light blue
          
      ##      if (res_e$padj >=0.1 & res_e$Rpadj < 0.1) x<-"#800080" #purpule
         
      ##      if (res_e$padj >=0.1 & res_e$Rpadj >=0.1) x<- "gray"
      ##      return(x)
      ## })

      cat(table(resJoin$ref_color),sep=":")
      resJoin$t<-resJoin$log2FoldChange/resJoin$lfcSE
      p2 <- resJoin %>% arrange(-padj) %>%
        #ggplot(aes(R.Log2FC,log2FoldChange,color=padj<0.1)) +
        # scale_color_manual(values=c("gray","black")) +
        #ggplot(aes(R.Log2FC,log2FoldChange),color="black") +
        ggplot(aes(Rt,t,color=ref_color)) +
        geom_point()+ #aes(colour = ref_color)) +
 ##       scale_color_manual(values=resJoin$ref_color)+
        geom_smooth(method=lm, se=FALSE,linetype = "dashed", color="black")+
        theme_bw()
      
      fname=paste0(outFolder,paste0(clust2Names[x],"Roger.png"))
      cat(fname,"\n")
      ggsave(fname,p2,width=6,height=4.5)
      
      #cor_spearman_DE,cor_pearson_DE)
      result<-c(spearman_all,pearson_all)
      #print(length(res))
    }
    return (result)
  })
  return (correlation_result)

}  



### call the function to calculate the correlation:   

correlation_result<-cor_with_ref(experiment="myometrium_term_TL-TNL_ALLList",ref_cutoff="All")



rownames(correlation_result)<-c("spearman_cor","spearman_pvalue","pearson_cor","pearson_pvalue")
colnames(correlation_result)<-as.character(clust2Names[unique(res$Cell_type)])
correlation_result<-correlation_result[,which(!is.na(correlation_result[1,]))]
correlation_result<-t(correlation_result)
cell_type<-rownames(correlation_result)
correlation_result_df<-as.data.frame(correlation_result)
correlation_result_df$cell_type<-cell_type
  
  library(reshape2)
  library(ggplot2)
  dat2 <- melt(correlation_result[,c("spearman_cor","pearson_cor")])
  dat3 <- melt(correlation_result[,c("spearman_pvalue","pearson_pvalue" )])

  colnames(dat2)<-c("cluster","cor","value")
  colnames(dat3)<-c("cluster","pvalue","value")
  
  dat2$pvalue=dat3$value
  dat2$pvalue<-sapply(dat2$pvalue, function(x){
    if (x<=0.0001) x="****"
    else if (x<=0.001)x="***" 
    else if (x<=0.01) x="**"
    else if (x<=0.05) x="*"
    else if (x>0.05) return ("ns")
    return (x)
  })
  dat2<-dat2 %>%filter(cor=="spearman_cor")
  dat2$clustercolor<-as.character(cluster.Colors[dat2$cluster])
  
  dat2$clusternumber<-as.numeric(sapply(as.character(dat2$cluster),function(x){return(unlist(strsplit(x,"_"))[1])}))
  #dat2$cluster <- factor(dat2$cluster,levels=unique(dat2$cluster))
  dat2<-dat2[order(dat2$clusternumber,decreasing = FALSE),]
  
  fname=paste0(outFolder,"barplot_cor_v2.pdf")
  pdf(fname,width=12,height=6)
  ggplot(data=dat2, aes(x=cluster, y=value,fill=cluster)) +
    geom_bar(stat="identity",position="stack")+
    geom_text(aes(label=pvalue,vjust = -sign(value)), vjust=1.6, color="black", size=3.5)+
    theme_bw()+
    scale_fill_manual("legend", values = c("0_Stromal-1"="#DF7D99" ,"1_Macrophage-2"="#838EDF" ,"2_Macrophage-1"="#4E65A6","3_Endothelial-1"="#FFC000" ,"4_Monocyte"="#2BA3D3", "5_CD4_T-cell"="#9ABF5C" ,"6_Decidual"="#D14357" ,"7_CD8_T-cell"="#329B2D","8_LED"="#D5438E","9_Stromal-2"="#ED4315" ,"10_ILC"="#76956C" ,"11_NK-cell"="#7BC791","12_Smooth muscle cells-1"="#CA8588" ,"13_Myofibroblast"="#F88091" , "14_Macrophage-3"="#72C6C8" ,"15_Endothelial-2"="#E4652C" ,"16_DC"="#9B91B9" ,"17_Smooth muscle cells-2"="#A37584" ,"18_EVT"="#2C3E18" ,"19_Plasmablast"="#745B48" ,"20_Smooth muscle cells-3"="#AA5485" ,"21_Macrophage-4"="#4E747A","22_B-cell"="#C59A89","23_Unciliated Epithelial"="#C9C76F"))+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    theme(legend.position="none")+
    xlab("")+
    ylab("Spearman correlation")
  dev.off() 
  
  

  ## outFolder<-"reference_Adi/"
  ## DriverMap<-readRDS(paste0(outFolder,"CELLECTA.rds"))
  ## dim(DriverMap[["LaborEffect"]])
  ## 
  ## RNASeq<-readRDS(paste0(outFolder,"RNASeq.rds"))
  ## dim(RNASeq[["LaborEffect"]])
  ## 
  ## PCR<-readRDS(paste0(outFolder,"PCR.rds"))
  ## dim(PCR[["LaborEffect"]])
  
  
  #ref_data<-load_ref_data(fl="myometrium_bulk")
  # correlation_result<-cor_with_ref("myometrium_bulk")
  
  # ref_data<-load_ref_data(fl="RNASeq.rds")
  # outFolder<-"./8_outputs_DESeq_Plots/RNASeq/"
  # correlation_result<-cor_with_ref("RNASeq")
  # 
  
################################################################################################################################
## overlap between DE genes between ref and our study
###############################################################################################################################

ref_data<-load_ref_data(fl="myometrium_bulk")
ref_data<-load_ref_data(fl="RNASeq.rds")
 
ref_data<-load_ref_data(fl="myometrium_term_TL-TNL_ALLList")

barplot_data<-sapply(unique(res$Cell_type),function(x){
  res4 <- res %>% filter(Cell_type==x)
  res4 <- res4 %>% filter(padj<0.1)
  res4 <- res4 %>% filter(!is.na(padj))
  ref_data<-ref_data %>% filter(Rpadj<0.1 )
  ref_data <- ref_data %>% filter(!is.na(Rpadj))
  resJoin <- res4 %>% inner_join(ref_data) 
  resJoin <- resJoin %>% filter(!is.na(padj))
  
  
  bulk<-length(which(!ref_data$ENTREZID %in% res4$ENTREZID))
  single_cell<-length(which(!res4$ENTREZID %in% ref_data$ENTREZID))
  both<-nrow(resJoin)
  res<-c(both,single_cell,bulk)
  return(res)
  
})


rw<-c("both","single_cell")
barplot_data<-barplot_data[-3,]
clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Name<-paste0(c(0:23),"_",clust2Name)
names(clust2Name)<-c(0:23)
cl<-clust2Name[colnames(barplot_data)]
barplot_data<-t(barplot_data)
colnames(barplot_data)<-rw
rownames(barplot_data)<-cl
barplot_data<-barplot_data[order(barplot_data[,"single_cell"]),]

cl<-rownames(barplot_data)
rownames(barplot_data)<-NULL
par(mar=c(4, 14 ,4.1 ,2.1))
outFolder<-"./8_outputs_DESeq_Plots/"

fname=paste0(outFolder,"comparison_with_bulk.pdf");
pdf(fname,width=3,height=4)
rownames(barplot_data)<-NULL

out<-barplot(t(barplot_data),beside=FALSE,horiz=TRUE,col = c("#0000EE","#BDD7EE"),axis.lty=1,las=1,xlab=seq(0,max(barplot_data),by=500),xlim=c(0,2500))
legend("bottomright", legend=c("Single cell and bulk analyses combined", "Single cell analysis"), fill=c("#0000EE","#BDD7EE"),bty = "n",title="",cex=1)
text(out, cl, pos=2, xpd=TRUE, cex=.8)
#out
#abline(v=8,col = "grey",lty = 5)
#axis(side=1, at=out, labels=seq(0,max(barplot_data),by=500))
dev.off()



library(reshape2)
library(ggplot2)
dat2 <- melt(barplot_data)
colnames(dat2)<-c("cluster","exp","value")
ggplot(dat2,aes(y=cluster,fill=exp,x=value)) +
#geom_bar(position="stack") +
geom_bar(position="stack", stat="identity")
#guides(colour = guide_legend(override.aes = list(size=5),title="# DE genes")) +
#scale_color_manual("# DE genes",values=c("both"="#333399","single_cell"="#A50021","bulk"="seagreen"))


############################################################
# intersected genes
# to do: meta plot and boxplot
############################################################
ref_data<-load_ref_data(fl="RNASeq.rds")
outFolder<-"./8_outputs_DESeq_Plots/RNASeq/"
correlation_result<-cor_with_ref("RNASeq")


ref_data<-load_ref_data(fl="CELLECTA.rds")
outFolder<-"./8_outputs_DESeq_Plots/DriverMap/"
# correlation_result<-cor_with_ref("CELLECTA")

intersected_genes<-c()
for (x in unique(res$Cell_type))
{
  res4 <- res %>% filter(Cell_type==x)
  res4 <- res4 %>% filter(padj<0.1)
  res4 <- res4 %>% filter(!is.na(padj))
  ref_data<-ref_data %>% filter(Rpadj<0.1 )
  ref_data <- ref_data %>% filter(!is.na(Rpadj))
  resJoin <- res4 %>% inner_join(ref_data) 
  resJoin <- resJoin %>% filter(!is.na(padj))
  intersected_genes<-rbind(intersected_genes,resJoin)
}
write.csv(intersected_genes,file=paste0(outFolder,"intersected_genes.csv"))

intersected_genes<-read.csv(paste0(outFolder,"intersected_genes.csv"),stringsAsFactors = FALSE)
