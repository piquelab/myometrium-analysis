############################################
# comparison with blood data based on negative and positive effects 
############################################
library(dplyr)
library(tidyverse)
library(mashr)
library(clusterProfiler)
#outFolder <- paste0("12_mash_analysis/non-shared/")
outFolder <- paste0("12_mash_analysis/")
system(paste0("mkdir -p ",outFolder))





#PosteriorMean_positive<-rownames(PosteriorMean)[which(PosteriorMean)>0]
#PosteriorMean_negative<-rownames(PosteriorMean)[which(PosteriorMean)<0]

outFolder2 <- paste0("12_mash_analysis/")
system(paste0("mkdir -p ",outFolder2))
write_rds(m.c, paste0(outFolder2,"m.c.rds"))
write_rds(PosteriorMean_positive, paste0(outFolder2,"PosteriorMean_positive.rds"))
write_rds(PosteriorMean_negative, paste0(outFolder2,"PosteriorMean_negative.rds"))



PosteriorMean<-m.c$result$PosteriorMean

PosteriorMean_positive<-apply(PosteriorMean, 1, function(x) {if (all(x>0)) return (TRUE) else return(FALSE) })
PosteriorMean_positive<-names(PosteriorMean_positive)[which(PosteriorMean_positive)]

PosteriorMean_negative<-apply(PosteriorMean, 1, function(x) {if (all(x<0)) return (TRUE) else return(FALSE) })
PosteriorMean_negative<-names(PosteriorMean_negative)[which(PosteriorMean_negative)]

consistent_effect_genes<-unique(c(PosteriorMean_negative,PosteriorMean_negative))
select_genes<-consistent_effect_genes


non_consistent_effect_genes<-rownames(m.c$result$PosteriorMean)[!rownames(m.c$result$PosteriorMean) %in% consistent_effect_genes]
select_genes<-non_consistent_effect_genes

load_ref_data<-function(fl)
{
    ref_data <- read.csv("TL-TNL_21vs28.csv",stringsAsFactors = FALSE)
    ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
    colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
    ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
    return(ref_data)
 
}




cor_with_ref<-function(experiment,ref_cutoff="All",plotcolor=TRUE,select_genes)
{
  
  
  res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
  res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
  
  
  eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  names(eg)[1]="gene_name"
  head(eg)
  e2g <- eg$gene_name
  names(e2g) <- eg$ENTREZID
  
  #ENTREZID
  res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))
  res <- res %>% filter(!is.na(padj) & !is.na(log2FoldChange)& !is.na(ENTREZID))
  res<-res %>% filter(kbid %in% select_genes)
  
  ref_data<-load_ref_data(fl="TL-TNL_21vs28") 
  
  correlation_result<-sapply(unique(res$Cell_type),function(x){
            
            result<-c(rep(NA,4))
            res <- res %>% filter(!is.na(padj) & !is.na(log2FoldChange)& !is.na(ENTREZID))
            res4 <- res %>% filter(padj<0.1)
           
            
            ref_data <- ref_data %>% filter(!is.na(Rpadj) & !is.na(R.Log2FC)& !is.na(ENTREZID))
            ref_data2 <- ref_data %>% filter(Rpadj<0.1)
            length(which(ref_data2$ENTREZID %in% unique(res4$ENTREZID)))
            
            
            res_rest<-res4 %>% filter(Cell_type!=x)
            
            res4 <- res4 %>% filter(Cell_type==x)
            
            if (ref_cutoff!="All")
              resJoin <- res4 %>% inner_join(ref_data2)  else resJoin <- res4 %>% inner_join(ref_data) 
            
            resJoin <- resJoin %>% filter(!is.na(padj))
            if(nrow(resJoin)>=50)  # 5 to 50
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
              
              
              
              if(ref_cutoff=="DE")
              {
                
                resall_celltype <- res4 %>% filter(Cell_type==x)
                resJoin <- resall_celltype %>% inner_join(ref_data) 
                cat(table(resJoin$ref_color),sep=":")
                resJoin$t<-resJoin$log2FoldChange/resJoin$lfcSE
                p2 <- resJoin %>% arrange(-padj) %>%
                  ggplot(aes(Rt,t)) +
                  geom_point(color="black")+ #aes(colour = ref_color)) +
                  ##       scale_color_manual(values=resJoin$ref_color)+
                  geom_smooth(method=lm, se=FALSE,linetype = "dashed", color="black")+
                  #scale_color_manual(name="Differentially expressed",values=c("None"="#CCCCCC","Only single cell"="#BDD7EE","Only bulk"="#DD99DD","Single cell and bulk"="#0000EE"))+
                  #scale_color_manual(name="Differentially expressed",values=c("#CCCCCC"="#CCCCCC","#BDD7EE"="#BDD7EE","#DD99DD"="#DD99DD","#0000EE"="#0000EE"))+
                  theme_bw()
                
                fname=paste0(outFolder,paste0(clust2Names[x],".png"))
                ggsave(fname,p2,width=6,height=4.5)
                
              }
              
              
              else
              {
                resall_celltype <- res %>% filter(Cell_type==x)
                resJoin <- resall_celltype %>% inner_join(ref_data) 
                # resJoin$ref_color <- "#CCCCCC" ## Light gray default no sig.
                # resJoin$ref_color[resJoin$padj<0.1] <- "#BDD7EE"  #light blue
                # resJoin$ref_color[resJoin$Rpadj<0.1] <- "#DD99DD" #purpule
                # resJoin$ref_color[resJoin$padj<0.1 & resJoin$Rpadj<0.1]="#0000EE" #blue
                
                print(length(!(resJoin$ENTREZID %in% res_rest$ENTREZID)))
                resJoin$ref_color <- "None" ## Light gray default no sig.
                resJoin$ref_color[resJoin$padj<0.1 ] <- "Only single cell"  #light blue
                resJoin$ref_color[resJoin$Rpadj<0.1 &  !(resJoin$ENTREZID %in% res_rest$ENTREZID)] <- "Only bulk" #purpule
                resJoin$ref_color[resJoin$padj<0.1 & resJoin$Rpadj<0.1]="Single cell and bulk" #blue
                
                
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
                  scale_color_manual(name="Differentially expressed",values=c("None"="#CCCCCC","Only single cell"="#BDD7EE","Only bulk"="#DD99DD","Single cell and bulk"="#0000EE"))+
                  #scale_color_manual(name="Differentially expressed",values=c("#CCCCCC"="#CCCCCC","#BDD7EE"="#BDD7EE","#DD99DD"="#DD99DD","#0000EE"="#0000EE"))+
                  theme_bw()
                
                fname=paste0(outFolder,paste0(clust2Names[x],".png"))
                ggsave(fname,p2,width=6,height=4.5)
                
              }

              result<-c(spearman_all,pearson_all)
            }
            return (result)
          })
          return (correlation_result)
          
}  


correlation_result<-cor_with_ref(experiment="TL-TNL_21vs28",ref_cutoff="All",plotcolor=TRUE,select_genes=select_genes)
### call the function to calculate the correlation:   

# correlation_result<-cor_with_ref(experiment="myometrium_term_TL-TNL_ALLList",ref_cutoff="All")
# correlation_result<-cor_with_ref(experiment="myometrium_term_TL-TNL_ALLList",ref_cutoff="DE")
# correlation_result<-cor_with_ref(experiment="TL-TNL_21vs28",ref_cutoff="All")



rownames(correlation_result)<-c("spearman_cor","spearman_pvalue","pearson_cor","pearson_pvalue")
colnames(correlation_result)<-as.character(clust2Names[unique(res$Cell_type)])
correlation_result<-correlation_result[,which(!is.na(correlation_result[1,]))]
correlation_result<-t(correlation_result)
cell_type<-rownames(correlation_result)
correlation_result_df<-as.data.frame(correlation_result)
correlation_result_df$cell_type<-cell_type

write.csv(correlation_result,file=paste0(outFolder,"correlation_result.csv"))


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

fname=paste0(outFolder,"barplot_cor.pdf")
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




