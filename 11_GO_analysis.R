###################################################
### pathway enrichment analysis

###################################################

library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)


outFolder <- paste0("11_pathway_enrichment_FC.0.5/whitebackground/")
system(paste0("mkdir -p ",outFolder))

######################################################################################################

pathway_enrich<-function(res_gene=res3,cname_select="CAM_T-cell_M",padj_cutoff=0.1,log2FoldChange_cutoff=0)
{
    cat ("=============== ",cname_select,"=============== ", "\n")
    pathway_enrich_cname_dir<-paste0(outFolder,cname_select,"/")
    system(paste0("mkdir -p ",pathway_enrich_cname_dir))
    result<-list()
    aux <- res_gene %>% filter(cname==cname_select)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    message(".................................")
    message("Number of DE genes: ",length(genes))
    #print(length(genes))
    
    message(".................................")
    message("enrichGO")
    ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
    print(head(ego))
    result$enrichGO<-ego
    #save(ego,file=paste0(pathway_enrich_cname_dir,"ego.RData"))
    #write.csv(ego,file=paste0(pathway_enrich_cname_dir,"ego.csv"))
    
    print(".................................")
    print("enrichKEGG")
    ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
    print(head(ekegg)) 
    result$enrichKEGG<-ekegg
    #save(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.RData"))
    #write.csv(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.csv"))
    
    
    message(".................................")
    message("enrichPathway")
    erpath <- enrichPathway(gene=genes,universe=geneUniv)
    print(head(erpath))
    result$enrichPathway<-erpath
    #save(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.RData"))
    # write.csv(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.csv"))
    
    message(".................................")
    message("gseGO")
    # BP: biological_process, CC: cellular_component, MF: molecular_function
    gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
    print(head(gseGO.res))
    result$gseGO<-gseGO.res
    #save(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.RData"))
    #write.csv(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.csv"))
    
    
    message(".................................")
    message("gsePathway")
    gseRPath.res <- gsePathway(geneList)
    print(head(gseRPath.res))
    result$gsePathway<-gseRPath.res
    return (result)
}



# load DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Cluster","Origin"),sep="_",remove=FALSE) #Cell_type


clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cluster]
res$cname<-paste0(res$Cell_type,"_",res$Origin)  
    

# Removing na pvalues
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 
res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    group_by(Cell_type,Origin) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))


#ENTREZID id 
eg = bitr(res2$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

# non-na ENTREZID were included
res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

#exploratory analysis
# # counting the number of DE genes per cname   
DE_per_cname<-sapply(unique(res3$cname), function(x,padj_cutoff=0.1,log2FoldChange_cutoff=0.5){
    aux <- res3 %>% filter(cname==x)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    length(genes)
    
})

cname_selected<-names(DE_per_cname)[which(DE_per_cname>0)]
result_pathway_en_list<-lapply(cname_selected, function(x) return(pathway_enrich(res3,x)))
names(result_pathway_en_list)<-cname_selected

save(result_pathway_en_list,file=paste0(outFolder,"pathwayEnrich_result.RData"))
which(DE_per_cname>0)


##########################################################################################  
#####                                  dot plot
##########################################################################################


load(paste0("11_pathway_enrichment_FC.0.5/pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)


DE_per_cname_select<-names(DE_per_cname)[which(DE_per_cname>=5)]
result_pathway_en_list<-result_pathway_en_list [DE_per_cname_select]
cname_selected<-names(result_pathway_en_list)
#gseGO
res_gseGO_list<-lapply(cname_selected, function(x)
{
    
    rs<-result_pathway_en_list[[x]]$gseGO@result %>% filter(qvalues<=0.05)
    dim1<-dim(rs)[1]

    if(min(dim1,5)>0)
    {
        res_en<-rs
        
        # to calculate GeneRatio=count/setSize
        
        #count
        gene_count<- res_en %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        
        ## merge with the original dataframe
        dot_df<- left_join(res_en, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
        
        dot_df<-dot_df[1:min(dim1,5),c("ID","Description" ,"enrichmentScore","p.adjust","GeneRatio")]
        dot_df$cname<-rep(x,min(dim1,5))
        dot_df
        }
    })
    

res_df_gseGO <- do.call(rbind,res_gseGO_list)

res_df_gseGO<-res_df_gseGO[1:15,]
pdf(paste0(outFolder,"gseGO_cname_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_gseGO, 
       aes(x = cname, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
    ylab(NULL)+ 
    xlab(NULL) +
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()


#enrichGO
res_enrichGO_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
    )  
    
res_df_enrichGO <- do.call(rbind,res_enrichGO_list)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO[1:15,]

pdf(paste0(outFolder,"enrichGO_cname_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = Description)) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()



####################

#enrichKEGG
res_enrichKEGG_list<-lapply(cname_selected, function(x)
{
    if(length(result_pathway_en_list[[x]]$enrichKEGG)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichKEGG@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
)  

res_df_enrichKEGG <- do.call(rbind,res_enrichKEGG_list)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})


res_df_enrichKEGG<-res_df_enrichKEGG[1:15,]
pdf(paste0(outFolder,"enrichKEGG_cname_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    #theme_black()+
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()


####################

#enrichPathway

res_enrichPathway_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichPathway)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichPathway@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
)  

res_df_enrichPathway <- do.call(rbind,res_enrichPathway_list)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichPathway<-res_df_enrichPathway[1:15,]
pdf(paste0(outFolder,"enrichPathway_cname_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
    aes(x = cname, y = Description)) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+
    xlab(NULL)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(axis.text=element_text(size=30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()

############################################################

### DE genes combined

############################################################


genes <- filter(res3,padj<0.1,abs(log2FoldChange)>0.5) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
geneUniv <- res3 %>% dplyr::select(ENTREZID) %>% unlist %>% unique

geneList <- -log10(res3$pvalue)
names(geneList) <- res3$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

message(".................................")
message("enrichGO")
ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
print(head(ego))

res_df_enrichGO<-ego@result %>% filter(qvalue<=0.05)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichGO<-res_df_enrichGO[1:15,]
pdf(paste0(outFolder,"enrichGO_combined_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()




print(".................................")
print("enrichKEGG")
ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
print(head(ekegg)) 

res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.05)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

res_df_enrichKEGG<-res_df_enrichKEGG[1:15,]
pdf(paste0(outFolder,"enrichKEGG.combined_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    theme_bw()+
    theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()

message(".................................")
message("enrichPathway")
erpath <- enrichPathway(gene=genes,universe=geneUniv)
print(head(erpath))

res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.05)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder,"enrichPathway.combined_DotPlot.pdf"),width=10,height=10)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
    aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    theme_bw()+
    theme(axis.text=element_text(size=30), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=30)) 
dev.off()













#