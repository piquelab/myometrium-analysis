
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ReactomePA)
######################################################################
# Using the ordered pvalue form the reference, and color each line by gene-category,
# The gene-category would be DEG for each cell-type in our data set.
######################################################################

#you may be able to this plot without the package, even, I think you just need to plot the rank of the p-value,  and caluclate the running enichment score. 
# but you may as well use the package, the only thing is that you will not use any database, you will have your list of genes from the single cell analysis.

outFolder <- paste0("./9_enrichment_Plots/")
system(paste0("mkdir -p ",outFolder))


# load data 


########################################################
# load single cell data 
########################################################
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")

res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)



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

####################################################################
# load reference data 
####################################################################
load_ref_data<-function(fl="CELLECTA.rds")
{
    if(fl=="TLvsTNL_blood_ENTREZ")
    {
        ref_data <- read.csv("TLvsTNL_blood_ENTREZ.csv",stringsAsFactors = FALSE)
        ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
        ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
        return(ref_data)
        
    }
    else 
    if (fl=="TL-TNL_21vs28")
    {
        ref_data <- read.csv("TL-TNL_21vs28.csv",stringsAsFactors = FALSE)
        ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
        ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
        return(ref_data)
        
    }
    else
    if (fl=="myometrium_term_TL-TNL_ALLList")
    {
        ref_data <- read.delim("myometrium_term_TL-TNL_ALLList.txt")
        ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t )
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
        ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
        ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
        return(ref_data)
    }
    else
        if(fl=="myometrium_bulk")
        {
            ref_data <- read_tsv("myo_bulk_TIN_TNL.txt")
            ref_data<-ref_data %>% dplyr::select(SYMBOL,FoldChange,pval.fdr,ENTREZ ,t)
            colnames(ref_data) <- c("R.gene_name","R.Log2FC","Rpadj","ENTREZID","Rt")
            ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
            ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
            return(ref_data)
        }
    
    else if (fl=="PMID31921132")
    {
        ref_data<-read.delim("PMID31921132.txt")
        ref_data<-ref_data %>% dplyr::select(SYMBOL,FC,P.Value,adj.P.Val,ENTREZ,t )
        colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
        ref_data <- ref_data %>% filter(!is.na(R.Log2FC) & !is.na(ENTREZID)  & !is.na(Rpadj))
        ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
        return (ref_data)
    }
    else 
    {
        
        d<-readRDS(paste0("reference_Adi/",fl,".rds"))
        ref_data<-d[["LaborEffect"]]
        if(fl=="PCR")
            
        {
            ref_data<-ref_data %>% dplyr::select(SYMBOL,logFC,P.Value,adj.P.Val,t)
            colnames(ref_data)<-c("gene_name","R.Log2FC","Rpvalue","Rpadj","Rt")
            #ref_data <- ref_data %>% left_join(eg) %>% filter(!is.na(ENTREZID))
            ref_data <- ref_data %>% left_join(eg) %>% filter(!is.na(ENTREZID))
            
            
        }
        
        else
        {
            ref_data<-ref_data %>% dplyr::select(SYMBOL,log2FoldChange,P.Value,adj.P.Val)
            colnames(ref_data)<-c("gene_name","R.Log2FC","Rpvalue","Rpadj")
            ref_data <- ref_data %>% left_join(eg) %>% filter(!is.na(ENTREZID))
        }
        return (ref_data)
    }
    
}


#############################################################################
# gene set enrichment analysis
#############################################################################
# reference data 
experiment<-"myometrium_term_TL-TNL_ALLList"
#experiment<-"TL-TNL_21vs28"
#experiment<-"CELLECTA"
#experiment<-"RNASeq"

experiment<-"TLvsTNL_blood_ENTREZ"
ref_data<-load_ref_data(fl=experiment) 



# colnames(ref_data)[which(colnames(ref_data)=="ENTREZID")]<-"R.ENTREZID"
# colnames(ref_data)[which(colnames(ref_data)=="R.gene_name")]<-"gene_name"



# cell_markers <- read_tsv('Human_cell_markers.txt') 
# cell_markers<-cell_markers%>%
#     tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
#     dplyr::select(cellMarker, geneID) %>%
#     dplyr::mutate(geneID = strsplit(geneID, ', '))
# cell_markers



resDE<-res %>% filter(padj<0.1 ) #single cell fdr 0.1
resDE<-resDE%>% dplyr::select(Cell_type, ENTREZID)
colnames(resDE)<-c("cellMarker","geneID")

#cell_markers<-resDE %>% dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers <- resDE %>% group_by(cellMarker) %>% summarize(geneID=list(unique(geneID)))

#genelist_origin <- -log10(ref_data$Rpadj)
genelist_origin <- -log10(ref_data$Rpvalue)
names(genelist_origin) <- ref_data$ENTREZID

genelist_df<-as.data.frame(cbind(names(genelist_origin),as.numeric(genelist_origin)))
colnames(genelist_df)<-c("id","value")
genelist_df<-genelist_df %>%   group_by(id) %>%   sample_n(1)

genelist<-as.numeric(genelist_df$value)
names(genelist)<-genelist_df$id
#genelist<-genelist_origin[which(table(names(genelist_origin))==1)]



genelist = sort(genelist, decreasing = TRUE)

gene <- names(genelist)[genelist > 2]



system(paste0("mkdir -p ",outFolder,experiment))

GSEA.res <- GSEA(genelist, TERM2GENE=cell_markers,minGSSize=5, maxGSSize=3000,eps =0,pvalueCutoff = 1)
res_df<-GSEA.res@result
#res_df<-res_df %>% filter(p.adjust<0.1)
fname<-paste0(outFolder,experiment,"/GSEA-cellmarker.png")
p1<-ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
           aes(x = enrichmentScore, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(min(res_df$p.adjust)-0.001,max(res_df$p.adjust)+0.001), low="red",high = "blue")+
    #scale_colour_gradient2(limits=c(0,1), low="red",high = "blue",midpoint=median(res_df$p.adjust),mid="gray")+
    #scale_color_gradient2(low = "red",  high = "blue", space = "Lab" ,mid = "gray",midpoint=0.1)+
    
    scale_color_gradient(limits=c(0,1),low = "red",  high = "blue", space = "Lab" )+
    #scale_color_gradient2(limits=c(min(res_df$p.adjust)-0.001,max(res_df$p.adjust)+0.001),low = "red",  high = "blue", space = "Lab" ,mid = "gray",midpoint=0.1)+
    
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
    ylab(NULL) 
ggsave(fname,p1,width=6,height=5)
save(GSEA.res,file=paste0(fname,"GSEA.res.RData"))


cnetplot(GSEA.res, foldChange=genelist)




# #wrong
# p11<-ggplot(res_df, aes(x=enrichmentScore, y=Description, fill=p.adjust)) +
#     #geom_tile() +
#     geom_point(aes(size = enrichmentScore)) +
#     theme_bw(base_size = 11) +
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     #scale_fill_gradientn(colours = c("red", "gray", "blue"), values = c(0, 0.025,0.050,0.075,0.1, 0.2,0.4,0.7,1))+
#     scale_fill_gradient(low="red",high="blue")+
#     labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)

fname<-paste0(outFolder,experiment,"/gseaplot.png")
p2 <- gseaplot2(GSEA.res, geneSetID = 1:12, subplots = 1,color=cluster.Colors[GSEA.res@result$ID])
ggsave(fname,p2,width=12,height=10)
#p3 <- gseaplot2(GSEA.res, geneSetID = 1:6, subplots = 1:2)
#cowplot::plot_grid(p2, p3, ncol=1)

write.csv(res_df,file=paste0(outFolder,experiment,"/gsea_df.csv"))


# enriche.res<- enricher(gene, TERM2GENE=cell_markers, minGSSize=1)
# res_df<-enriche.res@result
# #res_df<-res_df %>% filter(p.adjust<=0.05)
# res_df$GeneRatio<-sapply(res_df$GeneRatio, function(x){
#     numden<-unlist(strsplit(x,"/"))
#     return (as.numeric(numden[1])/as.numeric(numden[2]))
# })
# 
# fname<-paste0(outFolder,experiment,"/enrich-cellmarker.png")
# p3<-ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
#            aes(x = GeneRatio, y = Description)) +
#     geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
#     labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)+
#     xlab(NULL)+
#     theme_bw()+
#     theme(axis.text.y = element_text(hjust = 1))+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(fname,p3,width=6,height=5)
# 


#########################################################
# binomial test
#########################################################
resDE<-res %>% filter(padj<0.1 ) #single cell fdr 0.1
total<-table(resDE$Cell_type)

res_up<-resDE%>%filter(log2FoldChange>0)
upregulated<-rep(0,length(total))
names(upregulated)<-names(total)
upregulated[names(table(res_up$Cell_type))]<-table(res_up$Cell_type)


sc <- read_rds("4_harmony/sc.NormByLibrary.Harmony.StringentFiltering.res0.8.rds")
md <- read_rds("./4_harmony_cellClass_PBMC/sc.NormByLocation.ref.Anchors.rds") %>%
    as.data.frame %>%
    rownames_to_column("BARCODES") %>%
    dplyr::select(BARCODES,scLabor_ID=predicted.celltype.l2,scLabor_Score=predicted.celltype.l2.score)
md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
    left_join(md) 
identical(md$BARCODES,rownames(sc@meta.data))
sc@meta.data$cluster_name <- clust2Name[sc@meta.data$seurat_clusters]
cell_counts<-table(sc$cluster_name)



binom.test.res<-c()

for( x in unique(resDE$Cell_type))
{
   
    btest<-binom.test(upregulated[x],total[x],0.5)
    if(is.na(upregulated[x]))
        upregulated[x]<-0
    rb<-as.numeric(c(cell_counts[x],upregulated[x], (total[x]-upregulated[x]),total[x],btest$p.value))
    binom.test.res<-rbind(binom.test.res,rb)
    print(x)
    print(btest$p.value)
    }

rownames(binom.test.res)<-unique(resDE$Cell_type)
colnames(binom.test.res)<-c("Cell-counts","Up-regulated","Down-regulated","Total","P-value")
binom.test.res<-as.data.frame(binom.test.res)
binom.test.res$padj<-p.adjust(binom.test.res$`P-value`,"fdr")
binom.test.res<-binom.test.res[order(binom.test.res[,"padj"],decreasing = FALSE),]

write.csv(binom.test.res,file="8_outputs_DESeq_Plots/binom.test.res.csv")
# res_down<-resDE%>%filter(Cell_type =="3_Endothelial-1" & log2FoldChange<=0)
# binom.test(753,1653,0.5)

###############################################################################################
# # reference data 
# 
# experiment<-"myometrium_term_TL-TNL_ALLList"
# #experiment<-"TL-TNL_21vs28"
# 
# ref_data<-load_ref_data(fl=experiment) 
# 
# colnames(ref_data)[which(colnames(ref_data)=="ENTREZID")]<-"R.ENTREZID"
# colnames(ref_data)[which(colnames(ref_data)=="R.gene_name")]<-"gene_name"
# resDE<-res %>% filter(padj<0.1 ) #single cell fdr 0.1
# res_join<-ref_data %>% inner_join(resDE)
# 
# ref_data_DE<-ref_data %>% filter(Rpadj<0.1 )
# res_join_DE<-ref_data_DE %>% inner_join(resDE)
# 
# 
# #system(paste0("mkdir -p ",outFolder,experiment))
# system(paste0("mkdir -p ",outFolder,experiment,"/enrichGO"))
# system(paste0("mkdir -p ",outFolder,experiment,"/enrichKEGG"))
# system(paste0("mkdir -p ",outFolder,experiment,"/enrichPathway"))
# gseaplot_list<-list()
# for (i in 1: length(unique(res_join$Cell_type)))
# {
#     x<-unique(res_join$Cell_type)[i]
#     res_join_celltype<-res_join %>%filter(Cell_type %in% x)
#     #geneList <- -log10(res_join_celltype$Rpadj)
#     # geneList <- res_join_celltype$R.Log2FC
#     #names(geneList) <- res_join_celltype$ENTREZID
#     
#     ## geneList_norepeat<-tapply(geneList, as.factor(names(geneList)),function(x)
#     ##     {
#     ##     return (max(x))
#     ## })
#     ## 
#     ## geneList<-geneList_norepeat
#     ## geneList_norepeat<-as.numeric(geneList_norepeat)
#     ## names(geneList_norepeat)<-names(geneList)
#     ## geneList<-geneList_norepeat
#     
#     # geneList = sort(geneList, decreasing = TRUE)
#     # edo2 <- gseNCG(geneList, nPerm=10000)
#     # print(nrow(edo2@result))
#     # if(nrow(edo2@result)>0)
#     #     p1 <- gseaplot2(edo2, geneSetID = 1, subplots = 1)  #gseaplot_list[[i]]<-
# 
#    
#     
#     
#     
#     geneUniv <- res_join_celltype %>% dplyr::select(ENTREZID) %>% unlist
#     #res_join_DE_celltype<-res_join_DE %>%filter(Cell_type %in% x)
#     
#     geneList_DE <- ref_data$Rpadj
#     names(geneList_DE) <- ref_data$ENTREZID
#     genes<-names(geneList_DE)[geneList_DE<=0.1 ]
#         
#     message(".................................")
#     message("enrichGO")
#     ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
#     print(head(ego))
#     
#     if (length(ego)>0)
#     {
#         res_df_enrichGO<-ego@result %>% filter(qvalue<=0.05)
#         
#         if (nrow(res_df_enrichGO)>0 )
#         {
#             mx<-min(nrow(res_df_enrichGO),15)
#             res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
#                 numden<-unlist(strsplit(x,"/"))
#                 return (as.numeric(numden[1])/as.numeric(numden[2]))
#             })
#             
#             res_df_enrichGO<-res_df_enrichGO[1:mx,]
#             fname<-paste0(outFolder,experiment,"/enrichGO","/enrichGO_dotplot_",x,".png")
#             #pdf(paste0(outFolder,experiment,"/enrichGO","/enrichGO_dotplot_",x,".pdf"),width=9,height=5)
#             p1<-ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
#                    aes(x = GeneRatio, y = Description)) + 
#                 geom_point(aes(size = GeneRatio, color = p.adjust)) +
#                 theme_bw(base_size = 14) +
#                 scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#                 theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#                 labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#                 ylab(NULL)+ 
#                 xlab(NULL)+
#                 theme_bw()+
#                 theme(axis.text.y = element_text(hjust = 1))+
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#              ggsave(fname,p1,width=9,height=5)
#             
#             
#         }
#     }
#     
#         
#     print(".................................")
#     print("enrichKEGG")
#     ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
#     print(head(ekegg))
# 
#     if(length(ekegg)>0)
#     {
#         res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.05)
# 
#         if (nrow(res_df_enrichKEGG)>0 )
#         {
#             mx<-min(15,nrow(res_df_enrichKEGG))
#             res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
#                 numden<-unlist(strsplit(x,"/"))
#                 return (as.numeric(numden[1])/as.numeric(numden[2]))
#             })
# 
#             res_df_enrichKEGG<-res_df_enrichKEGG[1:mx,]
# 
#             #pdf(paste0(outFolder,experiment,"/enrichKEGG","/enrichKEGG_dotPlot_",x,".pdf"),width=9,height=5)
#             fname<-paste0(outFolder,experiment,"/enrichKEGG","/enrichKEGG_dotPlot_",x,".png")
#             p1<-ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
#                    aes(x = GeneRatio, y = Description)) +
#                 geom_point(aes(size = GeneRatio, color = p.adjust)) +
#                 theme_bw(base_size = 14) +
#                 scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#                 labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#                 ylab(NULL)+
#                 xlab(NULL)+
#                 theme_bw()+
#                 theme(axis.text.y = element_text(hjust = 1))+
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#             ggsave(fname,p1,width=9,height=5)
#             #dev.off()
#         }
#     }
#     # 
#     # 
#     # 
#     # 
#     # 
#     # 
#     message(".................................")
#     message("enrichPathway")
#     erpath <- enrichPathway(gene=genes,universe=geneUniv)
#     print(head(erpath))
# 
#     if(length(erpath)>0)
#     {
#         res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.05)
# 
#         if (nrow(res_df_enrichPathway)>0 )
#         {
# 
#             mx<-min(15,nrow(res_df_enrichPathway))
# 
#             res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
#                 numden<-unlist(strsplit(x,"/"))
#                 return (as.numeric(numden[1])/as.numeric(numden[2]))
#             })
# 
#             res_df_enrichPathway<-res_df_enrichPathway[1:mx,]
#             fname<-paste0(outFolder,experiment,"/enrichPathway","/enrichPathway_",x,".png")
#             #pdf(paste0(outFolder,experiment,"/enrichPathway","/enrichPathway_",x,".pdf"),width=9,height=5)
#             p1<-ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
#                    aes(x = GeneRatio, y = Description)) +
#                 geom_point(aes(size = GeneRatio, color = p.adjust)) +
#                 theme_bw(base_size = 11) +
#                 scale_color_gradient(low = "red",  high = "blue", space = "Lab")+
#                 theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
#                 labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
#                 ylab(NULL)+
#                 xlab(NULL)+
#                 theme_bw()+
#                 theme(axis.text.y = element_text(hjust = 1))+
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#             ggsave(fname,p1,width=9,height=5)
#             #dev.off()
#         }
# 
#     }
#     # 
#     
#     
# }
# 
# # edo <- enrichDGN(genes)
# # print(nrow(edo@result))
# # barplot(edo, showCategory=20)
# 
# 
# 
# #################  all genes in reference 
# 
# 
# 
# 
# 
# ########################################################################################################
# # ORA and GSEA dot plots
# 
# experiment<-"myometrium_term_TL-TNL_ALLList"
# # experiment<-"TL-TNL_21vs28"
# 
# ref_data<-load_ref_data(fl=experiment) 
# 
# colnames(ref_data)[which(colnames(ref_data)=="ENTREZID")]<-"R.ENTREZID"
# colnames(ref_data)[which(colnames(ref_data)=="R.gene_name")]<-"gene_name"
# 
# 
# geneList <- ref_data$R.Log2FC
# names(geneList) <- ref_data$R.ENTREZID
# geneList = sort(abs(geneList), decreasing = TRUE)
# 
# 
# message(".................................")
# message("gsePathway")
# gseRPath.res <- gsePathway(geneList)
# print(head(gseRPath.res))
# 
# #plot
# res_df<-gseRPath.res@result[1:20,]
# res_df<-res_df %>% filter(qvalues<=0.05)
# system(paste0("mkdir -p ",outFolder,experiment))
# pdf(paste0(outFolder,experiment,"/",experiment,".gseRPath.dotPlot.pdf"),width=30,height=20)
# ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
#        aes(x = enrichmentScore, y = Description)) +
#     geom_point(aes(size = enrichmentScore, color = p.adjust)) +
#     theme_bw(base_size = 11) +
#     #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#     scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
#     theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
#     labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
#     ylab(NULL)
# dev.off()
# 
# 
# 
# #edo2 <- gseNCG(geneList, nPerm=10000)
# 


####################### permutation based on phenotype

# package RGSEA

#https://www.bioconductor.org/packages/release/bioc/manuals/RGSEA/man/RGSEA.pdf

RGSEAfix(query, reference, queryclasses, refclasses, random = 5000, featurenum
         = 500, iteration = 100)


# package Category 
# https://www.rdocumentation.org/packages/Category/versions/2.38.0/topics/gseattperm

haveALL <- require("ALL")
gseattperm(eset, fac, mat, nperm)



#https://rdrr.io/github/GSEA-MSigDB/GSEA_R/man/GSEA.html

#https://rdrr.io/github/GSEA-MSigDB/GSEA_R/
    
install.packages("remotes")
remotes::install_github("GSEA-MSigDB/GSEA_R")
# this option is prefered
GSEA(input.ds, input.cls, input.chip = "NOCHIP", gene.ann = "", gs.db,
     gs.ann = "", output.directory = getwd(),
     doc.string = "gsea_result", reshuffling.type = "sample.labels",
     nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1,
     fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, topgs = 20,
     adjust.FDR.q.val = F, gs.size.threshold.min = 15,
     gs.size.threshold.max = 500, reverse.sign = F, preproc.type = 0,
     random.seed = as.integer(Sys.time()), perm.type = 0, fraction = 1,
     replace = F, collapse.dataset = FALSE,
     collapse.mode = "NOCOLLAPSE", save.intermediate.results = F,
     use.fast.enrichment.routine = T, gsea.type = "GSEA",
     rank.metric = "S2N")


# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading

geneset<-read_csv("9_enrichment_Plots/geneset.csv",col_names = FALSE)
geneset<-geneset[-1,]
geneset<-t(geneset)
write.csv(geneset,"9_enrichment_Plots/geneset_gmtformat.csv")
