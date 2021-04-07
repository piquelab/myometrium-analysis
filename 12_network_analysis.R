library(tidyverse)
##library(knitr)
library(STRINGdb)
library(igraph)
library(ggraph)
library(igraph)

outFolder <- paste0("./13_network_analysis/")
system(paste0("mkdir -p ",outFolder))

## selected genes

################################################
# load data
################################################
#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")

res<-res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-c(0:23)
res$Cell_type<-clust2Names[res$Cell_type]



################################################
# string db
################################################
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=200, input_directory="")
aliases<-string_db$get_aliases()
#aliases<-get_aliases(res$gene_name, db = org.Hs.eg.db::org.Hs.eg.db)
aliases_map<-aliases$STRING_id
names(aliases_map)<-aliases$alias
res$STRING_id<-aliases_map[res$gene_name]


cell_types<-unique(res$Cell_type)



selected_ensm_join_list<-list()
full.graph <- string_db$get_graph()
V(full.graph)$name<-names(aliases_map)[which(aliases_map %in% V(full.graph)$name)]


res_allgenes<-res$STRING_id[which(res$gene_name %in% V(full.graph)$name)]
res_allgenes<-unique(res_allgenes)
res_allgenes<-res_allgenes[!is.na(res_allgenes)]
subnetwork_res_allgenes<-string_db$get_subnetwork(res_allgenes)
#V(subnetwork_res_allgenes)$name<-names(aliases_map)[which(aliases_map %in% V(subnetwork_res_allgenes)$name)]
hs_allgenes <- hub_score(subnetwork_res_allgenes, weights=NA)$vector
hs_allgenes<-hs_allgenes[order(hs_allgenes,decreasing = TRUE)]
#names(hs_allgenes)[1:100]

# for global graph 
res_selected_ensm<-res %>% filter(padj <= 0.01 & log2FoldChange>2)

#adding hub genes to selected genes
selected_ensm<-unique(c(res_selected_ensm$STRING_id,names(hs_allgenes)[1:20]))
#neighbors<-string_db$get_neighbors(selected_ensm)
subnetwork<-string_db$get_subnetwork(selected_ensm)
V(subnetwork)$name<-names(aliases_map)[which(aliases_map %in% V(subnetwork)$name)]
nm<-V(subnetwork)$name
hs <- hub_score(subnetwork, weights=NA)$vector
hub_scores<-as.data.frame(hs)
hub_scores$gene_name<-rownames(hub_scores)
hub_scores<-hub_scores[order(hub_scores[,"hs"],decreasing = TRUE),]
as <- authority_score(subnetwork, weights=NA)$vector
authority_scores<-as.data.frame(as)
authority_scores$gene_name<-rownames(authority_scores)
authority_scores<-authority_scores[order(authority_scores[,"as"],decreasing = TRUE),]
node_size<-hub_scores$hs[hub_scores$gene_name %in% as.character(V(subnetwork)$name)]
names(node_size)<-hub_scores$gene_name[hub_scores$gene_name %in% as.character(V(subnetwork)$name)]
node_size<-node_size[as.character(V(subnetwork)$name)]
#### draw the global graph## 

nm[which(node_size<0.5)]<-""
myfill<-rep("yellow",length(node_size))
myfill[which(node_size>0.5)]<-"orange"

#logfc for mapping 
# color based on the celltypes
#more genes 


graph_tbl <- subnetwork %>% 
  as_tbl_graph() %>% 
  activate(nodes) %>% 
  mutate(degree  = centrality_degree()) 


########################################
# plot graph
########################################
ggraph(subnetwork,layout = "kk") +  #fr
  geom_edge_link(colour = "gray",show.legend = FALSE,alpha=0.2) +
  geom_node_point(aes(size = node_size),color=myfill,shape=16,show.legend = FALSE) + 
  #scale_size(range = c(20, 30)) +
  #scale_color()+
 # theme_graph()+
  #scale_size(range = c(1,6))+
  #geom_node_point(aes(size = centrality, colour = centrality)) + 
geom_node_text(aes(label = nm), colour = 'black', vjust = -0.8) 


                        

########################################
# systematic way of selecting genes to plot 
########################################

gene_pvalue_membership<-matrix(0,nrow=length(unique(res$gene_name)),ncol=length(cell_types))
rownames(gene_pvalue_membership)<-unique(res$gene_name)
colnames(gene_pvalue_membership)<-cell_types

gene_logfc_membership<-matrix(0,nrow=length(unique(res$gene_name)),ncol=length(cell_types))
rownames(gene_logfc_membership)<-unique(res$gene_name)
colnames(gene_logfc_membership)<-cell_types


for ( i in 1:length(cell_types))
{
  print(i)
  res_selected_ensm<-res %>% filter(Cell_type ==cell_types[i] & padj <= 0.05 & log2FoldChange>2)# 
  selected_ensm<-res_selected_ensm$STRING_id
  
  subnetwork<-string_db$get_subnetwork(selected_ensm)
  V(subnetwork)$name<-names(aliases_map)[which(aliases_map %in% V(subnetwork)$name)]
  hs <- hub_score(subnetwork, weights=NA)$vector
  hub_scores<-as.data.frame(hs)
  hub_scores$gene_name<-rownames(hub_scores)
  
  as <- authority_score(subnetwork, weights=NA)$vector
  authority_scores<-as.data.frame(as)
  authority_scores$gene_name<-rownames(authority_scores)


}




 
# for ( i in 1:length(cell_types))
# {
#   print(i)
#   res_selected_ensm<-res %>% filter(Cell_type ==cell_types[i])# & padj <= 0.05 & log2FoldChange>1
#   selected_ensm<-res_selected_ensm$STRING_id
#   #neighbors<-string_db$get_neighbors(selected_ensm)
#   subnetwork<-string_db$get_subnetwork(selected_ensm)
#   V(subnetwork)$name<-names(aliases_map)[which(aliases_map %in% V(subnetwork)$name)]
#   #string_db$plot_network( selected_ensm )
#   hs <- hub_score(subnetwork, weights=NA)$vector
#   hub_scores<-as.data.frame(hs)
#   hub_scores$gene_name<-rownames(hub_scores)
#   
#   as <- authority_score(subnetwork, weights=NA)$vector
#   authority_scores<-as.data.frame(as)
#   authority_scores$gene_name<-rownames(authority_scores)
#   
#   deg<-degree(subnetwork)
#   
#   degree_scores<-as.data.frame(deg)
#   degree_scores$gene_name<-rownames(degree_scores)
#   
#   selected_ensm_join<-res_selected_ensm %>% inner_join(degree_scores)%>% inner_join(hub_scores)%>% inner_join(authority_score)
#   
#   selected_ensm_join <- selected_ensm_join %>% filter(!is.na(pvalue)) %>%
#   arrange(desc(hs),desc(deg),desc(abs(log2FoldChange)),pvalue)
#   
#   
#   
#   selected_ensm_join$padj_map<-sapply(selected_ensm_join$padj, function(x){
#     if (x<=0.0001) x="****"
#     else if (x<=0.001)x="***" 
#     else if (x<=0.01) x="**"
#     else if (x<=0.05) x="*"
#     else if (x>0.05) return ("ns")
#     return (x)
#   })
#   
#   selected_ensm_join_list[[i]]<-selected_ensm_join
#   
#   #selected_ensm_join_plot<-selected_ensm_join[1:10,]
#   
#   # fname=paste0(outFolder,cell_types[i],"_barplot_hubgenes.pdf")
#   # pdf(fname,width=12,height=6)
#   # ggplot(data=selected_ensm_join_plot, aes(x=gene_name, y=log2FoldChange)) +
#   #   geom_bar(stat="identity",fill="blue")+
#   #   geom_text(aes(label=padj_map), vjust=1.6, color="white", size=3.5)+
#   #   theme_bw()+
#   #   #scale_fill_manual("legend", values = c("0_Stromal-1"="#DF7D99" ,"1_Macrophage-2"="#838EDF" ,"2_Macrophage-1"="#4E65A6","3_Endothelial-1"="#FFC000" ,"4_Monocyte"="#2BA3D3", "5_CD4_T-cell"="#9ABF5C" ,"6_Decidual"="#D14357" ,"7_CD8_T-cell"="#329B2D","8_LED"="#D5438E","9_Stromal-2"="#ED4315" ,"10_ILC"="#76956C" ,"11_NK-cell"="#7BC791","12_Smooth muscle cells-1"="#CA8588" ,"13_Myofibroblast"="#F88091" , "14_Macrophage-3"="#72C6C8" ,"15_Endothelial-2"="#E4652C" ,"16_DC"="#9B91B9" ,"17_Smooth muscle cells-2"="#A37584" ,"18_EVT"="2C3E18" ,"19_Plasmablast"="#745B48" ,"20_Smooth muscle cells-3"="#AA5485" ,"21_Macrophage-4"="#4E747A","22_B-cell"="#C59A89","23_Unciliated Epithelial"="#C9C76F"))+
#   #   theme(axis.text.x = element_text(angle = 45, hjust=1))+
#   #   theme(legend.position="none")+
#   #   xlab("")+
#   #   ylab("LogFC")+
#   #   ggtitle()
#   # dev.off()
#   #par(mfrow=c(1,2))
#   
#   #all_shortest_paths(subnetwork,from=x,to=s,mode="out")
#   
#   #plot(subnetwork, vertex.size=hs*50, main="Hubs")
#   
#   #plot(subnetwork, vertex.size=as*30, main="Authorities")
#   
#  # mean_distance(subnetwork, directed=F)
#   
#   }


#save(selected_ensm_join_list,file=paste0(outFolder,"selected_ensm_join_list.RData"))

#selected_ensm_join_df<-do.call(rbind,selected_ensm_join_list)
#save(selected_ensm_join_df,file=paste0(outFolder,"selected_ensm_join_df.RData"))
#ap<-all_shortest_paths(full.graph,from=,to=,mode="out")#,algorithm=c("bellman-ford"))


# ppi<-read.delim("protein.links.v11.0.txt.gz",sep=" ")
# #ppi<-read.delim("protein.aliases.v11.0.txt.gz",sep=" ")
# 
# 
# ppi<-ppi[which(ppi[,"combined_score"]>=950),] #confidence score >=95
# ppiBackup<-ppi
# 
# ppi[,"protein1"]<-as.character(ppi[,"protein1"])
# ppi[,"protein2"]<-as.character(ppi[,"protein2"])