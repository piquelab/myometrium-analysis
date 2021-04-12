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

res <- res %>% filter(!is.na(log2FoldChange))

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

# nm[which(node_size<0.5)]<-""
# myfill<-rep("yellow",length(node_size))
# myfill[which(node_size>0.5)]<-"orange"

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
ggraph(subnetwork,layout = "kk") +  #fr kk
  geom_edge_link(colour = "gray",show.legend = FALSE,alpha=0.2) +
  geom_node_point(aes(size = node_size),color=myfill,shape=16,show.legend = FALSE) + 
  #scale_size(range = c(20, 30)) +
  #scale_color()+
 # theme_graph()+
  #scale_size(range = c(1,6))+
  #geom_node_point(aes(size = centrality, colour = centrality)) + 
geom_node_text(aes(label = nm), colour = 'black', vjust = -0.8)

#scale_fill_manual(values =cluster.Colors[unique(res$Cell_type)], labels = as.character(unique(res$Cell_type)))+
#theme_graph()


##################################################################################
# systematic way of selecting genes to plot 
##################################################################################

gene_pvalue_membership<-matrix(-1,nrow=length(unique(res$gene_name)),ncol=length(cell_types))
rownames(gene_pvalue_membership)<-unique(res$gene_name)
colnames(gene_pvalue_membership)<-cell_types

gene_logfc_membership<-matrix(0,nrow=length(unique(res$gene_name)),ncol=length(cell_types))
rownames(gene_logfc_membership)<-unique(res$gene_name)
colnames(gene_logfc_membership)<-cell_types


gene_membership<-matrix(0,nrow=length(unique(res$gene_name)),ncol=length(cell_types))
rownames(gene_membership)<-unique(res$gene_name)
colnames(gene_membership)<-cell_types


selectedgenes_per_celltype<-array()

for ( i in 1:length(cell_types))
{
  print(i)
  
  #rank (positive)
  res_selected_ensm<-res %>% filter(Cell_type ==cell_types[i] & padj <= 0.1)# & padj <= 0.05 & log2FoldChange>2)# 
  res_celltype<-res %>% filter(Cell_type ==cell_types[i])
  res_celltype<-res_celltype[order(abs(res_celltype$log2FoldChange),-res_celltype$padj, decreasing = TRUE),]
  
  
  res_selected_ensm<-res_selected_ensm[order(abs(res_selected_ensm$log2FoldChange),-res_selected_ensm$padj, decreasing = TRUE),]
  selected_ensm<-res_selected_ensm$STRING_id
  selected_gene_name<-res_selected_ensm$gene_name
  
  gene_logfc_membership[selected_gene_name,cell_types[i]]<-res_selected_ensm$log2FoldChange [which(res_selected_ensm$gene_name %in% selected_gene_name)] 
  gene_pvalue_membership[selected_gene_name,cell_types[i]]<-res_selected_ensm$padj [which(res_selected_ensm$gene_name %in% selected_gene_name)] 
  gene_membership[selected_gene_name,cell_types[i]]<-c(1:length(selected_gene_name))
  
  
  #rank negative
  subnetwork<-string_db$get_subnetwork(res_celltype$STRING_id[res_celltype$Cell_type==cell_types[i]])
  #V(subnetwork)$name<-names(aliases_map)[which(aliases_map %in% V(subnetwork)$name)]
  hs <- hub_score(subnetwork, weights=NA)$vector
  hs<-hs[order(hs,decreasing = TRUE)]
  selected_hubs<-names(aliases_map)[which(aliases_map %in% names(hs))][1:20]
 
  gene_membership[selected_hubs,cell_types[i]]<--1*c(1:length(selected_hubs))
  gene_logfc_membership[selected_hubs,cell_types[i]]<-res_celltype$log2FoldChange [which(res_celltype$gene_name %in% selected_hubs)] 
  gene_pvalue_membership[selected_hubs,cell_types[i]]<-res_celltype$padj [which(res_celltype$gene_name %in% selected_hubs)] 
  
  sl<-c(selected_ensm[1:10],names(hs)[1:20])
  
  selectedgenes_per_celltype<-cbind(selectedgenes_per_celltype,sl)
}

colnames(gene_pvalue_membership)<-cell_types
colnames(gene_logfc_membership)<-cell_types
colnames(gene_membership)<-cell_types
rownames(gene_pvalue_membership)<-unique(res$gene_name)
rownames(gene_logfc_membership)<-unique(res$gene_name)
rownames(gene_membership)<-unique(res$gene_name)


selectedgenes_per_celltype<-selectedgenes_per_celltype[,-1]
colnames(selectedgenes_per_celltype)<-cell_types

selected_ensm<-unique(c(selectedgenes_per_celltype))

remove<-aliases_map[which(names(aliases_map) %in% c("GUCY1A3", "GPR97",   "GPR126" ,"GIG25" ))]
selected_ensm<-selected_ensm[!selected_ensm %in% remove]

subnetwork<-string_db$get_subnetwork(selected_ensm)



V(subnetwork)$name<-names(aliases_map)[which(aliases_map %in% V(subnetwork)$name)]
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






cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")

#names(cluster.Colors)<-paste0(c(0:23),"_",names(cluster.Colors))




myfill<-rep("#B0B0B0",length(node_size))
names(myfill)<-names(node_size)
#myfill[which(node_size>0.5)]<-"orange"


summary_membership<-matrix(nrow=length(node_size),ncol=1)
rownames(summary_membership)<-names(node_size)
hub_genes<-c()
log_fc_genes<-c()
names_fc<-c()
for(x in names(node_size))
  {
  
#  gene_pvalue_membership
 # gene_logfc_membership
  if( x %in% rownames(gene_membership))
  {
    celltype_exits<-names(which(gene_membership[x,]!=0))
    
    
    if(length(celltype_exits)==1) 
    {
      
      fc_celltype_exits<-gene_logfc_membership[x,celltype_exits]
      names(fc_celltype_exits)<-celltype_exits
      pvalue_celltype_exits<-gene_pvalue_membership[x,celltype_exits]
      #fc_celltype_exits<-fc_celltype_exits[which(pvalue_celltype_exits<0.1)]
      log_fc_genes<-c(log_fc_genes,as.numeric(fc_celltype_exits))
      myfill[x]<-cluster.Colors[celltype_exits]
      names_fc<-c(names_fc,x)
    }
    
    
    else if(length(celltype_exits)>1)
    {
      fc_celltype_exits<-gene_logfc_membership[x,celltype_exits]
      pvalue_celltype_exits<-gene_pvalue_membership[x,celltype_exits]
      #fc_celltype_exits<-fc_celltype_exits[which(pvalue_celltype_exits<0.1)]
      
      fc_celltype_exits<-fc_celltype_exits[order(abs(fc_celltype_exits),decreasing = TRUE)]
      summary_membership[x,1]<-paste(names(fc_celltype_exits),sep=", ",collapse = ", ")
      
      
      myfill[x]<-cluster.Colors[names(fc_celltype_exits)[1]]
      log_fc_genes<-c(log_fc_genes,as.numeric(fc_celltype_exits)[1])
      names_fc<-c(names_fc,x)
    }
  }
    
    
if( ! x %in% rownames(gene_membership) || length(celltype_exits)==0)
{
  hub_genes<-c(hub_genes,x)
  #which(res$gene_name ==x)
  fcs<-res$log2FoldChange[which(res$gene_name ==x)]
  adpvalues<-res$padj[which(res$gene_name ==x)]
  fcs<-fcs[!is.na(adpvalues)]
  adpvalues<-adpvalues[!is.na(adpvalues)]
  if(length(adpvalues)==0)
    log_fc_genes<-c(log_fc_genes,NA)
  else
    log_fc_genes<-c(log_fc_genes,fcs[which(adpvalues==min(adpvalues))])
  
  names_fc<-c(names_fc,x)
}
  
}


# choose the color of nodes


#logfc for mapping 
# color based on the celltypes
#more genes 


# graph_tbl <- subnetwork %>% 
#   as_tbl_graph() %>% 
#   activate(nodes) %>% 
#   mutate(degree  = centrality_degree()) 


Centrality<-node_size
## remove text for non-hub nodes
nm<-V(subnetwork)$name
nm[which(node_size<0.2)]<-""
########################################
# plot graph
########################################


# cell type colors
ggraph(subnetwork,layout = "kk") +  #fr  kk
  geom_edge_link(colour = "gray",show.legend = FALSE,alpha=0.2) +
  geom_node_point(aes(size = node_size),color=myfill,shape=16,show.legend = FALSE) + 
  #scale_size(range = c(20, 30)) +
  #scale_color()+
  # theme_graph()+
  #scale_size(range = c(1,6))+
  #geom_node_point(aes(size = centrality, colour = centrality)) + 
  geom_node_text(aes(label = nm), colour = 'black', vjust = -0.8,label.size = 0.85)+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),)
#+scale_color_manual(values =cluster.Colors[unique(res$Cell_type)], labels = as.character(unique(res$Cell_type)))
#+theme_graph()



# log2fc gradient
ggraph(subnetwork,layout = "kk") +  #fr  kk
  geom_edge_link(colour = "gray",show.legend = FALSE,alpha=0.2) +
  geom_node_point(aes(size = Centrality,color = log_fc_genes),shape=16,show.legend = TRUE,legend.title="centrality") +  #color=log_fc_genes
  
  scale_color_gradient(name = "log2FC",low = "blue", high = "red")+
  #scale_colour_gradient(low="red",high = "blue") + #limits=c(-4, 5)
  geom_node_text(aes(label = nm), colour = 'black', vjust = -0.8,label.size = 0.8)+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),)



myfill[which(myfill=="#B0B0B0")]<-"black"

#names(myfill)<-names(node_size)

# text color: cell type
# node color: log2fc
ggraph(subnetwork,layout = "fr") +  #fr  kk
  geom_edge_link(colour = "gray",show.legend = FALSE,alpha=0.2) +
  geom_node_point(aes(size = Centrality,color = log_fc_genes),shape=16,show.legend = TRUE,legend.title="centrality",stroke=2) +  #color=log_fc_genes
  
  scale_color_gradient2(name = "log2FC",low = "#333399", high = "#A50021",mid="white")+ 
  #scale_colour_gradient(low="red",high = "blue") + #limits=c(-4, 5)
  geom_node_text(aes(label = nm), colour = myfill, vjust = -1,label.size = 0.8)+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),)





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