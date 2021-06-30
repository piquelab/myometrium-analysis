outFolder="./feature_importance/"

scale.data<-sc@assays$RNA@scale.data
colnames(scale.data)<-sc$seurat_clusters[colnames(scale.data)]
rw<-anno$gene_name[which(anno$kbid %in% rownames(scale.data))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(scale.data))]
scale.data<-scale.data[names(rw),]
rownames(scale.data)<-as.character(rw)
scale.data2<-t(scale.data)
scale.data2<-as.data.frame(scale.data2)
scale.data2$target<-colnames(scale.data)
write.csv(scale.data2,file=paste0(outFolder,"scale.data2.csv"))

data<-sc@assays$RNA@data
colnames(data)<-sc$seurat_clusters[colnames(data)]
rw<-anno$gene_name[which(anno$kbid %in% rownames(data))]
names(rw)<-anno$kbid[which(anno$kbid %in% rownames(data))]
data<-data[names(rw),]
rownames(data)<-as.character(rw)

data2<-t(as.matrix(data))
data2<-as.data.frame(data2)
data2$target<-colnames(data)
write.csv(data2,file=paste0(outFolder,"data2.csv"))
