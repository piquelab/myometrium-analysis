################################################################
# Quality metrics
###############################################################

# load merged sc object


outFolder <- paste0("1_quality_metrics/")
#outFolder <- paste0("11_pathway_enrichment_revision/")
system(paste0("mkdir -p ",outFolder))


fname=paste0("3_MergeDemux_Output/scFilteredSeurat.Rdata")
cat(fname,"\n")

#sc,anno
load(fname)


pdf(paste0(outFolder,"quality_metrics.pdf"),width=25,height=8)
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
x<-VlnPlot(sc, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

pdf(paste0(outFolder,"scatter_plot_quality_metrics.pdf"),width=10,height=8)
FeatureScatter(sc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()



#number of genes detected per cell
nFeature_RNA<-VlnPlot(sc, group.by = "orig.ident", features = "nFeature_RNA", pt.size = 0.1, ncol = 3) + NoLegend()

#number of UMIs per cell
nCount_RNA<-VlnPlot(sc, group.by = "orig.ident", features = "nCount_RNA", pt.size = 0.1, ncol = 3) + NoLegend()
percent.mt<-VlnPlot(sc, group.by = "orig.ident", features = "percent.mt", pt.size = 0.1, ncol = 3) + NoLegend()

nFeature_RNA_metric<-tapply(nFeature_RNA$data$nFeature_RNA, nFeature_RNA$data$ident,mean)
nCount_RNA_metric<-tapply(nCount_RNA$data$nCount_RNA, nCount_RNA$data$ident,mean)
percent.mt_metric<-tapply(percent.mt$data$percent.mt, percent.mt$data$ident,mean)

qc<-cbind(nFeature_RNA_metric,nCount_RNA_metric,percent.mt_metric)
colnames(qc)<-c("nFeature_RNA","nCount_RNA","percent.mt")
write.csv(qc,file=paste0(outFolder,"quality_metrics.csv"))

doner_cellcounts<-table(sc$Pregnancy_ID)
write.csv(doner_cellcounts,file="1_quality_metrics/doner_cellcounts.csv")
