#Analysis work-flow for the article Targeted expression profiling by RNA-Seq improves detection of cellular dynamics during pregnancy and identifies a role for T cells in term parturition
#Targeted expression profiling by RNA-Seq improves detection of cellular dynamics during pregnancy and identifies a role for T cells in term parturition

#http://bioinformaticsprb.osf-demo.com/wp-content/uploads/2019/01/compPregnomics.pdf
library(pregnomics)
library(hta20sttranscriptcluster.db) #hta20sttranscriptcluster.db_8.3.1 
library(org.Hs.eg.db) #org.Hs.eg.db_3.2.3
library(hta20transcriptcluster.db)
library(EnsDb.Hsapiens.v75) # EnsDb.Hsapiens.v75_2.99.0 
library(annotate) #annotate_1.48.0
library(limma) 
library(DESeq2) 
library(UpSetR) 
#library(epiR) 
library(pROC) 
library(ROCR) 
library(gplots) 
#library(Heatplus) 
#library(marray) 
library(lme4) 
library(splines)
library(tidyverse)

outFolder="reference_Adi/"
#outFolder="reference_Adi/up-regulated/"
#outFolder="reference_Adi/down-regulated/"
system(paste0("mkdir -p ", outFolder))

############# 
data(ano32)
ano32$T=factor(ifelse(ano32$GA<37,"Preterm","Term")) #define gestational age interval anoALL<-ano32
anoALL<-ano32
head(anoALL,n=3)

anoGA=anoALL[anoALL$GAAnalysis==1,]
plot(as.numeric(as.factor(IndividualID))~GA,data=anoGA,col=ifelse(anoGA$Group=="TIL","red","black"),pch=19,xlab="GA",ylab="Patient")
legend("topleft",legend=c("TIL","TNL"),pch=c(19,19),col=c("red","black"))
abline(v=37,lty=2,col="grey")

anoLabor=anoALL[anoALL$LaborAnalysis==1,]
plot(as.numeric(as.factor(IndividualID))~GA,data=anoLabor,col=ifelse(anoLabor$Group=="TIL","red","black"),pch=19,xlab="GA",ylab="Patient",xlim=range(ano32$GA))
legend("topleft",legend=c("TIL","TNL"),pch=c(19,19),col=c("red","black"))


############################ 
## load data
############################
data(package="pregnomics",list=c("esetHTA","Rcount","Ccount","esetPCR"))

########################################################
# Comparison with other signatures
########################################################

################################
#microarray and qRT-PCR
################################
analyzeGA_limma=function(ano,eset){
  ano$ID=factor(ano$IndividualID)
  design <- model.matrix(~0+T+IndividualID,ano) 
  eset=eset[,rownames(ano)] 
  colnames(design)<-substr(colnames(design),2,100)
  fit <- lmFit(eset, design)
  cont.matrix <- makeContrasts( contrasts="Term-Preterm",levels=design) 
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  aT1<-topTable(fit2,coef=1, number=1000000, adjust="fdr") 
  aT1$FC=2^abs(aT1$logFC)*sign(aT1$logFC) #signed linear fold change 
  aT1$ID=rownames(aT1)
  
  aT1
}
analyzeLabor_limma=function(ano,eset){
  design <- model.matrix(~0+Group,ano) 
  colnames(design)<-gsub("Group","",colnames(design))
  fit <- lmFit(eset, design)
  cont.matrix <- makeContrasts( contrasts="TIL-TNL",levels=design) 
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  aT1<-topTable(fit2,coef=1, number=1000000, adjust="fdr") 
  aT1$FC=2^abs(aT1$logFC)*sign(aT1$logFC) #signed linear fold change 
  aT1$ID=rownames(aT1)
  aT1
}

#HTA microarray data

anpack="hta20sttranscriptcluster"
HTA=list()

################################
#HTA
################################
# aT1=analyzeGA_limma(anoGA,esetHTA[,rownames(anoGA)])
# head(aT1,n=3)
# aT1$SYMBOL<-unlist((lookUp(aT1$ID, anpack, 'SYMBOL')))
# aT1$ENTREZ<-unlist((lookUp(aT1$ID, anpack, 'ENTREZID')))
# aT1=aT1[!is.na(aT1$ENTREZ),]
# data(npspge) # based on detection above background from Affymetrix Transcriptome Analysis Consol e
# expressed=names(npspge)[npspge>0]
# aT1=aT1[rownames(aT1)%in%expressed,]
# aT1$adj.P.Val=p.adjust(aT1$P.Value,"fdr") # HTA[["GAEffect"]]<-aT1
# head(aT1,n=3)
# 
# aT1=analyzeLabor_limma(anoLabor,esetHTA[,rownames(anoLabor)]) 
# aT1$SYMBOL<-unlist((lookUp(aT1$ID, anpack, 'SYMBOL')))
# aT1$ENTREZ<-unlist((lookUp(aT1$ID, anpack, 'ENTREZID')))
# aT1=aT1[!is.na(aT1$ENTREZ),]
# aT1=aT1[rownames(aT1)%in%expressed,]
# aT1$adj.P.Val=p.adjust(aT1$P.Value,"fdr")
# HTA[["LaborEffect"]]<-aT1
# head(aT1,n=3)

################################
## PCR
################################
PCR<-list()
aT1=analyzeGA_limma(anoGA,esetPCR[,rownames(anoGA)])
aT1$SYMBOL=rownames(aT1)
aT1$Sig=(aT1$P.Value<0.05)
PCR[["GAEffect"]]<-aT1
head(aT1,n=3)

aT1=analyzeLabor_limma(anoLabor,esetPCR[,rownames(anoLabor)])
aT1$SYMBOL=rownames(aT1)
aT1$Sig=(aT1$P.Value<0.05)
PCR[["LaborEffect"]]<-aT1
head(aT1,n=3)

saveRDS(PCR,paste0(outFolder,"PCR.rds"))

################################
### sequencing platforms
################################

analyzeGA_DESeq=function(ano,anoall,countM){
  dds<- DESeqDataSetFromMatrix(countData= countM[,rownames(ano)],colData= ano,design=~T+IndividualID)
  dds<- DESeq(dds) 
  res<-results(dds,contrast=c("T","Term","Preterm"),independentFiltering=FALSE) 
  res=as.data.frame(res) 
  expressed=rownames(countM)[apply(countM[,rownames(anoall)]>=5,1,sum)>5] 
  res=res[rownames(res)%in%expressed,]
  res=res[!is.na(res$log2FoldChange),]
  res$logFC=res$log2FoldChange
  names(res)[names(res)=="pvalue"]<-"P.Value" 
  names(res)[names(res)=="padj"]<-"adj.P.Val"
  res
}

analyzeLabor_DESeq=function(ano,anoall,countM){
  dds<- DESeqDataSetFromMatrix(countData= countM[,rownames(ano)],colData= ano,design=~Group) 
  dds<- DESeq(dds) 
  res<-results(dds,contrast=c("Group","TIL","TNL"),independentFiltering=FALSE) 
  res=as.data.frame(res) 
  expressed=rownames(countM)[apply(countM[,rownames(anoall)]>=5,1,sum)>5] 
  res=res[rownames(res)%in%expressed,]
  res=res[!is.na(res$log2FoldChange),]
  res$logFC=res$log2FoldChange
  names(res)[names(res)=="pvalue"]<-"P.Value"
  names(res)[names(res)=="padj"]<-"adj.P.Val"
  res
}


edb <- EnsDb.Hsapiens.v75
Tx.ensemble <- transcripts(edb, columns = c("tx_id", "gene_id", "gene_name"),
                           return.type = "DataFrame")


################################
#RNASeq
################################
RNASeq<-list()

#GA effect
res=analyzeGA_DESeq(ano=anoGA,anoall=anoALL,countM=Rcount)
res$SYMBOL=Tx.ensemble[match(rownames(res),Tx.ensemble[,2]),3]
RNASeq[["GAEffect"]]<-res
#labor effect
res=analyzeLabor_DESeq(ano=anoLabor,anoall=anoALL,countM=Rcount)
res$SYMBOL=Tx.ensemble[match(rownames(res),Tx.ensemble[,2]),3]
RNASeq[["LaborEffect"]]<-res
head(res,n=3)

saveRDS(RNASeq,paste0(outFolder,"RNASeq.rds"))

################################
#DriverMap
################################
CELLECTA=list()
#GA effect
anoall=anoALL
anoall=anoall[anoall$SampleID!="Sample_26",] #remove the contaminated sample 
ano=anoGA
ano=ano[ano$SampleID!="Sample_26",] 
res=analyzeGA_DESeq(ano,anoall,countM=Ccount)
res$SYMBOL=rownames(res)
CELLECTA[["GAEffect"]]<-res

#labor effect
res=analyzeLabor_DESeq(ano=anoLabor,anoall,countM=Ccount)
res$SYMBOL=rownames(res)
CELLECTA[["LaborEffect"]]<-res
head(res,n=3)

saveRDS(CELLECTA,paste0(outFolder,"CELLECTA.rds"))


####################################################################
# Analysis of gene set signature expression
####################################################################
data(package="pregnomics",list=c("esetHTA","Rcount","Ccount","esetPCR"))

#outFolder="reference_Adi/"
CELLECTA<-read_rds("reference_Adi/CELLECTA.rds")
RNASeq<-read_rds("reference_Adi/RNASeq.rds")
PCR<-read_rds("reference_Adi/PCR.rds")

res <-read_tsv(paste0("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv"))
res$cluster<-unlist(strsplit(res$cname,"_"))[seq(1,2*length(res$cname),by=2)]
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
               "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
               "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(clust2Names)<-as.character(c(0:23))
res$Type<-clust2Names[res$cluster]
#res<-res  %>% dplyr::filter(log2FoldChange<0)

res<-res  %>% filter(log2FoldChange<=0)

SCGeneSets<-res %>% dplyr::select (gene_name,Type)
colnames(SCGeneSets)[1]<-"Symbol"

#data(SCGeneSets) 
#SCGeneSets=SCGeneSets[SCGeneSets$Symbol%in%comg,] 
#head(SCGeneSets)
table(SCGeneSets$Type)



data(ano32)
ano32$T=factor(ifelse(ano32$GA<37,"Preterm","Term")) #define gestational age interval anoALL<-ano32
anoALL<-ano32

Pr=esetPCR[,rownames(anoALL)]
dds<- DESeqDataSetFromMatrix(countData= Rcount[,rownames(anoALL)],colData= anoALL,design=~IndividualID)
dds=estimateSizeFactors(dds)
Rr=counts(dds, normalized=TRUE)
a=RNASeq[[2]] 
Rr=Rr[rownames(Rr)%in%rownames(a),] 
Rr=Rr[rownames(a),] 
Rr=Rr[!duplicated(a$SYMBOL),] 
rownames(Rr)=a[rownames(Rr),"SYMBOL"]
dds<- DESeqDataSetFromMatrix(countData= Ccount[,rownames(anoALL)],colData= anoALL,design=~IndividualID)
dds=estimateSizeFactors(dds)
Cr=counts(dds, normalized=TRUE)



#comg=c(rownames(Rr),rownames(Cr))
#comg=names(table(comg)[table(comg)==2])

#selection based on our genes:
#comg<-unique(res$gene_name)



nms=c("RNA-Seq","DriverMap")
ys=c("R","C")
allsig<-unique(SCGeneSets$Type)
for (sig in allsig)
{
  #sig="T cell"
  anoLabor=anoALL[anoALL$LaborAnalysis==1,]
  ano=anoLabor
  ano$Group0=factor(ifelse(ano$Group=="TIL","TIL","_TNL"))
  #ano$H=apply(Hr[SCGeneSets[SCGeneSets$Type==sig,"Symbol"],ano$SampleID],2,mean)
  gs<- SCGeneSets$Symbol[SCGeneSets$Type==sig]
   
  gsR<-gs[which(gs %in% rownames(Rr))]
  gsC<-gs[which(gs %in% rownames(Cr))]
  if (length(gsR)>1 & length(gsC)>1)
  {
   
    ano$R=apply(Rr[gsR,ano$SampleID],2,mean)
    ano$C=apply(Cr[gsC,ano$SampleID],2,mean)
    
    for(meths in 1:2){
      ano$Y=ano[,ys[meths]]
      lgFC=mean(ano$Y[ano$Group=="TIL"])-mean(ano$Y[ano$Group=="TNL"]) 
      pv=t.test(ano$Y[ano$Group=="TIL"],ano$Y[ano$Group=="TNL"])$p.value 
      cat(pv)
      #fname=paste0(paste0(outFolder,"down-regulated/"),sig, "_",nms[meths],"_boxplot.pdf")
      fname=paste0(paste0(outFolder,"up-regulated/"),sig, "_",nms[meths],"_boxplot.pdf")
      
      pdf(fname,width=10,height=6)
      #boxplot(Y~Group0,ano,ylab=sig,ylim=c(min(ano$Y)-0.9,max(ano$Y)),main=nms[meths],cex.axis=1.2) 
      
      
      TIL<-ano$Y[ano$Group=="TIL"]
      TNL<-ano$Y[ano$Group=="TNL"]
      
      data<-cbind(TIL, TNL)
      #boxplot(data,ylab=sig,ylim=c(min(data)-0.9,max(data)),main=nms[meths],cex.axis=1.2) 
      
      boxplot(Y~Group0,ano,ylab=sig,ylim=c(min(ano$Y)-0.9,max(ano$Y)),main=nms[meths],cex.axis=1.2) 
      
      legend("bottomleft",c(paste("log2FC=",round(lgFC,2)),paste("p=",round(pv,4))),cex=0.9,bty="n") 
      dev.off()
      } 
  }
  
}

####################################################
# gene boxplot 
####################################################
outFolder2<-"./8_outputs_DESeq_Plots/RNASeq/"
intersected_genes<-read.csv(paste0(outFolder2,"intersected_genes.csv"),stringsAsFactors = FALSE)
allsig<-unique(intersected_genes$gene_name)
#allsig<-unique(SCGeneSets$Type)
meths<-1 #"RNA-Seq"

outFolder2<-"./8_outputs_DESeq_Plots/DriverMap/"
intersected_genes<-read.csv(paste0(outFolder2,"intersected_genes.csv"),stringsAsFactors = FALSE)
allsig<-unique(intersected_genes$gene_name)
meths<-2

nms=c("RNA-Seq","DriverMap")
ys=c("R","C")

for (sig in allsig)
{
  #sig="T cell"
  anoLabor=anoALL[anoALL$LaborAnalysis==1,]
  ano=anoLabor
  ano$Group0=factor(ifelse(ano$Group=="TIL","TIL","_TNL"))
  #ano$H=apply(Hr[SCGeneSets[SCGeneSets$Type==sig,"Symbol"],ano$SampleID],2,mean)
  gs<- sig #SCGeneSets$Symbol[SCGeneSets$Type==sig]
  
  if (meths==1)
    gsR<-gs[which(gs %in% rownames(Rr))]
  if (meths==2)
    gsC<-gs[which(gs %in% rownames(Cr))]
  if (length(gsR)>0 & meths==1 || length(gsC)>0 & meths==2)
  {
    
    if (meths==1)
     ano$R=Rr[gsR,ano$SampleID] #apply(as.matrix(Rr[gsR,ano$SampleID]),2,mean)
     
    if (meths==2)
     ano$C=Cr[gsC,ano$SampleID] #apply(as.matrix(Cr[gsC,ano$SampleID]),2,mean)
    
      ano$Y=ano[,ys[meths]]
      lgFC=mean(ano$Y[ano$Group=="TIL"])-mean(ano$Y[ano$Group=="TNL"]) 
      pv=t.test(ano$Y[ano$Group=="TIL"],ano$Y[ano$Group=="TNL"])$p.value 
      print(pv)
      system(paste0("mkdir -p ", paste0(outFolder,nms[meths],"/")))
      fname=paste0(outFolder,nms[meths],"/",sig, "_boxplot.pdf")
      #fname=paste0(paste0(outFolder,"up-regulated/"),sig, "_",nms[meths],"_boxplot.pdf")
      
      pdf(fname,width=10,height=6)
      
      TIL<-ano$Y[ano$Group=="TIL"]
      TNL<-ano$Y[ano$Group=="TNL"]
      
      data<-cbind(TIL, TNL)
      #boxplot(data,ylab=sig,ylim=c(min(data)-0.9,max(data)),main=nms[meths],cex.axis=1.2) 
      
      boxplot(Y~Group0,ano,ylab=sig,ylim=c(min(ano$Y)-0.9,max(ano$Y)),main=nms[meths],cex.axis=1.2) 
      
      #if (lgFC>0)
      #legend("bottomleft",c(paste("log2FC=",round(lgFC,2)),paste("p=",round(pv,4))),cex=0.9,bty="n") 
      #else 
      legend("topleft",c(paste("log2FC=",round(lgFC,2)),paste("p=",round(pv,4))),cex=0.9,bty="n") 
      
      dev.off()
    
  }
  
}
