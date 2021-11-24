
####################################################################################################################

# Fig. 7. Comparative analysis of single-cell signatures from uterine tissues and bulk transcriptomics from maternal peripheral blood throughout gestation.

####################################################################################################################
  
library(lme4)
library(tidyverse)


outFolder<-"./13_Longitudinal_study/"
system(paste0("mkdir -p ", outFolder))


##########################################################
## cluster colors and labels 
##########################################################

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D",
                  "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48",
                  "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte",
                         "CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast",
                         "Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
names(cluster.Colors)<-paste0(c(0:23),"_",names(cluster.Colors))


## cluster labels 
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)


##########################################################
# load Data from Dr.Adi's team
##########################################################

load("CellularTransc260Longit.RData")

# aT1
# The data frame aT1 in the .RData file that Azam has holds the gene symbols and mapping to probesets (transcript clusters) which are the rows of the expression set.
# b (expression matrix)


##########################################################
# load myometrium single cell signatures 
##########################################################

SCGeneSets2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="avg_log2FC")]<-"avg_logFC"
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="gene")]<-"ENSG"
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="symbol")]<-"gene"
SCGeneSets<-SCGeneSets2

SCGeneSets$cluster<-clust2Names[as.character(SCGeneSets$cluster)]
SCGeneSet<-SCGeneSets



top20Sig <- SCGeneSet %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% ungroup()
table(top20Sig$cluster)
write.csv( top20Sig,file=paste0(outFolder,"top20Sig.csv"))


##########################################################
# analysis throughout gestation
##########################################################


signature_plots<-function(SCGeneSets,outFolder="./13_Longitudinal_study/")
{
  
  system(paste0("mkdir -p ", outFolder))
  pile=NULL
  for(a in unique(SCGeneSets$cluster)){
    tmp=SCGeneSets[SCGeneSets$cluster==a,]
    tmp=tmp[order(tmp$avg_logFC,decreasing=TRUE),]
    pile=rbind(pile,tmp[1:min(c(20,dim(tmp)[1])),c("gene","cluster")])
  }
  names(pile)<-c("Symbol","Type")
  genes=pile
  #genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("Tcells-activated","Tcells-resting"),"Symbol"],Type="Tcell"))
  #genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("5_CD4_T-cell" ,"7_CD8_T-cell"),"Symbol"],Type="Tcell"))
  
  celsig=list()
  for(cs in unique(genes$Type)){
    celsig[[cs]]<-intersect(as.character(genes$Symbol[genes$Type==cs]),aT1$SYMBOL)
  }
  
  myb=b[rownames(aT1),]
  
  esetM=apply(myb[,ano$GA>=37],1,mean)
  esetS=apply(myb[,ano$GA>=37],1,sd)
  eset2=apply(myb,2,function(x){(x-esetM)/esetS})
  #eset2=apply(myb,2,function(x){(x-esetM)})
  eset2=myb
  
  
  SYMBOLS<-aT1$SYMBOL
  
  toad=NULL
  names_celsig<-c()
  for (cs in names(celsig)){
    
    if (length(which(SYMBOLS%in%celsig[[cs]]))>1)
    {
      toad=rbind(toad,apply(as.matrix(eset2[SYMBOLS%in%celsig[[cs]],]),2,mean))
      names_celsig<-c(names_celsig,cs)
    }
    
  }
  rownames(toad)<-names_celsig #names(celsig)
  
  tg=data.frame(ID=rownames(toad),ID2=rownames(toad))
  
  f3=function(xx,tg,nm="_Sigs.pdf"){
    
    #pdf(paste("./",nm,sep=""))
    pdf(paste(outFolder,nm,sep=""))
    par(mfrow=c(3,3))
    for(an in tg$ID){
      ano2=ano
      ano2$Y=(xx[an,])
      
      
      
      mod1=lmer(Y~x+x2+x3+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      mod2=lmer(Y~1+(1|Main_Index),data = ano2,control=lmerControl(optimizer="bobyqa"),REML=FALSE)  
      
      p=anova(mod1,mod2)$"Pr(>Chisq)"[2]
      pred=expand.grid(GA=seq(10,40,by=0.1),Main_Index="newpoint")
      pred$x=pred$GA/40
      pred$x2=pred$x^2
      pred$x3=pred$x^3
      pred$Y=predict(mod1,pred,allow.new.levels=TRUE)
      FC=round(2^(max(pred$Y)-min(pred$Y)),1)
      
      #white 
      par(bg = "white")
      plot(0,0,xlim=c(9,40),cex.lab=1.2,cex.axis=1.5,ylim=range(ano2$Y),
           main=paste(an," p=",round(p,3),"Fold Change=",FC),ylab="Average log2 Intensity",
           xlab="Gestational age (weeks)",cex.main=0.7)  
      for(s in unique(ano2$Main_Index)){
        tmp=ano2[ano2$Main_Index==s,]
        tmp=tmp[order(tmp$GA),]
        points(Y~GA,tmp,type="l",col="darkgrey",lty=1)
        
        points(Y~GA,tmp,col="darkgrey")
        
      }
      points(Y~GA,data=pred[pred$GA<=40&pred$GA>=10,],lwd=3,type="l",col="blue",lty=1)
      
    }
    dev.off()
    
  }
  
  #f3(xx=toad,tg=tg[c(13,11,7),],nm="Sigs_sel.pdf")
  f3(xx=toad,tg=tg,nm="Sigs_sel.pdf")
  
}


##########################################################
# Make the signature analysis plots
##########################################################

signature_plots(SCGeneSets=SCGeneSet,outFolder="./13_Longitudinal_study/")




# script ends here
##############################################################################################################################################################################





