library(lme4)
#load expression data


outFolder="./13_Longitudinal_study/"
system(paste0("mkdir -p ", outFolder))



load("CellularTransc260Longit.RData")




# read myometrium single  cell sig (gene sets)
SCGeneSets2 = read_tsv("5_harmonyClustersDGE/ClusterDEG.tsv")
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="avg_log2FC")]<-"avg_logFC"
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="gene")]<-"ENSG"
colnames(SCGeneSets2)[which(colnames(SCGeneSets2)=="symbol")]<-"gene"
SCGeneSets<-SCGeneSets2

clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)
SCGeneSets$cluster<-clust2Names[as.character(SCGeneSets$cluster)]




#read gene sets
#SCGeneSets=read.table("MarkersAll.txt",sep=" ",header=TRUE,stringsAsFactors = FALSE)


pile=NULL
for(a in unique(SCGeneSets$cluster)){
  tmp=SCGeneSets[SCGeneSets$cluster==a,]
  tmp=tmp[order(tmp$avg_logFC,decreasing=TRUE),]
  pile=rbind(pile,tmp[1:min(c(20,dim(tmp)[1])),c("gene","cluster")])
}
names(pile)<-c("Symbol","Type")
genes=pile
#genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("Tcells-activated","Tcells-resting"),"Symbol"],Type="Tcell"))
genes=rbind(genes,data.frame(Symbol=genes[genes$Type%in%c("5_CD4_T-cell" ,"7_CD8_T-cell"),"Symbol"],Type="Tcell"))


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
for (cs in names(celsig)){
  toad=rbind(toad,apply(eset2[SYMBOLS%in%celsig[[cs]],],2,mean))
}
rownames(toad)<-names(celsig)

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


