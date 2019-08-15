### analyse CANDIDATES###
library(data.table)
library(dplyr)
for (pair in c("INECARHCADRG","INECAROBIGUN","INECARTKOHRA","INECARWILKAS","INECARZEPSUB","LACTISHCADRG","LACTISINECAR","LACTISOBIGUN","LACTISTKOHRA","LACTISWILKAS","LACTISZEPSUB","OBIGUNHCADRG","TKOHRAHCADRG","TKOHRAOBIGUN","TKOHRAWILKAS","TKOHRAZEPSUB","WILKASOBIGUN","WILKASZEPSUB","ZEPSUBHCADRG","ZEPSUBOBIGUN")) { ###
  setwd(paste('/home/aa/alpine/dmc/analysis_5/',pair,sep="")) ###
  i2<-dir(path = ".",include.dirs = F,recursive = F,pattern = "AL.G*") 
  for (id in i2){ #   id = "AL3G47770"
    setwd(paste("/home/aa/alpine/dmc/analysis_5/",pair,"/",id,sep="")) ###
    ar<-read.table(paste('dmcInput.txt',sep=''),sep = " ")
    colnames(ar) <- c("START","ID","ANN","AAS","p1","p2","p3","p4")
    
    ar1<-subset(ar,ar$ID %in% i2)
    left<-min(ar1$START)
    right<-max(ar1$START)
    
    ar$afd1<-abs(ar$p1-ar$p2)
    ar$afd2<-abs(ar$p3-ar$p4)
    ar$afd3<-abs(ar$p1-ar$p3)
    
    pdf(paste("AFD_dotplot.pdf",sep=""),height = 11,width = 11)
    par(oma=c(2,3,3,0),mar=c(2.1, 4.1, 0.1, 1.1),mfrow=c(3,1))
    smoothingSpline = smooth.spline(subset(ar,ar$afd1>=(quantile(ar$afd1,probs = 0.70)))[,1], subset(ar,ar$afd1>=(quantile(ar$afd1,probs = 0.70)))[,9], spar=0.3)
    plot(ar$afd1~ar$START,pch=1,ylab = "AFD - alpine-foothill 1",col="blue3",ylim=c(0,1.1),xlab="")
    lines(x = smoothingSpline$x,y = smoothingSpline$y,col="red")
    segments(x0 = left,y0 = 1.05,x1 = right,y1 =1.05,col = "grey50",lwd = 4)
    
    smoothingSpline = smooth.spline(subset(ar,ar$afd2>=(quantile(ar$afd2,probs = 0.70)))[,1], subset(ar,ar$afd2>=(quantile(ar$afd2,probs = 0.70)))[,10], spar=0.3)
    plot(ar$afd2~ar$START,pch=1,ylab = "AFD - alpine-foothill 2",col="blue3",ylim=c(0,1.1),xlab="")
    lines(x = smoothingSpline$x,y = smoothingSpline$y,col="red")
    segments(x0 = left,y0 = 1.05,x1 = right,y1 =1.05,col = "grey50",lwd = 4)
    smoothingSpline = smooth.spline(subset(ar,ar$afd3>=(quantile(ar$afd3,probs = 0.70)))[,1], subset(ar,ar$afd3>=(quantile(ar$afd3,probs = 0.70)))[,11], spar=0.3)
    plot(ar$afd3~ar$START,pch=1,ylab = "AFD - alpine-alpine",col="blue3",ylim=c(0,1.1),xlab="")
    lines(x = smoothingSpline$x,y = smoothingSpline$y,col="red")
    
    segments(x0 = left,y0 = 1.05,x1 = right,y1 =1.05,col = "grey50",lwd = 4)
    
    dev.off()
    
  }
}






  
  
   
  
#  points(x=td$start, y=td$TDI, pch=19, col="red", bg=NA, cex=1.2)
#  text(tdm$start,tdm$TDI,tdm$aas,cex=0.35, pos=1,col="blue3",srt=90,offset = 0.7)
#  plot(s4$TTE~s4$start,pch=1,ylab = "Tatry4x",col="grey50",ylim=c(0,1))
#  points(x=tt$start, y=tt$TTE, pch=19, col="red", bg=NA, cex=1.2)
#  text(ttm$start,ttm$TTE,ttm$aas,cex=0.35, pos=1,col="blue3",srt=90,offset = 1)
#  plot(s4$FAG~s4$start,pch=1,ylab = "Fagaras",col="grey50",ylim=c(0,1))
#  points(x=fa$start, y=fa$FAG, pch=19, col="red", bg=NA, cex=1.2)
#  text(fam$start,fam$FAG,fam$aas,cex=0.35, pos=1,col="blue3",srt=90,offset = 1)
#  plot(s4$ROD~s4$start,pch=1,ylab = "Rodna",col="grey50",ylim=c(0,1))
#  points(x=ro$start, y=ro$ROD, pch=19, col="red", bg=NA, cex=1.2)
#  text(rom$start,rom$ROD,rom$aas,cex=0.35, pos=1,col="blue3",srt=90,offset = 1)
#  plot(s4$ALP~s4$start,pch=1,ylab = "Alps",col="grey50",ylim=c(0,1))
#  points(x=al$start, y=al$ALP, pch=19, col="red", bg=NA, cex=1.2)
#  text(alm$start,alm$ALP,alm$aas,cex=0.35, pos=1,col="blue3",srt=90,offset = 1)
#  mtext(paste(i2,ann1$AT,ann1$annShort,sep=" - "), side = 3, line = 1, cex = 1.1,outer=TRUE)
#  mtext("Allele frequency difference between foothill and alpine populations", side = 2, line = 1, cex = 1,outer=TRUE)
#  mtext("Position", side = 1, line = 1, cex = 1,outer=TRUE)
#  #   dev.off()
#  
#  
#  ###Heatmap
#  s<-select(tot,ID,ann,aas,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,KAS,HOC,tot)
#  #repolarise
#  for (ii in  1:nrow(s)){ # i=1
#    if (s$tot[ii]>0.5)
#    {s[ii,4:21]<-1-s[ii,4:21]
#    } else {}}
#  #4. plot
#  ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt")
#  id<-i2 #   id = "AL3G47770"
#  if (nrow(subset(s,s$ID %in% id))>1)
#  {pops<-c("VEL", 'ZEP', 'TKO', 'TRT', 'LAC', 'BAL', 'INE','SCH','WIL','SUB', 'BAB','HRA','SPI','DRA','TIS','CAR','KAS','HOC')
#  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
#  s1<-subset(s,s$ID %in% id)
#  ann1<-subset(ann,ann$AL %in% id)
#  
#  if (nrow(s1)>110)
#  {p=0.75
#  } else {p=1}
#  #  df<-as.matrix(s1[,4:25],rownames = paste(s1$i.ann,s1$i.aas,sep="   "))
#  df<-as.matrix(s1[,4:21],rownames = paste(s1$aas,sep="   "))
#  # pdf(paste("heatmap_",id,".pdf",sep=""),height = 11,width = 9.3)
#  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,4,6#,7,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,11,13,15,16,18),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c#("orange","orange","orange","orange","orange","orange","orange","orange","orange","green","green","green","green","green","green"#,"green","green","green"),labCol=pops,colCol= c("orange","orange","orange","orange","orange","orange","orange","orange","orange"#,"green","green","green","green","green","green","green","green","green"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p#),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - "))
#  
#  } else {}
#  dev.off()
#  
#}
