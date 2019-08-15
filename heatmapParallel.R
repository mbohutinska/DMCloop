### VISUALIZE CANDIDATES###
### ID list in i2 <- tt2$Category[2:nrow(tt2)]
### ID list in i2 <- c("AL5G34820"),"AL5G34820")
library(data.table)
library(dplyr)
library(gplots)
library(RColorBrewer)

setwd("/home/aa/alpine/arenosaGenome/selScans/")
#1.extract genes from file
ar<-fread(paste('ann/ALLarenosa.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
#  ar<-fread(paste('ann/ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
ha<-fread(paste('ann/ALLhalleri.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ar),1))
colnames(ha) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ha),1))

for (pair in c("INECARHCADRG","INECAROBIGUN","INECARTKOHRA","INECARWILKAS","INECARZEPSUB","LACTISHCADRG","LACTISINECAR","LACTISOBIGUN","LACTISTKOHRA","LACTISWILKAS","LACTISZEPSUB","OBIGUNHCADRG","TKOHRAHCADRG","TKOHRAOBIGUN","TKOHRAWILKAS","TKOHRAZEPSUB","WILKASOBIGUN","WILKASZEPSUB","ZEPSUBHCADRG","ZEPSUBOBIGUN")) { ###
  setwd(paste('/home/aa/alpine/dmc/analysis_5/',pair,sep="")) ###
  i2<-dir(path = ".",include.dirs = F,recursive = F,pattern = "AL.G*") ###
  
  
ar1<-subset(ar,ar$ID %in% i2)
ha1<-subset(ha,ha$ID %in% i2) ###halleri

   #ar1<-subset(ar,ar$scaff %in% ha$scaff & ar$start %in% ha$start)
#write.table("all.txt",cbind(ar1,ha),col.names = F,row.names = F,quote = F)
#   ar2<-subset(ar,paste(ar$scaff,ar$start,sep="") %in% paste(ar1$scaff,ar1$start,sep=""))
#   ar1<-ar2

#2. filter too low freq
i<-which(colnames(ar1) %in% "10")
y<-which(colnames(ar1) %in% "151")
afh1<-ar1[,i:as.numeric(y)]
ar1$ACh1<-rowSums(afh1,na.rm = T)
ar1$NAh1<-apply(is.na(afh1), 1, sum)
ar1$ANh1<-504
ar1$tot<-ar1$ACh1/ar1$ANh1
ar2<-subset(ar1,ar1$tot > 24/504)
ar3<-subset(ar2,ar2$tot < as.numeric(1-(24/504)))

  #2. filter too low freq
  i<-which(colnames(ha1) %in% "10")
  y<-which(colnames(ha1) %in% "41")
  afh1<-ha1[,i:as.numeric(y)]
  ha1$ACh1<-rowSums(afh1,na.rm = T)
  ha1$NAh1<-apply(is.na(afh1), 1, sum)
  ha1$ANh1<-64
  ha1$tot<-ha1$ACh1/ha1$ANh1
  ha2<-subset(ha1,ha1$tot > 3.2/64)
  ha4<-subset(ha2,ha2$tot < as.numeric(1-(3.2/64)))

 # ar3<-ar1 ### DMC
ha3<-subset(ha,ha$ID %in% i2)
ar3$end<-ar3$start
ha3$end<-ha3$start
setkey(ha3, scaff, start,end)
a<-foverlaps(x = ar3, y = ha3, type="any")
#3. calculate AF for each lineage
i<-which(colnames(a) %in% "10")
afh1<-a[,i:as.numeric((i-1)+8)]
a$ACh1<-rowSums(afh1,na.rm = T)
a$NAh1<-apply(is.na(afh1), 1, sum)
a$ANh1<-(8-a$NAh1)*2
a$OBI<-a$ACh1/a$ANh1
afh2<-a[,as.numeric(i+8):as.numeric((i-1)+16)]
a$ACh2<-rowSums(afh2,na.rm = T)
a$NAh2<-apply(is.na(afh2), 1, sum)
a$ANh2<-(8-a$NAh2)*2
a$GUN<-a$ACh2/a$ANh2
aff1<-a[,as.numeric(i+16):as.numeric((i-1)+24)]
a$ACf1<-rowSums(aff1,na.rm = T)
a$NAf1<-apply(is.na(aff1), 1, sum)
a$ANf1<-(8-a$NAf1)*2
a$HCA<-a$ACf1/a$ANf1
aff2<-a[,as.numeric(i+24):as.numeric((i-1)+32)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*2
a$DRG<-a$ACf2/a$ANf2
i<-which(colnames(a) %in% "i.10")
afh1<-a[,i:as.numeric((i-1)+8)]
a$ACh1<-rowSums(afh1,na.rm = T)
a$NAh1<-apply(is.na(afh1), 1, sum)
a$ANh1<-(8-a$NAh1)*4
a$ING<-a$ACh1/a$ANh1
afh2<-a[,as.numeric(i+8):as.numeric((i-1)+16)]
a$ACh2<-rowSums(afh2,na.rm = T)
a$NAh2<-apply(is.na(afh2), 1, sum)
a$ANh2<-(8-a$NAh2)*2
a$VEL<-a$ACh2/a$ANh2
aff1<-a[,as.numeric(i+16):as.numeric((i-1)+24)]
a$ACf1<-rowSums(aff1,na.rm = T)
a$NAf1<-apply(is.na(aff1), 1, sum)
a$ANf1<-(8-a$NAf1)*2
a$ZEP<-a$ACf1/a$ANf1
aff2<-a[,as.numeric(i+24):as.numeric((i-1)+32)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*2
a$SUB<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+32):as.numeric((i-1)+40)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*2
a$BAB<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+40):as.numeric((i-1)+47)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(7-a$NAf2)*4
a$SCH<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+47):as.numeric((i-1)+55)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$WIL<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+55):as.numeric((i-1)+63)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$KAS<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+63):as.numeric((i-1)+71)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$LAC<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+71):as.numeric((i-1)+79)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$BAL<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+79):as.numeric((i-1)+87)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$DRA<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+87):as.numeric((i-1)+95)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$TIS<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+95):as.numeric((i-1)+102)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(7-a$NAf2)*4
a$INE<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+102):as.numeric((i-1)+110)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$CAR<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+110):as.numeric((i-1)+118)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$TKO<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+118):as.numeric((i-1)+126)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$TRT<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+126):as.numeric((i-1)+134)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$HRA<-a$ACf2/a$ANf2
aff2<-a[,as.numeric(i+134):as.numeric((i-1)+142)]
a$ACf2<-rowSums(aff2,na.rm = T)
a$NAf2<-apply(is.na(aff2), 1, sum)
a$ANf2<-(8-a$NAf2)*4
a$SPI<-a$ACf2/a$ANf2
s4<-select(a,i.ID,i.ann,i.aas,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG,tot)
s6<-select(a,i.ID,i.ann,i.aas,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,tot)

  i<-which(colnames(ha4) %in% "10")
  afh1<-ha4[,i:as.numeric((i-1)+8)]
  ha4$ACh1<-rowSums(afh1,na.rm = T)
  ha4$NAh1<-apply(is.na(afh1), 1, sum)
  ha4$ANh1<-(8-ha4$NAh1)*2
  ha4$OBI<-ha4$ACh1/ha4$ANh1
  afh2<-ha4[,as.numeric(i+8):as.numeric((i-1)+16)]
  ha4$ACh2<-rowSums(afh2,na.rm = T)
  ha4$NAh2<-apply(is.na(afh2), 1, sum)
  ha4$ANh2<-(8-ha4$NAh2)*2
  ha4$GUN<-ha4$ACh2/ha4$ANh2
  aff1<-ha4[,as.numeric(i+16):as.numeric((i-1)+24)]
  ha4$ACf1<-rowSums(aff1,na.rm = T)
  ha4$NAf1<-apply(is.na(aff1), 1, sum)
  ha4$ANf1<-(8-ha4$NAf1)*2
  ha4$HCA<-ha4$ACf1/ha4$ANf1
  aff2<-ha4[,as.numeric(i+24):as.numeric((i-1)+32)]
  ha4$ACf2<-rowSums(aff2,na.rm = T)
  ha4$NAf2<-apply(is.na(aff2), 1, sum)
  ha4$ANf2<-(8-ha4$NAf2)*2
  ha4$DRG<-ha4$ACf2/ha4$ANf2
  s5<-select(ha4,ID,ann,aas,OBI,GUN,HCA,DRG,tot)


# s4<-select(a,i.ID,i.ann,i.aas,i.end,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG,tot) ###DMC
# write.table (s4,"s4.txt",quote = F,col.names = T,row.names = F) ###DMC

### MORE POSSIBILITIES ###
s<-subset(s4,!s4$i.ann %in% "intragenic_variant" & !s4$i.ann %in% "downstream_gene_variant" & !s4$i.ann %in% "upstream_gene_variant")

s5<-subset(s5,!s5$ann %in% "intragenic_variant" & !s5$ann %in% "downstream_gene_variant" & !s5$ann %in% "upstream_gene_variant")

#s<-subset(s4,!s4$i.ann %in% "intragenic_variant")

#s<-subset(s4,s4$i.ann %like% "missense_variant")
#s<-subset(s4,s4$i.ann %like% "downstream_gene_variant")

# s<-subset(s4,s4$i.ann %like% "intragenic_variant")[0:2000,]
# s<-subset(s4,s4$i.ann %like% "intragenic_variant")[2001:4000,]
# s<-subset(s4,s4$i.ann %like% "intragenic_variant")[4001:6000,]
# s<-subset(s4,s4$i.ann %like% "intragenic_variant")[6001:9483,]

#s<-s4

#repolarise
for (i in  1:nrow(s)){ # i=1
  if (s$tot[i]>0.5)
  {s[i,4:25]<-1-s[i,4:25]
  } else {}}

#4. plot
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20171231.txt")
for (id in i2){ #   id = "AL3G47770"
  setwd(paste("/home/aa/alpine/dmc/analysis_5/",pair,"/",id,sep="")) ###
  if (nrow(subset(s,s$i.ID %in% id))>1)
  {pops<-c("VEL", 'ZEP', 'TKO', 'TRT', 'LAC', 'BAL', 'INE','SCH','WIL','SUB', 'BAB','HRA','SPI','DRA','TIS','CAR','ING','KAS','OBI','HCA','GUN','DRG')
  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
  s1<-subset(s,s$i.ID %in% id)
  ann1<-subset(ann,ann$AL %in% id)
  if (nrow(s1)>110)
  {p=0.75
  } else {p=1}
#  df<-as.matrix(s1[,4:25],rownames = paste(s1$i.ann,s1$i.aas,sep="   "))
  df<-as.matrix(s1[,4:25],rownames = paste(s1$i.aas,sep="   "))
  pdf(paste("heatmap_",id,".pdf",sep=""),height = 11,width = 9.3)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,4,6,7,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,11,13,15,16,18,19,20,21),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","black","black","black","black","black","black","black","black","black","grey","grey","black","black"),labCol=pops,colCol= c("#ED7C7F","#ED7C7F","#960C13","#960C13","#065570","#065570","#068FBD","#E0912F","#E0912F","#ED7C7F","#ED7C7F","#960C13","#960C13","#065570","#065570","#068FBD","#E0912F","#E0912F","#9acd32","#90EE90","#9acd32","#90EE90"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - ")) ###figuresGenes/
#  dev.off() #comment if uncommenting next part!!!!!!!!
  } else {} #comment if uncommenting next part!!!!!!!!
  
#  ###Arenosa
#  #repolarise
#  for (i in  1:nrow(s6)){ # i=1
#    if (s6$tot[i]>0.5)
#    {s6[i,4:21]<-1-s6[i,4:21]
#    } else {}}
#    if (nrow(subset(s6,s6$i.ID %in% id))>1)
#    {pops<-c("VEL", 'ZEP', 'TKO', 'TRT', 'LAC', 'BAL', 'INE','SCH','WIL','SUB', 'BAB','HRA','SPI','DRA','TIS','CAR','ING','KAS')
#    my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
#    s1<-subset(s6,s6$i.ID %in% id)
#    if (nrow(s1)>110)
#    {p=0.75
#    } else {p=1}
#    #  df<-as.matrix(s1[,4:25],rownames = paste(s1$i.ann,s1$i.aas,sep="   "))
#    df<-as.matrix(s1[,4:21],rownames = paste(s1$i.aas,sep="   "))
#    heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,4,6,7,9#,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,11,13,15,16,18),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("grey"#,"grey","grey","grey","grey","grey","grey","grey","grey","black","black","black","black","black","black","black","black","black"#),labCol=pops,colCol= c("grey","grey","grey","grey","grey","grey","grey","grey","grey","black","black","black","black","black","black"#,"black","black","black"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c#(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - "))
#    } else {}
#    
  ###Halleri
  #repolarise
  for (i in  1:nrow(s5)){ # i=1
    if (s5$tot[i]>0.5)
    {s5[i,4:7]<-1-s5[i,4:7]
    } else {}}
    if (nrow(subset(s5,s5$ID %in% id))>1)
    {pops<-c('OBI','GUN','HCA','DRG')
    my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
    s1<-subset(s5,s5$ID %in% id)
    if (nrow(s1)>110)
    {p=0.75
    } else {p=1}
    #  df<-as.matrix(s1[,4:25],rownames = paste(s1$i.ann,s1$i.aas,sep="   "))
    df<-as.matrix(s1[,4:7],rownames = paste(s1$aas,sep="   "))
    heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,4),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = c("grey","black","grey","black"),labCol=pops,colCol= c("#9acd32","#9acd32","#90EE90","#90EE90"),offsetRow = c(0.05),offsetCol= c(0.05),cexRow = c(p),cexCol = c(2),margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,ann1$AT,ann1$annShort,sep=" - "))
    dev.off()
    } else {}
    
}
}
