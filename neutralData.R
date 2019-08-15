### 2. Generate neutral allele frequencies
library(data.table)
library(dplyr)
setwd("/home/aa/alpine/arenosaGenome/neutralomePM/")
#1.extract genes from file
ar3<-fread(paste('ALL.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
ha3<-fread(paste('ALH.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
colnames(ar3) <- c("pop","ploidy","scaff","start","AN","DP",seq(7,ncol(ar3),1))
colnames(ha3) <- c("pop","ploidy","scaff","start","AN","DP",seq(7,ncol(ha3),1))
#2.combine together
ar3$end<-ar3$start
ha3$end<-ha3$start
setkey(ha3, scaff, start,end)
a<-foverlaps(x = ar3, y = ha3, type="any")
a<-subset(a,!is.na(a$pop))
#3. calculate AF for each lineage
i<-which(colnames(a) %in% "7")
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
i<-which(colnames(a) %in% "i.7")
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
s1<-dplyr::select(a,VEL, ZEP, TKO, TRT, LAC, BAL, INE, SCH, WIL, SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG,start)
#export
p1<-"LAC"
p2<-"TIS"
p3<-"HCA"
p4<-"DRG"
sdmc<-dplyr::select(s1,start,p1,p2,p3,p4)
ss<-dplyr::select(s1,p1,p2,p3,p4)
sdmc$tot<-apply(ss,1,sum)
sdmc<-subset(sdmc,!sdmc$tot==0)
sdmc<-subset(sdmc,!sdmc$tot==4)
setwd("/home/aa/alpine/dmc/neutralData/")
write.table(sdmc[,2:5],paste(p1,p2,p3,p4,".txt",sep=""),col.names = F,row.names = F,quote = F)
