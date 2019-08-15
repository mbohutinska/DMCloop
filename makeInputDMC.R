#library(ggplot2,lib.loc="../programs/Rpackages/"
library(optparse)
library(data.table)
library(dplyr)
option_list = list(make_option(c("-g","--genes"), action="store", default=NA, type='character',help="character vector of gene names, separated by comma without space behind"),
                   make_option(c("-a","--p1"), action="store", type='character',default = NA,help = "population name"),
                   make_option(c("-b","--p2"), action="store", type='character',default = NA,help = "population name"),
                   make_option(c("-c","--p3"), action="store", type='character',default = NA,help = "population name"),
                   make_option(c("-d","--p4"), action="store", type='character',default = NA,help = "population name"))
opt = parse_args(OptionParser(option_list = option_list))
genes <-strsplit(opt$genes,split = ', ')[[1]]
p1<-opt$p1
p2<-opt$p2
p3<-opt$p3
p4<-opt$p4
dir.create(paste(p1,p2,p3,p4,sep="")) 
setwd(paste(p1,p2,p3,p4,sep=""))

ar<-fread(paste('../ALLarenosa.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
ha<-fread(paste('../ALLhalleri.table.recode.txt',sep=''),h=F,na.strings = "-9",nThread = 3)
colnames(ar) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ar),1))
colnames(ha) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",seq(10,ncol(ha),1))

#1.extract genes from file
for (id in genes){ #   id = "AL1G48930"
  dir.create(paste(id,sep=""))
  setwd(paste(id,sep=""))
  #start to loop over genes
  ar3<-subset(ar,ar$ID %in% id)
  mid<-min(ar3$start)+round(abs(ar3$start[1]-ar3$start[nrow(ar3)])/2,0)
  left<-mid-25000 # should be changed for each new species; there is formula saying how big window it should be, for A.lyrata recomb. rate and A. arenosas Ne is 50 kb fine
  right<-mid+25000
  ar3<-subset(ar,scaff %in% ar3$scaff[1] & start>=left & start<=right)
  ha3<-subset(ha,scaff %in% ar3$scaff[1] & start>=left & start<=right) 
  #2.combine together
  ar3$end<-ar3$start
  ha3$end<-ha3$start
  setkey(ha3, scaff, start,end)
  a<-foverlaps(x = ar3, y = ha3, type="any")
  setkey(ar3, scaff, start,end)
  b<-foverlaps(x = ha3, y = ar3, type="any")
  bb<-cbind(b[,1],b[,154:195],b[,2:153])
  colnames(bb)<-colnames(a)
  cc<-rbind(a,bb)
  for (i in  1:nrow(cc)){ # i=1
    if (is.na(cc$start[i]))
    {cc[i,4]<-cc[i,46]
    } else {}}
  for (i in  1:nrow(cc)){ # i=1
    if (is.na(cc$ID[i]))
    {cc[i,7]<-cc[i,49]
    } else {}}
  for (i in  1:nrow(cc)){ # i=1
    if (is.na(cc$ann[i]))
    {cc[i,8]<-cc[i,50]
    } else {}}
  for (i in  1:nrow(cc)){ # i=1
    if (is.na(cc$aas[i]))
    {cc[i,9]<-cc[i,51]
    } else {}}
  cc<-cc[ order(cc[,1], cc[,4]), ]
  cc<-cc[!duplicated(cc[,c('scaff','start')]),]
  a<-cc
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
  s1<-dplyr::select(a,ID,ann,aas,VEL, ZEP, TKO, TRT, LAC, BAL, INE,SCH,WIL,SUB, BAB,HRA,SPI,DRA,TIS,CAR,ING,KAS,OBI,HCA,GUN,DRG,start)
  for (i in  1:nrow(s1)){ # i=1
    if (is.na(s1[i,4]) & is.na(s1[i,5]) & is.na(s1[i,6]) & is.na(s1[i,7]) & is.na(s1[i,8]) & is.na(s1[i,9]) & is.na(s1[i,10]) & is.na(s1[i,11]) & is.na(s1[i,12]) & is.na(s1[i,13]))
    {s1[i,4:21]<-0
    } else {}}
  for (i in  1:nrow(s1)){ # i=1
    if (is.na(s1[i,22]) & is.na(s1[i,23]) & is.na(s1[i,24]) & is.na(s1[i,25]))
    {s1[i,22:25]<-0
    } else {}}
  #export
  sdmc<-dplyr::select(s1,start,ID,ann,aas,p1,p2,p3,p4)
  ss<-dplyr::select(s1,p1,p2,p3,p4)
  sdmc$tot<-apply(ss,1,sum)
  sdmc<-subset(sdmc,!sdmc$tot==0)
  sdmc<-subset(sdmc,!sdmc$tot==4)
  write.table(sdmc[,1:8],paste("dmcInput.txt",sep=""),col.names = F,row.names = F,quote = F)
  setwd("../")
}

