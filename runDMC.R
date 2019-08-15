library("MASS")
library(optparse)
option_list = list(make_option(c("--pops"), action="store", default=NA, type='character',help="character vector of gene names, separated by comma without space behind"),
                   make_option(c("--sampleSizes"), action="store", type='character',default = NA,help = ""),
                   make_option(c("--times_sv"), action="store", type='character',default = NA,help = "minimum time is split of population pairs"),
                   make_option(c("--times_sv_source"), action="store", type='character',default = NA,help = "maximum time is split of population pairs"),
                   make_option(c("--times_stagSweeps"), action="store", type='character',default = NA,help = "maximum time is split of sisters"))
opt = parse_args(OptionParser(option_list = option_list))
pops<-opt$pops
sampleSizes<-as.numeric(strsplit(opt$sampleSizes,split = ',')[[1]])
times_sv = as.numeric(strsplit(opt$times_sv,split = ',')[[1]])
times_sv.source = as.numeric(strsplit(opt$times_sv_source,split = ',')[[1]])
times_stagSweeps = as.numeric(strsplit(opt$times_stagSweeps,split = ',')[[1]])

rec = 3.7e-08 
Ne = 800000
numPops = 4
numBins = 1000
sels = c(1e-4, 1e-3, 0.01, 0.05)
gs = c(1/(2*Ne), 10^-(4:3))
migs = c(10^-(seq(5, 1, by = -2)))

# ###
# times_sv = c(1e5,2e6) 
# times_sv.source = c(0,1e6)
# times_stagSweeps = 0
# pops = "TKOHRAHCADRG"
# sampleSizes<-as.numeric(strsplit("32,32,16,16",split = ',')[[1]])

## GENERATE NEUTRAL MATRIX - F
allFreqs<-as.matrix(t(read.table(paste(pops,".txt",sep=""))))
neutralF_filename = paste(pops,"_neutralF",sep="")
source("calcNeutralF.R")
F_estimate[F_estimate<0]=0

#plot(F_estimate)
#png(paste(pops,"_neutralF.png",sep=""))
#heatmap(F_estimate)
#dev.off()

for (gene in dir(path = ".",include.dirs = F,recursive = F,pattern = "AL.G*")) { # gene = "AL1G49940"
  print(gene) ########
setwd(paste(gene,sep=""))
a<-read.table("dmcInput.txt",h=F,sep = " ")
aa<-subset(a,a$V2 %in% gene)
left<-min(aa$V1)
right<-max(aa$V1)
se<-as.matrix(t(a[,5:8]))
positions = a$V1
selSite = seq(left-4000, right+4000, length.out = 8)
#1. Selection in alpine 
selPops = c(1, 3)
sources = selPops
## GENERATE MATRICES WITH selectION - F(s)
# 1 - independent mutations, 2 - migration, 3 - standing variation from source pop, 4 - sv, 5 - stagered sweeps
F_estimate = readRDS(paste("../",pops,"_neutralF.RDS",sep=""))
source("../genSelMatrices_individualModes.R")
# model 1
FOmegas_ind = lapply(sels, function(sel) {
  calcFOmegas_indSweeps(sel)
})
saveRDS(FOmegas_ind, "FOmegas_ind.RDS")
# model 2
FOmegas_mig = lapply(sels ,function(sel) {
  lapply(migs, function(mig) {
    lapply(sources, function(my.source) {
      calcFOmegas_mig(sel, mig, my.source)
    })
  })
})
saveRDS(FOmegas_mig, "FOmegas_mig.RDS")
# model 3
FOmegas_sv.source = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_sv.source, function(time) {
      lapply(sources, function(my.source) {
        calcFOmegas_stdVar.source(sel, g, time, my.source)
      })
    })
  })
})
saveRDS(FOmegas_sv.source, "FOmegas_sv_source.RDS")
# model 4
FOmegas_sv = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_sv, function(time) {
      calcFOmegas_stdVar(sel, g, time)
    })
  })
})
saveRDS(FOmegas_sv, "FOmegas_sv.RDS")
# model 5
#FOmegas_mig.stagSweeps = lapply(sels, function(sel) {
#  lapply(gs, function(g) {
#    lapply(times_stagSweeps, function(time_stagSweeps) { ### time
#      lapply(sources, function(my.source) {
#        calcFOmegas_mig.stagSweeps(sel, g, time_stagSweeps, my.source) ### time
#      })
#    })
#  })
#})
#saveRDS(FOmegas_mig.stagSweeps, "FOmegas_mig_stagSweeps.RDS")

## GENERATE INVERSES AND DETERMINANTS FOR F(s) MATRICES

## Neutral model
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)
det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(det_FOmegas_neutral, "det_FOmegas_neutral.RDS")
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(inv_FOmegas_neutral, "inv_FOmegas_neutral.RDS")
## Model 1
FOmegas_ind = readRDS("FOmegas_ind.RDS")
det_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
  lapply(sel, function(dist) {
    det(dist)
  })
})
saveRDS(det_FOmegas_ind, "det_FOmegas_ind.RDS")
inv_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
  lapply(sel, function(dist) {
    ginv(dist)
  })
})
saveRDS(inv_FOmegas_ind, "inv_FOmegas_ind.RDS")
## Model 2
FOmegas_mig = readRDS("FOmegas_mig.RDS")
det_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(source) {
      lapply(source, function(dist) {
        det(dist)
      })
    })
  })
})
saveRDS(det_FOmegas_mig, "det_FOmegas_mig.RDS")
inv_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(source) {
      lapply(source, function(dist) {
        ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_mig, "inv_FOmegas_mig.RDS")
## Model 3
FOmegas_sv.source = readRDS("FOmegas_sv_source.RDS")
det_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          det(dist)
        })
      })
    })
  })
})

# this returns list of length times_sv
# each element is a vector whose length is all parameter combinations
# that a given time parameter is included in 
# if all determinants for a given time are <0 (whole vector = T),
# then the given time parameter is the problem
det_byTime = lapply(1:length(times_sv.source), function(time){
  unlist(lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      lapply(1:length(sources),function(src){
        any(unlist(det_FOmegas_sv.source[[s]][[g]][[time]][[src]]) < 0)
      })
    })
  }))
})
# vector of length times_sv, returns TRUE for each 
# time that is the cause of negative determinants
to_keep = sapply(det_byTime, function(i){
  !any(i)
})
# to_keep = c()
# for(i in det_byTime){
#   to_keep = c(to_keep, !any(i))
# }

if(any(to_keep==F)){ # if we want to remove a time, do the removal, otherwise move on
  # REMOVE BAD TIMES
  FOmegas_sv.source = lapply(1:length(sels), function(sel) {
    lapply(1:length(gs), function(g) {
      lapply((1:length(times_sv.source))[to_keep],function(time){
        FOmegas_sv[[sel]][[g]][[time]]
      })
    })
  })
  
  print(paste('standing variant source model -- excluded times: ',paste(times_sv.source[!to_keep],collapse=',')))
  times_sv.source = times_sv.source[to_keep]
  
  det_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
    lapply(sel, function(g) {
      lapply(g, function(time) {
        lapply(time, function(my.source) {
          lapply(my.source, function(dist) {
            det(dist)
          })
        })
      })
    })
  })
}
saveRDS(det_FOmegas_sv.source, "det_FOmegas_sv_source.RDS")


inv_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          ginv(dist)
        })
      })
    })
  })
})
saveRDS(inv_FOmegas_sv.source, "inv_FOmegas_sv_source.RDS")
## Model 4
FOmegas_sv = readRDS("FOmegas_sv.RDS")
det_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        det(dist)
      })
    })
  })
})


# this returns list of length times_sv
# each element is a vector whose length is all parameter combinations
# that a given time parameter is included in 
# if all determinants for a given time are <0 (whole vector = T),
# then the given time parameter is the problem
det_byTime = lapply(1:length(times_sv), function(time){
  unlist(lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      any(unlist(det_FOmegas_sv[[s]][[g]][[time]]) < 0)
    })
  }))
})
# vector of length times_sv, returns TRUE for each 
# time that is the cause of negative determinants
to_keep = sapply(det_byTime, function(i){
  !any(i)
})
# to_keep = c()
# for(i in det_byTime){
#   to_keep = c(to_keep, !any(i))
# }


if(any(to_keep==F)){
  # REMOVE BAD TIMES
  FOmegas_sv = lapply(1:length(sels), function(sel) {
    lapply(1:length(gs), function(g) {
      FOmegas_sv[[sel]][[g]][to_keep]
    })
  })
  
  print(paste('standing variant model -- excluded times: ',paste(times_sv[!to_keep],collapse=',')))
  times_sv = times_sv[to_keep]
  
  
  # recalculate determinants
  det_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
    lapply(sel, function(g) {
      lapply(g, function(time) {
        lapply(time, function(dist) {
          det(dist)
        })
      })
    })
  })
}

saveRDS(det_FOmegas_sv, "det_FOmegas_sv.RDS")

inv_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_sv, "inv_FOmegas_sv.RDS")




### Model 5 
#FOmegas_mig.stagSweeps = readRDS("FOmegas_mig_stagSweeps.RDS")
#det_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel) {
#  lapply(sel, function(g) {
#    lapply(g, function(time) {
#      lapply(time, function(my.source) {
#        lapply(my.source, function(dist) {
#          det(dist)
#        })
#      })
#    })
#  })
#})
#saveRDS(det_FOmegas_mig.stagSweeps, "det_FOmegas_mig_stagSweeps.RDS")
#inv_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel) {
#  lapply(sel, function(g) {
#    lapply(g, function(time) {
#      lapply(time, function(my.source) {
#        lapply(my.source, function(dist) {
#          ginv(dist)
#        })
#      })
#    })
#  })
#})
#saveRDS(inv_FOmegas_mig.stagSweeps, "inv_FOmegas_mig_stagSweeps.RDS")

## CALCULATE COMPOSITE LIKELIHOODS
freqs_notRand = se
randFreqs = apply(freqs_notRand, 2, function(my.freqs) {
  if(runif(1) < 0.5) {
    my.freqs = 1 - my.freqs
  }
  my.freqs
})
saveRDS(randFreqs, "selectedRegionAlleleFreqsRand.RDS")
freqs = readRDS("selectedRegionAlleleFreqsRand.RDS")
source("../calcCompositeLike.R")
## Neutral model
det_FOmegas_neutral = readRDS("det_FOmegas_neutral.RDS")
inv_FOmegas_neutral = readRDS("inv_FOmegas_neutral.RDS")
compLikelihood_neutral = lapply(1 : length(selSite), function(j) { # j=1
  calcCompLikelihood_neutral(j, det_FOmegas_neutral, inv_FOmegas_neutral)
})
saveRDS(compLikelihood_neutral, "compLikelihood_neutral.RDS")
## Model 1
det_FOmegas_ind = readRDS("det_FOmegas_ind.RDS")
inv_FOmegas_ind = readRDS("inv_FOmegas_ind.RDS")
compLikelihood_ind = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) calcCompLikelihood_1par(j, det_FOmegas_ind,
                                                                 inv_FOmegas_ind, sel))
})
saveRDS(compLikelihood_ind, "compLikelihood_ind.RDS")
## Model 2
det_FOmegas_mig = readRDS("det_FOmegas_mig.RDS")
inv_FOmegas_mig = readRDS("inv_FOmegas_mig.RDS")
compLikelihood_mig = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(migs), function(mig) {
      lapply(1 : length(sources), function(my.source) {
        calcCompLikelihood_3par(j, det_FOmegas_mig, inv_FOmegas_mig, sel, mig,
                                my.source)
      })
    })
  })
})
saveRDS(compLikelihood_mig, "compLikelihood_mig.RDS")
## Model 3
det_FOmegas_sv.source = readRDS("det_FOmegas_sv_source.RDS")
inv_FOmegas_sv.source = readRDS("inv_FOmegas_sv_source.RDS")
compLikelihood_sv.source = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_sv.source), function(t) {
        lapply(1: length(sources), function(my.source) {
          calcCompLikelihood_4par(j, det_FOmegas_sv.source, inv_FOmegas_sv.source, sel, g, t,
                                  my.source)
        })
      })
    })
  })
})
saveRDS(compLikelihood_sv.source, "compLikelihood_sv_source.RDS")

## Model 4
det_FOmegas_sv = readRDS("det_FOmegas_sv.RDS")
inv_FOmegas_sv = readRDS("inv_FOmegas_sv.RDS")
compLikelihood_sv = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_sv), function(t) {
        calcCompLikelihood_3par(j, det_FOmegas_sv, inv_FOmegas_sv, sel, g, t)
      })
    })
  })
})
saveRDS(compLikelihood_sv, "compLikelihood_sv.RDS")

## Model 5
#det_FOmegas_mig.stagSweeps = readRDS("det_FOmegas_mig_stagSweeps.RDS")
#inv_FOmegas_mig.stagSweeps = readRDS("inv_FOmegas_mig_stagSweeps.RDS")
#compLikelihood_mig.stagSweeps = lapply(1 : length(selSite), function(j) {
#  lapply(1 : length(sels), function(sel) {
#    lapply(1 : length(gs), function(g) {
#      lapply(1 : length(times_stagSweeps), function(t) {
#        lapply(1: length(sources), function(my.source) {
#          calcCompLikelihood_4par(j, det_FOmegas_mig.stagSweeps, inv_FOmegas_mig.stagSweeps, sel, g, t,my.source)
#        })
#      })
#    })
#  })
#})
#saveRDS(compLikelihood_mig.stagSweeps, "compLikelihood_mig_stagSweeps.RDS")

###GET MCL ESTIMATES FOR PARAMETERS
source("../getMCLE.R")
## Model 1
one<-getMCLEind(compLikelihood_ind, selSite, sels)
## Model 2
two<-getMCLEmig(compLikelihood_mig, selSite, sels, migs, sources)
## Model 3
three<-getMCLEsv_source(compLikelihood_sv, selSite, sels, gs, times_sv.source, sources)
## Model 4
four<-getMCLEsv(compLikelihood_sv, selSite, sels, gs, times_sv)
## Model 5
#five<-getMCLEmig_stagSweeps(compLikelihood_mig.stagSweeps, selSite, sels, gs, times_stagSweeps, sources)

write.table(paste("alpine de-novo",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(one),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(paste("alpine migration",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(two),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(paste("alpine standing source",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(three),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(paste("alpine standing ancestral",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(four),"maxParam.txt",append=T,col.names = T,row.names =F)
#write.table(as.data.frame(five),"maxParam.txt",append=T,col.names = T,row.names =F)
## PLOT
#read in composite likelihood files and calculate max for all proposed selected sites
compLikelihood_neutral = readRDS("compLikelihood_neutral.RDS")
compLikelihood_neutral_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_neutral[[i]]))
})
compLikelihood_ind = readRDS("compLikelihood_ind.RDS")
compLikelihood_ind_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_ind[[i]]))
})
compLikelihood_mig = readRDS("compLikelihood_mig.RDS")
compLikelihood_mig_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_mig[[i]]))
})
compLikelihood_sv.source = readRDS("compLikelihood_sv_source.RDS")
compLikelihood_sv.source_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_sv.source[[i]]))
})
compLikelihood_sv = readRDS("compLikelihood_sv.RDS")
compLikelihood_sv_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_sv[[i]]))
})
#compLikelihood_mig.stagSweeps = readRDS("compLikelihood_mig_stagSweeps.RDS")
#compLikelihood_mig.stagSweeps_site = sapply(1 : length(selSite), function(i) {  
#  max(unlist(compLikelihood_mig.stagSweeps[[i]]))
#})

plot_range = range(c((compLikelihood_ind_site - compLikelihood_neutral_site),
                     (compLikelihood_mig_site - compLikelihood_neutral_site), 
                     (compLikelihood_sv.source_site - compLikelihood_neutral_site),
                     (compLikelihood_sv_site - compLikelihood_neutral_site)))

write.table(as.data.frame(compLikelihood_ind_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(compLikelihood_mig_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(compLikelihood_sv.source_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(compLikelihood_sv_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
#write.table(as.data.frame(compLikelihood_mig.stagSweeps_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)

pdf("compLikel.pdf",width = 9,height = 9,pointsize = 15)
plot(selSite, compLikelihood_ind_site - compLikelihood_neutral_site, type = "b",
     ylim = c(plot_range[1] - 50, plot_range[2] + 50),
     xlab = "Proposed position selected site",
     ylab = "Composite log-likelihood (model - neutral)",main=gene,sub = "Selection in alpine")
lines(selSite, compLikelihood_mig_site - compLikelihood_neutral_site, col = "red",
      type = "b")
lines(selSite, compLikelihood_sv.source_site - compLikelihood_neutral_site, col = "green3",
      type = "b")
lines(selSite, compLikelihood_sv_site - compLikelihood_neutral_site, col = "blue",
      type = "b")
#lines(selSite, compLikelihood_mig.stagSweeps_site - compLikelihood_neutral_site, col = "orange2",type = "b")
legend("topright", col = c("black", "red","green3", "blue"),
       lty = 1, c("De-novo mut.","Migration","SV source", "SV ancestral"),
       cex = 0.5)
segments(x0 = left,y0 = plot_range[2]+40,x1 = right,y1 = plot_range[2]+40,col = "grey50",lwd = 4)
abline(h = 0,col="grey90",lty=2)

#2. Selection in foothill 
#To account for sweep happening in lowland repeatedly
selPops = c(2, 4)
sources = selPops
source("../genSelMatrices_individualModes.R")
# model 1
FOmegas_ind = lapply(sels, function(sel) {
  calcFOmegas_indSweeps(sel)
})
saveRDS(FOmegas_ind, "FOmegas_ind.RDS")
# model 2
FOmegas_mig = lapply(sels ,function(sel) {
  lapply(migs, function(mig) {
    lapply(sources, function(my.source) {
      calcFOmegas_mig(sel, mig, my.source)
    })
  })
})
saveRDS(FOmegas_mig, "FOmegas_mig.RDS")
# model 3
FOmegas_sv.source = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_sv.source, function(time) {
      lapply(sources, function(my.source) {
        calcFOmegas_stdVar.source(sel, g, time, my.source)
      })
    })
  })
})
saveRDS(FOmegas_sv.source, "FOmegas_sv_source.RDS")
# model 4
FOmegas_sv = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times_sv, function(time) {
      calcFOmegas_stdVar(sel, g, time)
    })
  })
})
saveRDS(FOmegas_sv, "FOmegas_sv.RDS")
# model 5
#FOmegas_mig.stagSweeps = lapply(sels, function(sel) {
#  lapply(gs, function(g) {
#    lapply(times_stagSweeps, function(time_stagSweeps) { ### time
#      lapply(sources, function(my.source) {
#        calcFOmegas_mig.stagSweeps(sel, g, time_stagSweeps, my.source) ### time
#      })
#    })
#  })
#})
#saveRDS(FOmegas_mig.stagSweeps, "FOmegas_mig_stagSweeps.RDS")

## GENERATE INVERSES AND DETERMINANTS FOR F(s) MATRICES

## Neutral model
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)
det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(det_FOmegas_neutral, "det_FOmegas_neutral.RDS")
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(inv_FOmegas_neutral, "inv_FOmegas_neutral.RDS")
## Model 1
FOmegas_ind = readRDS("FOmegas_ind.RDS")
det_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
  lapply(sel, function(dist) {
    det(dist)
  })
})
saveRDS(det_FOmegas_ind, "det_FOmegas_ind.RDS")
inv_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
  lapply(sel, function(dist) {
    ginv(dist)
  })
})
saveRDS(inv_FOmegas_ind, "inv_FOmegas_ind.RDS")
## Model 2
FOmegas_mig = readRDS("FOmegas_mig.RDS")
det_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(source) {
      lapply(source, function(dist) {
        det(dist)
      })
    })
  })
})
saveRDS(det_FOmegas_mig, "det_FOmegas_mig.RDS")
inv_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(source) {
      lapply(source, function(dist) {
        ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_mig, "inv_FOmegas_mig.RDS")
## Model 3
FOmegas_sv.source = readRDS("FOmegas_sv_source.RDS")
det_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          det(dist)
        })
      })
    })
  })
})

# this returns list of length times_sv
# each element is a vector whose length is all parameter combinations
# that a given time parameter is included in 
# if all determinants for a given time are <0 (whole vector = T),
# then the given time parameter is the problem
det_byTime = lapply(1:length(times_sv.source), function(time){
  unlist(lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      lapply(1:length(sources),function(src){
        any(unlist(det_FOmegas_sv.source[[s]][[g]][[time]][[src]]) < 0)
      })
    })
  }))
})
# vector of length times_sv, returns TRUE for each 
# time that is the cause of negative determinants
to_keep = sapply(det_byTime, function(i){
  !any(i)
})
# to_keep = c()
# for(i in det_byTime){
#   to_keep = c(to_keep, !any(i))
# }

if(any(to_keep==F)){ # if we want to remove a time, do the removal, otherwise move on
  # REMOVE BAD TIMES
  FOmegas_sv.source = lapply(1:length(sels), function(sel) {
    lapply(1:length(gs), function(g) {
      lapply((1:length(times_sv.source))[to_keep],function(time){
        FOmegas_sv[[sel]][[g]][[time]]
      })
    })
  })
  
  print(paste('standing variant source model -- excluded times: ',paste(times_sv.source[!to_keep],collapse=',')))
  times_sv.source = times_sv.source[to_keep]
  
  det_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
    lapply(sel, function(g) {
      lapply(g, function(time) {
        lapply(time, function(my.source) {
          lapply(my.source, function(dist) {
            det(dist)
          })
        })
      })
    })
  })
}
saveRDS(det_FOmegas_sv.source, "det_FOmegas_sv_source.RDS")


inv_FOmegas_sv.source = lapply(FOmegas_sv.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(my.source) {
        lapply(my.source, function(dist) {
          ginv(dist)
        })
      })
    })
  })
})
saveRDS(inv_FOmegas_sv.source, "inv_FOmegas_sv_source.RDS")
## Model 4
FOmegas_sv = readRDS("FOmegas_sv.RDS")
det_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        det(dist)
      })
    })
  })
})


# this returns list of length times_sv
# each element is a vector whose length is all parameter combinations
# that a given time parameter is included in 
# if all determinants for a given time are <0 (whole vector = T),
# then the given time parameter is the problem
det_byTime = lapply(1:length(times_sv), function(time){
  unlist(lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      any(unlist(det_FOmegas_sv[[s]][[g]][[time]]) < 0)
    })
  }))
})
# vector of length times_sv, returns TRUE for each 
# time that is the cause of negative determinants
to_keep = sapply(det_byTime, function(i){
  !any(i)
})
# to_keep = c()
# for(i in det_byTime){
#   to_keep = c(to_keep, !any(i))
# }


if(any(to_keep==F)){
  # REMOVE BAD TIMES
  FOmegas_sv = lapply(1:length(sels), function(sel) {
    lapply(1:length(gs), function(g) {
      FOmegas_sv[[sel]][[g]][to_keep]
    })
  })
  
  print(paste('standing variant model -- excluded times: ',paste(times_sv[!to_keep],collapse=',')))
  times_sv = times_sv[to_keep]
  
  
  # recalculate determinants
  det_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
    lapply(sel, function(g) {
      lapply(g, function(time) {
        lapply(time, function(dist) {
          det(dist)
        })
      })
    })
  })
}

saveRDS(det_FOmegas_sv, "det_FOmegas_sv.RDS")

inv_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_sv, "inv_FOmegas_sv.RDS")




### Model 5 
#FOmegas_mig.stagSweeps = readRDS("FOmegas_mig_stagSweeps.RDS")
#det_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel) {
#  lapply(sel, function(g) {
#    lapply(g, function(time) {
#      lapply(time, function(my.source) {
#        lapply(my.source, function(dist) {
#          det(dist)
#        })
#      })
#    })
#  })
#})
#saveRDS(det_FOmegas_mig.stagSweeps, "det_FOmegas_mig_stagSweeps.RDS")
#inv_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel) {
#  lapply(sel, function(g) {
#    lapply(g, function(time) {
#      lapply(time, function(my.source) {
#        lapply(my.source, function(dist) {
#          ginv(dist)
#        })
#      })
#    })
#  })
#})
#saveRDS(inv_FOmegas_mig.stagSweeps, "inv_FOmegas_mig_stagSweeps.RDS")

## CALCULATE COMPOSITE LIKELIHOODS
freqs_notRand = se
randFreqs = apply(freqs_notRand, 2, function(my.freqs) {
  if(runif(1) < 0.5) {
    my.freqs = 1 - my.freqs
  }
  my.freqs
})
saveRDS(randFreqs, "selectedRegionAlleleFreqsRand.RDS")
freqs = readRDS("selectedRegionAlleleFreqsRand.RDS")
source("../calcCompositeLike.R")
## Neutral model
det_FOmegas_neutral = readRDS("det_FOmegas_neutral.RDS")
inv_FOmegas_neutral = readRDS("inv_FOmegas_neutral.RDS")
compLikelihood_neutral = lapply(1 : length(selSite), function(j) { # j=1
  calcCompLikelihood_neutral(j, det_FOmegas_neutral, inv_FOmegas_neutral)
})
saveRDS(compLikelihood_neutral, "compLikelihood_neutral.RDS")
## Model 1
det_FOmegas_ind = readRDS("det_FOmegas_ind.RDS")
inv_FOmegas_ind = readRDS("inv_FOmegas_ind.RDS")
compLikelihood_ind = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) calcCompLikelihood_1par(j, det_FOmegas_ind,
                                                                 inv_FOmegas_ind, sel))
})
saveRDS(compLikelihood_ind, "compLikelihood_ind.RDS")
## Model 2
det_FOmegas_mig = readRDS("det_FOmegas_mig.RDS")
inv_FOmegas_mig = readRDS("inv_FOmegas_mig.RDS")
compLikelihood_mig = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(migs), function(mig) {
      lapply(1 : length(sources), function(my.source) {
        calcCompLikelihood_3par(j, det_FOmegas_mig, inv_FOmegas_mig, sel, mig,
                                my.source)
      })
    })
  })
})
saveRDS(compLikelihood_mig, "compLikelihood_mig.RDS")
## Model 3
det_FOmegas_sv.source = readRDS("det_FOmegas_sv_source.RDS")
inv_FOmegas_sv.source = readRDS("inv_FOmegas_sv_source.RDS")
compLikelihood_sv.source = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_sv.source), function(t) {
        lapply(1: length(sources), function(my.source) {
          calcCompLikelihood_4par(j, det_FOmegas_sv.source, inv_FOmegas_sv.source, sel, g, t,
                                  my.source)
        })
      })
    })
  })
})
saveRDS(compLikelihood_sv.source, "compLikelihood_sv_source.RDS")

## Model 4
det_FOmegas_sv = readRDS("det_FOmegas_sv.RDS")
inv_FOmegas_sv = readRDS("inv_FOmegas_sv.RDS")
compLikelihood_sv = lapply(1 : length(selSite), function(j) {
  lapply(1 : length(sels), function(sel) {
    lapply(1 : length(gs), function(g) {
      lapply(1 : length(times_sv), function(t) {
        calcCompLikelihood_3par(j, det_FOmegas_sv, inv_FOmegas_sv, sel, g, t)
      })
    })
  })
})
saveRDS(compLikelihood_sv, "compLikelihood_sv.RDS")

## Model 5
#det_FOmegas_mig.stagSweeps = readRDS("det_FOmegas_mig_stagSweeps.RDS")
#inv_FOmegas_mig.stagSweeps = readRDS("inv_FOmegas_mig_stagSweeps.RDS")
#compLikelihood_mig.stagSweeps = lapply(1 : length(selSite), function(j) {
#  lapply(1 : length(sels), function(sel) {
#    lapply(1 : length(gs), function(g) {
#      lapply(1 : length(times_stagSweeps), function(t) {
#        lapply(1: length(sources), function(my.source) {
#          calcCompLikelihood_4par(j, det_FOmegas_mig.stagSweeps, inv_FOmegas_mig.stagSweeps, sel, g, t,my.source)
#        })
#      })
#    })
#  })
#})
#saveRDS(compLikelihood_mig.stagSweeps, "compLikelihood_mig_stagSweeps.RDS")

###GET MCL ESTIMATES FOR PARAMETERS
source("../getMCLE.R")
## Model 1
one<-getMCLEind(compLikelihood_ind, selSite, sels)
## Model 2
two<-getMCLEmig(compLikelihood_mig, selSite, sels, migs, sources)
## Model 3
three<-getMCLEsv_source(compLikelihood_sv, selSite, sels, gs, times_sv.source, sources)
## Model 4
four<-getMCLEsv(compLikelihood_sv, selSite, sels, gs, times_sv)
## Model 5
#five<-getMCLEmig_stagSweeps(compLikelihood_mig.stagSweeps, selSite, sels, gs, times_stagSweeps, sources)
write.table(paste("foothill de-novo",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(one),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(paste("foothill migration",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(two),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(paste("foothill standing source",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(three),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(paste("foothill standing ancestral",sep=""),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(four),"maxParam.txt",append=T,col.names = T,row.names =F)
#write.table(as.data.frame(five),"maxParam.txt",append=T,col.names = T,row.names =F)
## PLOT
#read in composite likelihood files and calculate max for all proposed selected sites
compLikelihood_neutral = readRDS("compLikelihood_neutral.RDS")
compLikelihood_neutral_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_neutral[[i]]))
})
compLikelihood_ind = readRDS("compLikelihood_ind.RDS")
compLikelihood_ind_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_ind[[i]]))
})
compLikelihood_mig = readRDS("compLikelihood_mig.RDS")
compLikelihood_mig_site = sapply(1 : length(selSite), function(i) {
  max(unlist(compLikelihood_mig[[i]]))
})
compLikelihood_sv.source = readRDS("compLikelihood_sv_source.RDS")
compLikelihood_sv.source_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_sv.source[[i]]))
})
compLikelihood_sv = readRDS("compLikelihood_sv.RDS")
compLikelihood_sv_site = sapply(1 : length(selSite), function(i) {  
  max(unlist(compLikelihood_sv[[i]]))
})
#compLikelihood_mig.stagSweeps = readRDS("compLikelihood_mig_stagSweeps.RDS")
#compLikelihood_mig.stagSweeps_site = sapply(1 : length(selSite), function(i) {  
#  max(unlist(compLikelihood_mig.stagSweeps[[i]]))
#})

plot_range = range(c((compLikelihood_ind_site - compLikelihood_neutral_site),
                     (compLikelihood_mig_site - compLikelihood_neutral_site), 
                     (compLikelihood_sv.source_site - compLikelihood_neutral_site),
                     (compLikelihood_sv_site - compLikelihood_neutral_site)))

write.table(as.data.frame(compLikelihood_ind_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(compLikelihood_mig_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(compLikelihood_sv.source_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
write.table(as.data.frame(compLikelihood_sv_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)
#write.table(as.data.frame(compLikelihood_mig.stagSweeps_site - compLikelihood_neutral_site),"maxParam.txt",append=T,col.names = T,row.names =F)

#pdf("compLikel.pdf",width = 9,height = 9,pointsize = 15)
plot(selSite, compLikelihood_ind_site - compLikelihood_neutral_site, type = "b",
     ylim = c(plot_range[1] - 50, plot_range[2] + 50),
     xlab = "Proposed position selected site",
     ylab = "Composite log-likelihood (model - neutral)",main=gene,sub = "Selection in foothill")
lines(selSite, compLikelihood_mig_site - compLikelihood_neutral_site, col = "red",
      type = "b")
lines(selSite, compLikelihood_sv.source_site - compLikelihood_neutral_site, col = "green3",
      type = "b")
lines(selSite, compLikelihood_sv_site - compLikelihood_neutral_site, col = "blue",
      type = "b")
#lines(selSite, compLikelihood_mig.stagSweeps_site - compLikelihood_neutral_site, col = "orange2",type = "b")
legend("topright", col = c("black", "red","green3", "blue"),
       lty = 1, c("De-novo mut.","Migration","SV source", "SV ancestral"),
       cex = 0.5)
segments(x0 = left,y0 = plot_range[2]+40,x1 = right,y1 = plot_range[2]+40,col = "grey50",lwd = 4)
abline(h = 0,col="grey90",lty=2)
dev.off()


setwd("../")
}
