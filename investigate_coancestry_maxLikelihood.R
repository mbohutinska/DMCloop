setwd('/home/aa/alpine/dmc/Convergence_Arabidopsis/')
F_estimate = readRDS('ZEPSUBOBIGUN_neutralF.RDS')
#midDistances = read.table('midDist.txt')$V1

numPops=ncol(F_estimate)
selPops=c(1,3)
nonSelPops=c(2,4)
gene<-"AL2G16960"
pops_byLoc = c("ZEP","SUB","OBI","GUN")
pops_byEnv = c("Alpine Arenosa VT", "Foothill Arenosa VT", "Alpine Halleri HN", "Foothill Halleri HN")
sampleSizes = c(rep(16, 2),rep(16,2))
rec = 3.7e-08 
Ne = 800000
numPops = 4
selPops = c(1, 3)
numBins = 1000
sources = selPops
#sels = c(1e-4, 1e-3, 0.01, 0.05, 0.5)
#times_sv = c(1e4,1e5,1e6) ## minimum time is split of population pairs
#times_sv.source = c(0,1000,5000,1e4)  ## maximum time is split of population pairs
#times_stagSweeps = c(0,50, 500, 1000,5000) ## maximum time is split of sisters
#gs = c(1/(2*Ne), 10^-(4:1))
#migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)



a<-read.table("dmcInput.txt",h=F,sep = " ")
aa<-subset(a,a$V2 %in% gene)
left<-min(aa$V1)
right<-max(aa$V1)
se<-as.matrix(t(a[,5:8]))
positions = a$V1
selSite = seq(left-4000, right+4000, length.out = 6)

source('genSelMatrices_individualModes.R')

# all functions below generate a list of length midDistances, in which 
# each element is a matrix describing the probability of coalescing
# for all population pairs under a given model, provided max likelihood params, 
# and distance from the selected site

genetic_dist = midDistances*rec # for plotting

# IND SWEEPS (de novo)
F_ind = calcFOmegas_indSweeps(sel=0.01)
length(F_ind)
# plot decay in prob of coalescing as we move away from selected site, for all pairs of populations
pdf('ZEPSUBOBIGUN/IndSweeps.pdf')
# within populations first
for(i in 1:numPops){
  # get within pop prob of coalescing at all genetic distances
  F_ii = sapply(1:length(midDistances), function(dist) F_ind[[dist]][i,i] ) # dist = 1000
  plot(F_ii~genetic_dist,ylim=c(-0.3,1),main=paste('Within',pops_byEnv[i]),type='l',lwd=2)
  abline(h=F_estimate[i,i],lty=2,col='red')
}
# between populations
for(i in 1:numPops){
  for(j in 1:numPops){
    if(j > i){
      F_ij = sapply(1:length(midDistances), function(dist) F_ind[[dist]][i,j] )
      plot(F_ij~genetic_dist,ylim=c(-0.3,1),main=paste('Between',pops_byEnv[i],'&',pops_byEnv[j]),type='l',lwd=2)
      abline(h=F_estimate[i,j],lty=2,col='red')
    }
  }
}
dev.off()

# CONC SWEEPS
F_mig_conc = calcFOmegas_mig(sel=1e-04, mig=1e-05, my.source=1)
# plot decay in prob of coalescing as we move away from selected site, for all pairs of populations
pdf('Mig_ConcSweeps.pdf')
# within populations first
for(i in 1:numPops){
  # get within pop prob of coalescing at all genetic distances
  F_ii = sapply(1:length(midDistances), function(dist) F_mig_conc[[dist]][i,i] )
  plot(F_ii~genetic_dist,ylim=c(-0.3,1),main=paste('Within',pops_byEnv[i]),type='l',lwd=2)
  abline(h=F_estimate[i,i],lty=2,col='red')
}
# between populations
for(i in 1:numPops){
  for(j in 1:numPops){
    if(j > i){
      F_ij = sapply(1:length(midDistances), function(dist) F_mig_conc[[dist]][i,j] )
      plot(F_ij~genetic_dist,ylim=c(-0.3,1),main=paste('Between',pops_byEnv[i],'&',pops_byEnv[j]),type='l',lwd=2)
      abline(h=F_estimate[i,j],lty=2,col='red')
    }
  }
}
dev.off()

# SV SOURCE
F_sv_src = calcFOmegas_stdVar.source(sel=1e-03, g=6.25e-07, time=0, my.source=1)
# plot decay in prob of coalescing as we move away from selected site, for all pairs of populations
pdf('ZEPSUBOBIGUN/StdVar_source.pdf')
# within populations first
for(i in 1:numPops){
  # get within pop prob of coalescing at all genetic distances
  F_ii = sapply(1:length(midDistances), function(dist) F_sv_src[[dist]][i,i] )
  plot(F_ii~genetic_dist,ylim=c(-0.3,1),main=paste('Within',pops_byEnv[i]),type='l',lwd=2)
  abline(h=F_estimate[i,i],lty=2,col='red')
}
# between populations
for(i in 1:numPops){
  for(j in 1:numPops){
    if(j > i){
      F_ij = sapply(1:length(midDistances), function(dist) F_sv_src[[dist]][i,j] )
      plot(F_ij~genetic_dist,ylim=c(-0.3,1),main=paste('Between',pops_byEnv[i],'&',pops_byEnv[j]),type='l',lwd=2)
      abline(h=F_estimate[i,j],lty=2,col='red')
    }
  }
}
dev.off()


# SV
F_sv = calcFOmegas_stdVar(sel=0.001, g=0.001, time=1e+05)
# plot decay in prob of coalescing as we move away from selected site, for all pairs of populations
pdf('StdVar.pdf')
# within populations first
for(i in 1:numPops){
  # get within pop prob of coalescing at all genetic distances
  F_ii = sapply(1:length(midDistances), function(dist) F_sv[[dist]][i,i] )
  plot(F_ii~genetic_dist,ylim=c(-0.3,1),main=paste('Within',pops_byEnv[i]),type='l',lwd=2)
  abline(h=F_estimate[i,i],lty=2,col='red')
}
# between populations
for(i in 1:numPops){
  for(j in 1:numPops){
    if(j > i){
      F_ij = sapply(1:length(midDistances), function(dist) F_sv[[dist]][i,j] )
      plot(F_ij~genetic_dist,ylim=c(-0.3,1),main=paste('Between',pops_byEnv[i],'&',pops_byEnv[j]),type='l',lwd=2)
      abline(h=F_estimate[i,j],lty=2,col='red')
    }
  }
}
dev.off()

# STAG SWEEPS
F_mig_stag = calcFOmegas_mig.stagSweeps(sel=0.01, G=0.001, standing.time=10000,my.source=3)
# plot decay in prob of coalescing as we move away from selected site, for all pairs of populations
pdf('Mig_StagSweeps.pdf')
# within populations first
for(i in 1:numPops){
  # get within pop prob of coalescing at all genetic distances
  F_ii = sapply(1:length(midDistances), function(dist) F_mig_stag[[dist]][i,i] )
  plot(F_ii~genetic_dist,ylim=c(-0.3,1),main=paste('Within',pops_byEnv[i]),type='l',lwd=2)
  abline(h=F_estimate[i,i],lty=2,col='red')
}
# between populations
for(i in 1:numPops){
  for(j in 1:numPops){
    if(j > i){
      F_ij = sapply(1:length(midDistances), function(dist) F_mig_stag[[dist]][i,j] )
      plot(F_ij~genetic_dist,ylim=c(-0.3,1),main=paste('Between',pops_byEnv[i],'&',pops_byEnv[j]),type='l',lwd=2)
      abline(h=F_estimate[i,j],lty=2,col='red')
    }
  }
}
dev.off()



