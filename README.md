# DMCloop 
### A set of scripts to run DMC analysis in a loop over multiple parallel population pairs and multiple loci
#### Original method by Kristin Lee https://github.com/kristinmlee/dmc, modified and extended by Sivan Yair (ssyair@ucdavis.edu) 
#### Wrapper and post-processing scripts by M. Bohutinska

## 1. run makeinputDMC 
#makes folder per lineage pair, all outliers inside (cca 3 min per gene)
#manually add pop names and outliers
qsub makeInputDMC_INECARHCADRG.sh 

## 2. make neutral matrices locally, upload 
#use neutralData.R or DMC part of selScans.R 

## 3. runDMC 
#add different parameters for witnin and between species
#between species:
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "16,16,16,16" --times_sv "1e6,2e6,4e6" --times_sv_source "1e5,1e6" --times_stagSweeps "0,50,500,5e3,1e4"
#within species:
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "16,16,16,16" --times_sv "5e3,1e4,5e4,1e5" --times_sv_source "100,1e3,5e3,1e4" --times_stagSweeps "0,50,500,1e3,1e4"

#times_sv minimum time is split of population pairs
#times_sv.source maximum time is split of population pairs
#times_stagSweeps maximum time is split of sisters

find ./ -name "*maxParam.txt" -type f -exec rm {} \;
qsub runDMC_INECARHCADRG.sh

## 4. Download composite likelihoods
rsync -avm --include='compLikel.pdf' -f 'hide,! */' . analysis_5
rsync -avm --include='maxParam.txt' -f 'hide,! */' . analysis_5
rsync -avm --include='dmcInput.txt' -f 'hide,! */' . analysis_5
rsync -avm --include='*_neutralF.png' -f 'hide,! */' . analysis_5
scp -r holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/DMC/analysis_5 .

## 5. Process the results
#make AFD dotplots
dotplot.R
#make AF heatmap
heatmapParallel.R
#check the coancestry for the ML model, needs its own genSelMatrices_individualModes.R function
investigate_coancestry_maxLikelihood.R


