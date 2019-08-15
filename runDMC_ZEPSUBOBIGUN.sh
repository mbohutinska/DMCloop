#!/bin/bash -e
#PBS -N runDMC
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=20gb
#PBS -m abe
#PBS -j oe

module add R-3.4.3-gcc
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 

pops='ZEPSUBOBIGUN'

DATADIR="/storage/plzen1/home/holcovam/DMC/`echo $pops`"
cp $DATADIR/../runDMC.R $SCRATCHDIR || exit 1
cp $DATADIR/../neutralData/`echo $pops`.txt $SCRATCHDIR || exit 1
cp $DATADIR/../scripts/genSelMatrices_individualModes.R $SCRATCHDIR || exit 1
cp $DATADIR/../scripts/calcCompositeLike.R $SCRATCHDIR || exit 1
cp $DATADIR/../scripts/getMCLE.R $SCRATCHDIR || exit 1
cp $DATADIR/../scripts/calcNeutralF.R $SCRATCHDIR || exit 1
cp -r $DATADIR/* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

export R_LIBS="/auto/plzen1/home/holcovam/programs/Rpackages/"
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "16,16,16,16" --times_sv "1e6,2e6,4e6" --times_sv_source "1e5,1e6" --times_stagSweeps "0,50,500,5e3,1e4"
## times_sv minimum time is split of population pairs
## times_sv.source maximum time is split of population pairs
##times_stagSweeps maximum time is split of sisters

rm *.txt
rm *.R
cp -r $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
