#!/bin/bash -e
#PBS -N makeInputDMC
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=20gb
#PBS -m abe
#PBS -j oe

module add R-3.4.3-gcc
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/DMC"
cp $DATADIR/annotatedData/ALL* $SCRATCHDIR || exit 1
cp $DATADIR/makeInputDMC.R $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

export R_LIBS="/auto/plzen1/home/holcovam/programs/Rpackages/"
Rscript --vanilla makeInputDMC.R --genes "AL1G18450, AL1G18470, AL1G22390, AL1G24590, AL1G30800, AL1G48930, AL2G16960, AL2G22180, AL2G22560, AL3G12310, AL3G53290, AL4G40380, AL4G42370, AL4G42990, AL4G47380, AL4G47600, AL4G47650, AL5G10130, AL5G22700, AL5G29990, AL6G21830, AL6G38530, AL6G42430, AL7G13680, AL7G19110, AL7G20350, AL7G21580, AL7G25530, AL7G30950, AL7G31090, AL7G34350, AL8G13260, AL8G13270, AL8G16970, AL8G16980, AL8G23730, AL8G24230, AL8G26390, AL8G36020, AL8G39200" --p1 "LAC" --p2 "TIS" --p3 "ZEP" --p4 "SUB"
rm ALL*
rm makeInputDMC.R
cp -r $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
