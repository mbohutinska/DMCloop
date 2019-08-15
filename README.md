# DMCloop - A set of scripts to run DMC analysis in a loop over multiple parallel population pairs and multiple loci
## Original method by Kristin Lee https://github.com/kristinmlee/dmc 
## Modified by Sivan Yair
1. run makeinputDMC - makes folder per lineage pair, all outliers inside (cca 3 min per gene)
cp makeInputDMC.sh makeInputDMC_INECARHCADRG.sh
cp makeInputDMC.sh makeInputDMC_INECAROBIGUN.sh
cp makeInputDMC.sh makeInputDMC_INECARTKOHRA.sh
cp makeInputDMC.sh makeInputDMC_INECARWILKAS.sh
cp makeInputDMC.sh makeInputDMC_INECARZEPSUB.sh
cp makeInputDMC.sh makeInputDMC_LACTISHCADRG.sh
cp makeInputDMC.sh makeInputDMC_LACTISINECAR.sh
cp makeInputDMC.sh makeInputDMC_LACTISOBIGUN.sh
cp makeInputDMC.sh makeInputDMC_LACTISTKOHRA.sh
cp makeInputDMC.sh makeInputDMC_LACTISWILKAS.sh
cp makeInputDMC.sh makeInputDMC_LACTISZEPSUB.sh
cp makeInputDMC.sh makeInputDMC_OBIGUNHCADRG.sh
cp makeInputDMC.sh makeInputDMC_TKOHRAHCADRG.sh
cp makeInputDMC.sh makeInputDMC_TKOHRAOBIGUN.sh
cp makeInputDMC.sh makeInputDMC_TKOHRAWILKAS.sh
cp makeInputDMC.sh makeInputDMC_TKOHRAZEPSUB.sh
cp makeInputDMC.sh makeInputDMC_WILKASHCADRG.sh
cp makeInputDMC.sh makeInputDMC_WILKASOBIGUN.sh
cp makeInputDMC.sh makeInputDMC_WILKASZEPSUB.sh
cp makeInputDMC.sh makeInputDMC_ZEPSUBHCADRG.sh
cp makeInputDMC.sh makeInputDMC_ZEPSUBOBIGUN.sh

##manually add pop names and outliers
qsub makeInputDMC_INECARHCADRG.sh
qsub makeInputDMC_INECAROBIGUN.sh
qsub makeInputDMC_INECARTKOHRA.sh
qsub makeInputDMC_INECARWILKAS.sh
qsub makeInputDMC_INECARZEPSUB.sh
qsub makeInputDMC_LACTISHCADRG.sh
qsub makeInputDMC_LACTISINECAR.sh
qsub makeInputDMC_LACTISOBIGUN.sh
qsub makeInputDMC_LACTISTKOHRA.sh
qsub makeInputDMC_LACTISWILKAS.sh
qsub makeInputDMC_LACTISZEPSUB.sh
qsub makeInputDMC_OBIGUNHCADRG.sh
qsub makeInputDMC_TKOHRAHCADRG.sh
qsub makeInputDMC_TKOHRAOBIGUN.sh
qsub makeInputDMC_TKOHRAWILKAS.sh
qsub makeInputDMC_TKOHRAZEPSUB.sh
#qsub makeInputDMC_WILKASHCADRG.sh
qsub makeInputDMC_WILKASOBIGUN.sh
qsub makeInputDMC_WILKASZEPSUB.sh
qsub makeInputDMC_ZEPSUBHCADRG.sh
qsub makeInputDMC_ZEPSUBOBIGUN.sh

2. make neutral matrices locally, upload

3. run DMC: manually add following
#between species:
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "16,16,16,16" --times_sv "1e5,1e6,2e6,4e6" --times_sv_source "0,1e3,1e4,1e5,1e6" --times_stagSweeps "0,50,500,5e3,1e4"
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "32,32,16,16" --times_sv "1e5,1e6,2e6,4e6" --times_sv_source "0,1e3,1e4,1e5,1e6" --times_stagSweeps "0,50,500,5e3,1e4"
#within species:
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "16,16,16,16" --times_sv "5e3,1e4,5e4,1e5" --times_sv_source "0,1e3,5e3,1e4" --times_stagSweeps "0,50,500,1e3,1e4"
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "32,32,16,16" --times_sv "5e3,1e4,5e4,1e5" --times_sv_source "0,1e3,5e3,1e4" --times_stagSweeps "0,50,500,1e3,1e4"
Rscript --vanilla runDMC.R --pops `echo $pops` --sampleSizes "32,32,32,32" --times_sv "5e3,1e4,5e4,1e5" --times_sv_source "0,1e3,5e3,1e4" --times_stagSweeps "0,50,500,1e3,1e4"
Rscript runDMC.R --pops TKOHRAHCADRG --sampleSizes "32,32,16,16" --times_sv "1e5" --times_sv_source "0" --times_stagSweeps "0"
## times_sv minimum time is split of population pairs
## times_sv.source maximum time is split of population pairs
## times_stagSweeps maximum time is split of sisters

find ./ -name "*maxParam.txt" -type f -exec rm {} \;

cp runDMC.sh runDMC_INECARHCADRG.sh
cp runDMC.sh runDMC_INECAROBIGUN.sh
cp runDMC.sh runDMC_INECARTKOHRA.sh
cp runDMC.sh runDMC_INECARWILKAS.sh
cp runDMC.sh runDMC_INECARZEPSUB.sh
cp runDMC.sh runDMC_LACTISHCADRG.sh
cp runDMC.sh runDMC_LACTISINECAR.sh
cp runDMC.sh runDMC_LACTISOBIGUN.sh
cp runDMC.sh runDMC_LACTISTKOHRA.sh
cp runDMC.sh runDMC_LACTISWILKAS.sh
cp runDMC.sh runDMC_LACTISZEPSUB.sh
cp runDMC.sh runDMC_OBIGUNHCADRG.sh
cp runDMC.sh runDMC_TKOHRAHCADRG.sh
cp runDMC.sh runDMC_TKOHRAOBIGUN.sh
cp runDMC.sh runDMC_TKOHRAWILKAS.sh
cp runDMC.sh runDMC_TKOHRAZEPSUB.sh
cp runDMC.sh runDMC_WILKASHCADRG.sh
cp runDMC.sh runDMC_WILKASOBIGUN.sh
cp runDMC.sh runDMC_WILKASZEPSUB.sh
cp runDMC.sh runDMC_ZEPSUBHCADRG.sh
cp runDMC.sh runDMC_ZEPSUBOBIGUN.sh

find ./ -name "*maxParam.txt" -type f -exec rm {} \;
qsub runDMC_INECARHCADRG.sh
qsub runDMC_INECAROBIGUN.sh
qsub runDMC_INECARTKOHRA.sh
qsub runDMC_INECARWILKAS.sh
qsub runDMC_INECARZEPSUB.sh
qsub runDMC_LACTISHCADRG.sh
qsub runDMC_LACTISINECAR.sh
qsub runDMC_LACTISOBIGUN.sh
qsub runDMC_LACTISTKOHRA.sh
qsub runDMC_LACTISWILKAS.sh
qsub runDMC_LACTISZEPSUB.sh
qsub runDMC_OBIGUNHCADRG.sh
qsub runDMC_TKOHRAHCADRG.sh
qsub runDMC_TKOHRAOBIGUN.sh
qsub runDMC_TKOHRAWILKAS.sh
qsub runDMC_TKOHRAZEPSUB.sh
#qsub runDMC_WILKASHCADRG.sh
qsub runDMC_WILKASOBIGUN.sh
qsub runDMC_WILKASZEPSUB.sh
qsub runDMC_ZEPSUBHCADRG.sh
qsub runDMC_ZEPSUBOBIGUN.sh



4. Download composite likelihoods
rsync -avm --include='compLikel.pdf' -f 'hide,! */' . analysis_5
rsync -avm --include='maxParam.txt' -f 'hide,! */' . analysis_5
rsync -avm --include='dmcInput.txt' -f 'hide,! */' . analysis_5


rsync -avm --include='*_neutralF.png' -f 'hide,! */' . analysis_3

find ./ -name "*maxParam.txt" -type f -exec rm {} \;

scp -r holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/DMC/analysis_5 .
scp /home/aa/alpine/arenosaGenome/selScans/ann/ALLarenosa.table.recode.txt holcovam@nympha.metacentrum.cz:/storage/plzen1/home/holcovam/DMC/annotatedData
