#!/bin/sh
#$ -l mem=5G,time=1:00:00
cd odtr
Rscript=/nfs/apps/R/4.0.3/bin/Rscript
export R_LIBS_USER=/ifs/home/msph/epi/ntw2117/R_4.0
${Rscript} _research/sim/sim2.R

exit 0
