#!/bin/sh
#$ -l mem=3G,time=:15:
cd odtr
Rscript=/nfs/apps/R/4.0.3/bin/Rscript
export R_LIBS_USER=/ifs/home/msph/epi/ntw2117/R_4.0
${Rscript} scripts/simulation/sim.R
