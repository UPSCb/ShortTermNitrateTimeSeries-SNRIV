#!/bin/bash
#SBATCH -p main
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -A u2022003

R=/mnt/picea/storage/singularity/R-4.4.0.sif

Rlib=/mnt/picea/Modules/apps/compilers/R/4.4.0/lib/R/library:/usr/local/lib/R/library

myscript=$(realpath src/R/Network_Clustering_with_Leiden.R)

in=$(realpath data/seidr/results/backbone/backbone-1-percent.tsv)
# in=$(realpath data/seidr/results/backbone/backbone-9-percent.tsv)
# in=$(realpath data/seidr/results/backbone/backbone-5-percent.tsv)

out=$(realpath data/seidr/clustering)
[[ ! -d $out ]] && mkdir -p $out

export R_LIBS=/mnt/picea/Modules/apps/compilers/R/4.4.0/lib/R/library

# apptainer exec -B /mnt:/mnt -B $Rlib $R Rscript --vanilla $myscript $in $out/backbone-9-percent.tsv
# apptainer exec -B /mnt:/mnt -B $Rlib $R Rscript --vanilla $myscript $in $out/backbone-5-percent.tsv
apptainer exec -B /mnt:/mnt -B $Rlib $R Rscript --vanilla $myscript $in $out/backbone-1-percent.tsv
