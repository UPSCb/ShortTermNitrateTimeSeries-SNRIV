#!/bin/bash -l 

in=$(realpath data/seidr/results/backbone/backbone-9-percent.tsv)
out=$(realpath data/seidr/clustering/backbone-9-percent.tsv)

# singularity
singr=$(realpath singularity/R-4.4.0.sif)

[[ ! -d $out ]] && mkdir -p $out

sbatch -t 12:00:00 -p main -n 2 --mem=16GB -A u2022003 \
-o $out/cluster_bb9.log -e $out/cluster_bb9.err \
src/R/Network_Clustering_with_Leiden.R $singr $in $out
# apptainer exec -B /mnt:/mnt -B /mnt/picea/Modules/apps/compilers/R/4.3.1/lib/R/library:/usr/local/lib/R/library /mnt/picea/storage/singularity/R-4.3.1.sif Rscript --vanilla correl.R 
