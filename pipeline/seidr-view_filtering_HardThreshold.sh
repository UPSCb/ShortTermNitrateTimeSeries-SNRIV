#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH -p main -n 1
#SBATCH --mem=16GB
#SBATCH -A u2023011

singularity=/mnt/picea/storage/singularity/seidr_0.14.2.sif
in=$(realpath ../data/seidr/results/aggregate/aggregated.sf)
out=$(realpath ../data/seidr/results/HardThreshold)

# singularity exec -B /mnt:/mnt $singularity seidr view --column-headers \
# -t 0.444 $in -o $out/filtered_0444.tsv

singularity exec -B /mnt:/mnt $singularity seidr view --binary \
-t 0.444 $in -o $out/filtered_0444.sf 
