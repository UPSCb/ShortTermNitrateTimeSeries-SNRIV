#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH -p main -n 20
#SBATCH --mem=16GB
#SBATCH -A u2022003

singularity=$(realpath singularity/seidr_0.14.2.sif)
in=$(realpath data/seidr/results/aggregate/aggregated.sf)
out=$(realpath data/seidr/results/aggregate)

singularity exec -B /mnt:/mnt $singularity seidr view --column-headers \
--filter "irp >= 0.444" $in -o $out/filtered_0.444.tsv
