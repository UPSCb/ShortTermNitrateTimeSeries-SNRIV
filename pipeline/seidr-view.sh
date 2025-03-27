#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH -p main -n 1
#SBATCH --mem=16GB
#SBATCH -A u2022003
#SBATCH --output data/seidr/results/seidr-view.log
#SBATCH --error data/seidr/results/seidr-view.err

singularity=$(realpath singularity/seidr_0.14.2.sif)


for thres in 1 5 9
do singularity exec -B /mnt:/mnt $singularity seidr view --column-headers \
data/seidr/results/backbone/backbone-"$thres"-percent.sf -o data/seidr/results/backbone/backbone-"$thres"-percent.tsv
done
