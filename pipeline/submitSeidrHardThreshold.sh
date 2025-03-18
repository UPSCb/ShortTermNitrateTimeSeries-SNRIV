#!/bin/bash

#singularity container
singularity=$(realpath singularity/seidr_0.14.2.sif)
account="u2022003"

# process the argument
input="/mnt/picea/projects/aspseq/htuominen/ShortTermNitrateTimeSeries-SNRIV/seidr/results/aggregate/aggregated.sf"
out="/mnt/picea/projects/aspseq/htuominen/ShortTermNitrateTimeSeries-SNRIV/seidr/results/HardThreshold"
mkdir -p $out

# submit
sbatch -A $account -e $out/thresholded.err -o $out/thresholded.out \
  -J thresholdedLD1 pipeline/runSeidrHardThreshold.sh $singularity $input $out/thresholded.tsv
