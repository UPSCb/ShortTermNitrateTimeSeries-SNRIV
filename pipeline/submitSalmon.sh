#!/bin/bash
set -eu 

in=$(realpath ../data/raw)
inx=$(realpath ../data/genome/index/salmon)
out=$(realpath ../data)/salmon/

# singularity
export SINGULARITY_BINDPATH="/mnt:/mnt"
sings=$(realpath ../../single_cell_analysis_poplar/singularity/salmon-1.6.0.simg)

[[ ! -d $out ]] && mkdir -p $out

for f in $(find $in -name "*_S*_L*_R1_001.fastq.gz"); do
  sname=$(basename ${f/_S*_L*_R1_001.fastq.gz/})
  sbatch -A u2022003 -o $out/$sname.out -e $out/$sname.err \
  $(realpath ../UPSCb-common/pipeline/runSalmon.sh) \
  $sings $inx $f ${f/_R1_001.fastq.gz/_R2_001.fastq.gz} $out
done
