#!/bin/bash -l

account=u2024010

# process the argument
pgs=/mnt/reference/goldStandard/Potra02_KEGG-based-positive-gold-standard.tsv
ngs=/mnt/reference/goldStandard/Potra02_KEGG-based-negative-gold-standard.tsv

bb=$(realpath data/seidr/results/aggregate/aggregated.sf)
ht=$(realpath data/seidr/results/HardThreshold/thresholded.tsv)
indir=$(realpath data/seidr/results/backbone)
out=$(realpath data/seidr/results/ResRoc)
singularity=$(realpath singularity/seidr_0.14.2.sif)

if [ ! -d $out ]; then
  mkdir -p $out
fi

bbnam=$(basename $bb)
htnam=$(basename $ht)
if [ ! -h $indir/$bbnam ]; then
  ln -sf -t $indir $bb 
fi

if [ ! -h $indir/$htnam ]; then
  ln -sf -t $indir $ht 
fi

# find the network files
# Running with the negative gold standard
for f in $(find $indir  -name "*.tsv" -o -name "*.sf"); do
  fnam=$(basename ${f/.sf/})

  # run the roc on all
  sbatch -A $account \
  -o $out/${fnam}_roc_WithNegative.out -e $out/${fnam}_roc_WithNegative.err \
  pipeline/runSeidrRoc.sh $singularity $f $pgs $ngs \
  $out/${fnam}_roc_WithNegative.tsv
  # Rerun seidr roc using no negative edges information
  sbatch -A $account \
  -o $out/${fnam}_roc_NoNegative.out -e $out/${fnam}_roc_NoNegative.err \
  pipeline/runSeidrRocNoNegative.sh $singularity $f $pgs \
  $out/${fnam}_roc_NoNegative.tsv
done
