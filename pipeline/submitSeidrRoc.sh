#!/bin/bash -l

account=u2024010


# process the argument
pgs=CHANGEME
ngs=CHANGEME
# example: 
# pgs=$(realpath ../goldStandard/Picea-abies_KEGG-based-positive-gold-standard.tsv)
# ngs=$(realpath ../goldStandard/Picea-abies_KEGG-based-negative-gold-standard.tsv)

bb=$(realpath data/seidr/results/aggregate/aggregated.sf)
indir=$(realpath data/seidr/results/backbone)
out=$(realpath data/seidr/results/roc)
singularity=$(realpath singularity/seidr_0.14.2.sif)

if [ ! -d $out ]; then
  mkdir -p $out
fi

bbnam=$(basename $bb)
if [ ! -h $indir/$bbnam ]; then
  ln -sf -t $indir $bb 
fi

# find the network files
for f in $(find $indir -name "*.sf"); do
  fnam=$(basename ${f/.sf/})

  # run the roc on all
  sbatch -A $account --mail-user=$email \
  -o $out/${fnam}_roc.out -e $out/${fnam}_roc.err \
  pipeline/runSeidrRoc.sh $singularity $f $pgs $ngs \
  $out/${fnam}_roc.tsv
done
