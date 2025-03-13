#!/bin/bash -l

account=u2022003

# define the thresholds
thresholds=( 2.33 2.05 1.88 1.75 1.64 1.55 1.48 1.41 1.34 1.28 )

# process the argument
network=$(realpath ../data/seidr/results/aggregate/aggregated.sf)
out=$(realpath ../data/seidr/results/backbone)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
for i in {0..9}; do
  j=$(expr $i + 1)
  sbatch -A $account \
  -o $out/backbone-${j}-percent.out \
  -e $out/backbone-${j}-percent.err -J bb-${j}  \
  ./runSeidrBackbone.sh $network ${thresholds[$i]} \
  $out/backbone-${j}-percent.sf
done
