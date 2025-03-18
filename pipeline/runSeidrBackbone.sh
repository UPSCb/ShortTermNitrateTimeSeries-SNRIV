#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 2
#SBATCH --partition main

set -ex

# variables
# we use 2 because of hyperthreading
# to force the use of a virtual core instead of a logical one, we could
# sbatch -n 1 -c 1 OMP_NUM_THREADS=1 seidr backbone -O 1
CPU=2


# helper functions
source /mnt/picea/home/schoudhary/shruti/ShortTermNitrateTimeSeries-SNRIV/UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Usage: $0 <seidr file> <threshold> <output filename>
  
  Note: the threshold is the quantile value from a normal distribution,
  so a backbone of 1% is qnorm(0.99) = 2.33. 10% is 1.28. etc.
"

# sanity

if [ $# -ne 4 ]; then
  abort "This script expects 4 arguments"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be an existing file"
fi

if [ ! -d $(dirname $4) ]; then
  abort "The fourth argument directory needs to exist"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing singularity container"
fi
# run
export OMP_NUM_THREADS=$CPU
singularity exec -B /mnt:/mnt $1 seidr backbone -f -F $3 -o $4 $2

