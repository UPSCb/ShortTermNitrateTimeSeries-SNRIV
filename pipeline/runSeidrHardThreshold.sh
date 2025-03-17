#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -p main -n 2

set -ex

# variables
CPU=2

# helper functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Usage: $0 <seidr file> <output folder>
  
"

# sanity


if [ $# -ne 3 ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing singularity container"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be an existing network file"
fi

if [ ! -d $(dirname $3) ]; then
  abort "The third argument directory needs to exist"
fi

# run
export OMP_NUM_THREADS=$CPU
singularity exec -B /mnt:/mnt $1 seidr threshold -f --in-file $2 -m 0.1 -M 0.9 -o $3

