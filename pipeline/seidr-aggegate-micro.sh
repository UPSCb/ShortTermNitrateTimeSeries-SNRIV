#!/bin/bash -l 
#SBATCH -A u2024010
#SBATCH -t 12:00:00
#SBATCH -n 28
#SBATCH -o data/seidr/results/aggregate/aggregate.out
#SBATCH -e data/seidr/results/aggregate/aggregate.err
#SBATCH -J aggregate

set -ex
DIRECTIONALITY="-k"
FORCE=
METHOD="-m irp"

# usage
USAGETXT=\
"
  Usage: $0 [options] <out dir> <sf file> [sf file] ... [sf file]

  Options:
          -f  force overwrite output
          -k  keep the directionality, default: true
          -m  method, default: irp
"

CPU=28

# source
source ${SLURM_SUBMIT_DIR:-$(pwd)}/UPSCb-common/src/bash/functions.sh

# singularity container
singularity="/mnt/picea/storage/singularity/seidr_0.14.2.sif"

# Get the options
while getopts fkm: option
do
    case "$option" in
        f) FORCE="-f";;
        k) DIRECTIONALITY=;;
        k) METHOD="-m $OPTARG";;
        \?) ## unknown flag
		    abort;;
    esac
done
shift `expr $OPTIND - 1`

OPTIONS="$FORCE $DIRECTIONALITY $METHOD"

if [ $# -lt 2 ]; then
  abort "This script expects at least 2 arguments"
fi

if [ ! -d $1 ]; then
  abort "The first argument needs to be an existing directory"
fi

# run
cd $1
shift
export OMP_NUM_THREADS=$CPU
#rm -f aggregated.sf

singularity exec -B /mnt:/mnt $singularity seidr aggregate $OPTIONS -O $CPU $@
