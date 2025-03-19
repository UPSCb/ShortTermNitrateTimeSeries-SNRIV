#!/bin/bash

mkdir -p data/seidr/results/aggregate


paths=$(realpath data/seidr/results/*/*.sf | grep -v aggregated.sf)
sbatch /mnt/picea/home/schoudhary/shruti/ShortTermNitrateTimeSeries-SNRIV/pipeline/seidr-aggregate-micro.sh -f data/seidr/results/aggregate $paths

