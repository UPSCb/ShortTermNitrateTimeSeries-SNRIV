#!/bin/bash

mkdir -p data/seidr/results/aggregate

paths=$(realpath data/seidr/results/*/*.sf)
sbatch pipeline/seidr-aggegate-micro.sh data/seidr/results/aggregate $paths

