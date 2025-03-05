#!/bin/bash

singularity=$(realpath singularity/seidr_0.14.2.sif)
headless=$(realpath data/seidr/headless.tsv)
genes=$(realpath data/seidr/genes.tsv)

bash pipeline/seidr-wrapper-micro.sh $singularity $headless $genes
