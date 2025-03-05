#!/bin/bash

headless=$(realpath data/seidr/headless.tsv)
genes=$(realpath data/seidr/genes.tsv)

bash pipeline/seidr-wrapper-micro.sh \
/mnt/picea/storage/singularity/seidr_0.14.2.sif \
$headless $genes
