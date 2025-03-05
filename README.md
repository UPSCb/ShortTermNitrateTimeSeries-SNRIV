# ShortTermNitrateTimeSeries-SNRIV

## Abstract
Hybrid aspen was grown under controlled nitrogen conditions. These trees were treated with nitrate and samples were taken for RNA sequencing at specific time points to trace gene expression changes in wood.

## Setup
```bash
ln -s /mnt/picea/projects/aspseq/htuominen/ShortTermNitrateTimeSeries-SNRIV data
ln -s /mnt/picea/storage/singularity/ .
```


## Nextflow

```bash
 nextflow run nf-core/rnaseq -r 3.14.0 \
  -profile upscb,singularity -work-dir data/work \
  -params-file nextflow/nf-params.json -c nextflow/upscb.config
```


