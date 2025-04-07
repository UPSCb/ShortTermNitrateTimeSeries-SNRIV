#!/bin/bash -l 

out=$(realpath data/seidr/clustering)
[[ ! -d $out ]] && mkdir -p $out

# sbatch -t 12:00:00 -p main -n 20 --mem=16GB -A u2022003 \
# -o $out/cluster_bb9.log -e $out/cluster_bb9.err \
# $(realpath pipeline/runNetworkCluster.sh) 

sbatch -t 12:00:00 -p main -n 20 --mem=16GB -A u2022003 \
-o $out/cluster_bb1.log -e $out/cluster_bb1.err \
$(realpath pipeline/runNetworkCluster.sh)

# sbatch -t 12:00:00 -p main -n 20 --mem=16GB -A u2022003 \
# -o $out/cluster_bb5.log -e $out/cluster_bb5.err \
# $(realpath pipeline/runNetworkCluster.sh) 
