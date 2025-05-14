mkdir /pfs/stor10/users/home/s/shruti/project/ShortTermNitrateTimeSeries-SNRIV/data/blastml bioinfo-tools

ml GCC/13.2.0  OpenMPI/4.1.6 SeqKit/2.8.2
seqkit grep -r -p nitResponsiveTFsShoot.txt -f nitResponsiveTFShoot.fasta

ml GCC/13.2.0  OpenMPI/4.1.6
ml BLAST+/2.16.0
cd /pfs/stor10/users/home/s/shruti/project/ShortTermNitrateTimeSeries-SNRIV/data/blast
make blastdb -in Potra02_proteins.fasta -out Potra02
blastp -db Potra02 -dbtype prot -query nitResponsiveTFShoot.fasta \
-out nitResponsiveTFPotra.txt -num_threads 4 -qcov_hsp_perc 50 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' 

