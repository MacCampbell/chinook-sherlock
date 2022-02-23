# Thinking up a way to compute tab delimited coverage....


# Can consider DepthOfCoverage (BETA) [options] -L meta/GCF_002872995.1.bed

# This would be a one liner
#gatk \
#   DepthOfCoverage \
# --disable-tool-default-read-filters \
#   -R reference.fasta \
#   -O file_name_base \
#  -I input_bams.list
#   -L ../meta/GCF_002872995.1.bed

# maccamp@farm:~/shernook/filter-sorted$ cat first-batch.txt | perl -pe 's/\n/.bam\n/g' > first-batch-bams.txt 

# srun -t 00:08:00 -p med --mem=8G gatk DepthOfCoverage -R $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna -O first-batch -I first-batch-bams.txt -L ../meta/GCF_002872995.1.bed
# Some sort of disagreement between contigs.  *.bams apparently don't have contig names (?)
# maccamp@farm:~/shernook/filter-sorted$ samtools view -H F020137COL_5_S27_L004_R1_001.bam  | grep NC
#@SQ	SN:NC_037097.1	LN:96198142
#@SQ	SN:NC_037098.1	LN:57406636
#@SQ	SN:NC_037099.1	LN:81079776

# A USER ERROR has occurred: Input files reference and reads have incompatible contigs: No overlapping contigs found.
#  reference contigs = [NC_037097.1, NC_037098.1, NC_037099.1, 
#  reads contigs = []

## For existing data
#maccamp@farm:~/google/chinook_WGS_processed$ ls | grep rmdup | grep -v bai | perl -pe 's/.bam//g' > list.txt
#maccamp@farm:~/google/chinook_WGS_processed$ bash ~/shernook/102-compute-coverage.sh list.txt ../meta/GCF_002872995.1.bed


# samtools depth -a -b bedFileOfRegions -o outfile infile
# bedfile is meta/GCF_002872995.1.bed

# USAGE
# ../102-compute-coverage.sh first-batch.txt ../meta/GCF_002872995.1.bed

list=$1 
bedf=$2

#Say we are reading in something tab delimited.
IFS=$'\t'

while read name
do
   echo $name
   
  echo "#!/bin/bash -l
  samtools depth -a -b $bedf $name.bam | awk '{sum+="\$3"} END {print sum/NR}' > $name.cov" > $name-depth.sh
  sbatch -t 1:00:00 --mem=8G $name-depth.sh

done < $list

