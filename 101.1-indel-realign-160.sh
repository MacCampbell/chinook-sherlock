#! /bin/bash

#for samples of the form DPCh_plate2_H02_S154.rmdup.bam
# maccamp@farm:~/google/chinook_WGS_processed$ bash ~/shernook/101.1-indel-realign-160.sh 160.list 
# $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna
# deduplicated but contain improperly paired reads?

list=$1 
ref=$2

#Say we are reading in something tab delimited.
IFS=$'\t'

while read name
do
   echo $name
   
  echo "#!/bin/bash -l
  samtools view -f 0x2 -b $name.rmdup.bam -o $name.sort.flt.bam
  gatk LeftAlignIndels --disable-tool-default-read-filters true -R $ref -I $name.sort.flt.bam -O $name.bam
  samtools index $name.bam
  samtools depth -a -b $HOME/shernook/meta/GCF_002872995.1.bed $name.bam | awk '{sum+="\$3"} END {print sum/NR}' > $name.cov" > $name-gatk.sh
  sbatch -t 7-0:00:00 --mem=8G $name-gatk.sh

done < $list
