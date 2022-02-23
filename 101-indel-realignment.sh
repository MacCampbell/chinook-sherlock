# I should make this a proper script for the farm
# Previously, this worked for me.

# indexed salmo
#samtools faidx salmo.fna
#then 
#java -jar /home/mac/picard/build/libs/picard.jar CreateSequenceDictionary R=salmo.fna O=salmo.dict

#/usr/local/bin/gatk/GenomeAnalysisTK.jar
#$ for i in $(ls -d sample*/); do echo $i;cd $i; java -Xmx2g -jar /usr/local/bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../salmo.fna -I $(basename ${i}).rmdup.addgp.bam -o $(basename ${i}).rmdup.addgp.intervals;java -Xmx4g -jar /usr/local/bin/gatk/GenomeAnalysisTK.jar -I $(basename ${i}).rmdup.addgp.bam -R ../salmo.fna -T IndelRealigner -targetIntervals $(basename ${i}).rmdup.addgp.intervals  -o $(basename ${i}).realigned.bam;cd ..;done

### Now, it looks kuje gatk 4.1.8.1 does something else
# /home/maccamp/gatk-4.1.8.1/gatk

# Created a .dict file
# maccamp@farm:~/genomes/chinook$ gatk CreateSequenceDictionary  -R GCF_002872995.1_Otsh_v1.0_genomic.fna

# Working command
# gatk LeftAlignIndels -R ~/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna -I W011934USR_41_S72_L004_R1_001.sort.flt.bam -O W011934USR_41_S72_L004_R1_001.bam
# output empty, all reads filtered out, disabling read filter
# gatk LeftAlignIndels --disable-tool-default-read-filters true -R ~/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna -I W011934USR_41_S72_L004_R1_001.sort.flt.bam -O W011934USR_41_S72_L004_R1_001.bam
# now have negative index error, is bug
#[Exists in version 4.1.7.0 and 4.1.8.1 ]
# [No issue in version 4.1.4.1 ]
# upgrading to 4.1.9.0


#gatk LeftAlignIndels \
#   -R reference.fasta \
#   -I input.bam \
#   -O output.bam

# We can copy the 100-align-and-sort.sh script basic


# Works like so ../101-indel-realignment.sh first-batch.txt  $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna

list=$1 
ref=$2

#Say we are reading in something tab delimited.
IFS=$'\t'

while read name
do
   echo $name
   
  echo "#!/bin/bash -l
  gatk LeftAlignIndels --disable-tool-default-read-filters true -R $ref -I $name.sort.flt.bam -O $name.bam
  samtools index $name.bam" > $name-gatk.sh
  sbatch -t 2-0:00:00 --mem=8G $name-gatk.sh

done < $list

