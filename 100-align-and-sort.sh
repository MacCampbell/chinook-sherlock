#!/bin/bash -l

#This script came from Mike Miller and works on the UC Davis farm.
#I'm using this for whole-genome resequencing and think I'll expand the -t and flag below.
#It should produced a sorted and deduplicated bam file for each set of paired sequences.

#Requires a tab delimited list of paired end files (list, $1)
#SRR1613242_1  SRR1613242_2

#Requires a path to indexed reference genome (ref, $2)

#Usage (from data/ dir I suppose.)
#bash ../100-align-and-sort.sh tester.txt $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna.gz

list=$1
ref=$2

wc=$(wc -l ${list} | awk '{print $1}')

x=1
while [ $x -le $wc ] 
do
        string="sed -n ${x}p ${list}" 
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1, $2}')   
        set -- $var
        c1=$1
        c2=$2

       echo "#!/bin/bash -l
       bwa mem $ref ${c1}.fastq.gz ${c2}.fastq.gz | samtools view -Sb | samtools sort - -o ${c1}.sort.bam
       samtools view -f 0x2 -b ${c1}.sort.bam | samtools rmdup - ${c1}.sort.flt.bam
       samtools index ${c1}.sort.flt.bam" > ${c1}.sh
       sbatch -t 10-10:00:00 --mem=8G ${c1}.sh

       x=$(( $x + 1 ))

done

# bwa mem $ref ${c1}.fastq ${c2}.fastq | samtools view -Sb - | samtools sort - ${c1}.sort
#Changed by Mac on 03182020