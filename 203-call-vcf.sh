#! /bin/bash -l

# Creating a .vcf with bcftools
# Bounding roughly around peak of associaiton with NC_037124.1 like so:
#-r NC_037124.1:9233657-13595086 
# Right now just trying out at the command line like so


bcftools mpileup -Ou -f reference.fa alignments.bam
filtering quality: bcftools view -i '%QUAL>=20' calls.bcf


#Previously, I did do something like this:
#i call variants
#bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Oz -o calls.vcf.gz
#bcftools index calls.vcf.gz

#ii normalize indels
#bcftools norm -f reference.fa calls.vcf.gz -Ob -o calls.norm.bcf

#iii filter adjacent indels within 5bp
#bcftools filter –IndelGap 5 calls.norm.bcf -Ob -o calls.norm.flt-indels.bcf

srun -p med -t 8:00:00 --nodes=1 bcftools mpileup -Ou --fasta-ref ~/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna -r NC_037124.1:9233657-13595086 -b bamlists/bamlist56.bamlist | \
bcftools call -mv -Ob -o outputs/200/bamlist56-chrom28.bcf

srun -p med -t 8:00:00 --nodes=1 bcftools mpileup -Ou --fasta-ref ~/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna -r NC_037112.1:22424576-28717394 -b bamlists/bamlist56.bamlist | \
bcftools call -mv -Ob -o outputs/202/bamlist56-chrom16.bcf