# shernook
Chinook WGS for SHERLOCK assays   
 
### General Plan

__1__  Align short read data to the reference genome (GCF_002872995.1_Otsh_v1.0) with bwa mem     
__2__  Sort short read alignments and remove PCR duplicates with samtools     
__3__  Realign indels with GATK "Indel Realigner"   
__4__  Identify run associated genetic variation with ANGSD doAsso     
__5__  Create a VCF that defintiely includes indels to find absolute allele frequency differences
__6__  Create some absolute frequency difference measurements of all variants around SNPs after imputation

### Scripts should generally follow the plan and are numbered consecutively

100-align-and-sort.sh       
101-indel-realignment.sh - looks like the tool has changed.      
102-compute-coverage.sh     

etc.


#### Color Scheme
Winter&Spring=purple AA3377      

Fall&LateFall=yellow CCBB44      

Winter-cyan 66CCEE     

Spring-green 228833     

Fall-red EE6677      

LateFall-grey BBBBBB     


### Thompson, Anderson et al. data    

`maccamp@farm:~/data/chinook-wgs-raw$ wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA667732' -O - | tee SraRunInfo.csv`

Locations of files:
https://sra-download.ncbi.nlm.nih.gov/traces/sra67/SRR/012498/SRR12798274
https://sra-download.ncbi.nlm.nih.gov/traces/sra37/SRR/012498/SRR12798433


Field 10 of SraRunInfo.csv

`maccamp@farm:~/data/chinook-wgs-raw$ cut -f 10 -d ',' SraRunInfo.csv | while read line; do wget $line; done;`

SraRunInfo.csv in /meta/ now.     

Now we need to dump.

Working comand:     
srun -p high -t 24:00:00 fastq-dump --outdir ./ --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR11744922 > dump.out 2>dump.err &     
Need to speed this all up...


----Bash Script
list=$1 

#Say we are reading in something tab delimited.
IFS=$'\t'

while read name
do
   echo $name
   
  echo "#!/bin/bash -l
  fastq-dump --outdir ./ --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip $name " > $name.sh
  sbatch -p high -t 2-0:00:00 --mem=8G $name.sh

done < $list
 ----
