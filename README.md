# shernook
Chinook WGS for SHERLOCK assays

# Color Scheme
Winter&Spring=purple AA3377      

Fall&LateFall=yellow CCBB44      

Winter-cyan 66CCEE     

Spring-green 228833     

Fall-red EE6677      

LateFall-grey BBBBBB     

 
### General Plan

__1__  Align short read data to the reference genome (GCF_002872995.1_Otsh_v1.0) with bwa mem     
__2__  Sort short read alignments and remove PCR duplicates with samtools     
__3__  Realign indels with GATK "Indel Realigner"   
__4__  Identify run associated genetic variation with ANGSD doAsso     
__5__  Create a VCF that defintiely includes indels to find absolute allele frequency differences
__6__  Create some absolute frequency difference measurements of all variants around SNPs identified by ANGSD     

### Tester
I'm going to get some Chinook WGS data. There are some samples available on SRA from AK, SRX733554: Oncorhynchus tshawytscha Genome sequencing SRR1613247 & SRR1613242.     

`wget --recursive --no-parent --execute robots=off --no-host-directories ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/002/SRR1613242/`    

`wget --recursive --no-parent --execute robots=off --no-host-directories ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/007/SRR1613247/`     

### Thompson, Anderson et al. data    
maccamp@farm:~/data/chinook-wgs-raw$ wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA667732' -O - | tee SraRunInfo.csv

Locations of files:
https://sra-download.ncbi.nlm.nih.gov/traces/sra67/SRR/012498/SRR12798274
https://sra-download.ncbi.nlm.nih.gov/traces/sra37/SRR/012498/SRR12798433


Field 10 of SraRunInfo.csv

maccamp@farm:~/data/chinook-wgs-raw$ cut -f 10 -d ',' SraRunInfo.csv | while read line; do wget $line; done;

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
### Scripts should generally follow the plan in term of numbering.

100-align-and-sort.sh       
101-indel-realignment.sh - looks like the tool has changed.      
102-compute-coverage.sh     


### Short workup notes
While waiting for all the mapping to complete, I ran a "first-batch.txt" set through.
It has 80 samples. I'm creating a new subdir "filter-sorted" to work in with symlinks to try not to screw everything up.
`maccamp@farm:~/shernook/filter-sorted$ ls | grep -v bai | perl -pe 's/.sort.flt.bam//g' > first-batch.txt`    
`maccamp@farm:~/shernook/filter-sorted$ ../101-indel-realignment.sh first-batch.txt  $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna`     
Seems to work.

Also a couple files were very large and the clock time was not specified right. So I have to rerun 
W011924USR_41_S70_L004_R1_001
F080016TOU_2_S63_L004_R1_001
Increased run time to 10 days here.

`maccamp@farm:~/shernook/data$ bash ../100-align-and-sort.sh batch1-redo.txt $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna.gz`

In filter sorted, I'm am re-analyzing "second-batch.txt to include the winter fish we have that finished mapping.
`maccamp@farm:~/shernook/filter-sorted$ ls | grep sort.flt.bam | grep -v bai |  perl -pe 's/.sort.flt.bam//g' > second-batch.txt`

We are now up to 90 samples.
`maccamp@farm:~/shernook/filter-sorted$ ../101-indel-realignment.sh second-batch.txt $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna`     
`maccamp@farm:~/shernook/filter-sorted$ bash ../102-compute-coverage.sh second-batch.txt ../meta/GCF_002872995.1.bed`


Now all samples have completed mapping. To include the final six samples:    
F080016TOU_2_S63_L004_R1_001 
F010240MER_19_S61_L004_R1_001
S980004YUB_43_S99_L004_R1_001
S990001YUB_43_S54_L004_R1_001
S020035MIL_43_S51_L004_R1_001
W011924USR_41_S70_L004_R1_001

`../101-indel-realignment.sh third-batch.txt $HOME/genomes/chinook/GCF_002872995.1_Otsh_v1.0_genomic.fna`   
` bash ../102-compute-coverage.sh third-batch.txt ../meta/GCF_002872995.1.bed`

## How about those Thompson et al. seqs?
I'm trying first an analysis of the existing *.bams.    
We don't need the Klamath fish I suppose.
