# rnavariantcalling
Tan

Requirement: 
python 2.7+


#usage: 
               pipe.py [-h] [--reads READS [READS ...]] [--outdir OUTDIR]
               [-r1 1.fastq.gz] [-r2 2.fastq.gz] [--ThreadsN N]

#Automatically generate SNPs variant for given RNA short reads

#optional arguments:
  -h, --help            show this help message and exit
  --reads READS [READS ...], -U READS [READS ...]
                        Input RNA unpaired reads
  --outdir OUTDIR, -o OUTDIR
                        Where the final result will be stored
  -r1 1.fastq.gz        First reads (In case of paired-end reads)
  -r2 2.fastq.gz        Second reads
  --ThreadsN N          Number of threads
