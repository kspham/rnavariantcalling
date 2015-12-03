# RNAvariantcalling

RNAvariantcalling is a powerful tool for calling variants from RNA sequences to provide high quality variants, which were filter by a recently-published database about editing sites in human genome. 

  - RNA sequences were mapped by 2 mapper STAR and HISAT2, which can handle splice junction problem very well
  - Variants were called by Freebayes, then only variants called from both STAR and HISAT2 were stored for the next step
  - Such reliable variants were filtered one more time with a database of humnan editing sites

### Version
1.0
### Author
ptdtan@gmail.com

### Dependencies

RNAvariantcalling uses a number of materials to work properly:

* [hg19 genome indexed by STAR] 
* [hg19 genome indexed by HISAT2] 
* [hg19 genome un-indexed] - for variant calling
* [hg19 gene annotation file]
* Python 2.7+
* vcftools
* STAR aligner
* HISAT2 aligner
* Freebayes variant calling

These will be downloaded automatically at the very first time you run rnavariantcalling.

And also, rnavariantcalling need these tools below have to be installed on your computer:
* samtools 
* tabix
* bgzip

### Installation
```
1. $git clone https://github.com/kspham/rnavariantcalling.git
2. $cd rnavariantcalling
3. $python setup.py install
4. $./configure
At the very first time you run rnavariantcalling, you have to download the indexed human genome for STAR Aligner and HISAT2 Aligner, and also the gene annotation of human genome.

3. $./initial.sh

4.Done! Now you rnavariantcalling can work properly.

```
### Usage 

```
$ rnavariantcalling [-h] [--ThreadsN N] [--reads read1 read2 ... readN] [--outdir OUTDIR] --config yamlFile
```
    -h, --helps help message and exit 
    --reads READS [READS ...], -U READS [READS ...] Input RNA reads 
    ---outdir OUTDIR, -o OUTDIR Where the final result will be stored 
    --ThreadsN N Number of threads
    --config yamlFile configuration file as yaml format
### Example
    - Paired-ends reads: 
```
$ rnavariantcalling.py --ThreadsN 32 --reads 1.fastq.gz 2.fastq.gz -o /output/directory/ --config /path/to/your/config.yaml
```

    -Single-end read: 
```
$ rnavariantcalling.py --ThreadsN 32 -U unpaired_read.fastq.gz -o /output/directory/ --config /path/to/your/config.yaml
```

License
----

MIT


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [hg19 genome indexed by STAR]: <https://www.encodeproject.org/files/ENCFF069ZCO/@@download/ENCFF069ZCO.tar.gz>
   [hg19 genome indexed by HISAT2]: <ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch37.tar.gz>
   [hg19 genome un-indexed]:<http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit>
   [hg19 gene annotation file]:<ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz>

