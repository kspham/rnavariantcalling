""" Pipeline calling variant from RNA-seq"""

import os,sys
import argparse

FILE = os.path.dirname(os.path.realpath(__file__))

TEMP= FILE +"/../tmp"
STAR = FILE + "/../bin/STARaligner/STAR"
HISAT2_BUILD = FILE + "/../bin/hisat2-2.0.1-beta/hisat2-build"
HISAT2 = FILE + "/../bin/hisat2-2.0.1-beta/hisat2"
HISAT2out = FILE +"/../HISAT2"
STARout = FILE +"/../STAR"
VCFPERL = FILE +"/../bin/vcftools_0.1.13/perl"
VCFTOOLS =  FILE+"/../bin/vcftools_0.1.13/bin/vcftools"
VCF_MERGE = FILE+"/../bin/vcftools_0.1.13/bin/vcf-merge"
FREEBAYES = FILE + "/../bin/freebayes/bin/freebayes"
SRC = FILE+"/"
CHRO = FILE +"/../lib/chr_coordinates.txt"
HUMAN_EDITING_SITES = FILE+"/../lib/human_edit.txt"
GTF = FILE+"/../lib/Homo_sapiens.GRCh37.75.gtf"
REFERENCE = FILE+"/../reference"
REF = FILE+"/../lib/GRCh37.fa"
SCRIPTS = FILE+"/../scripts/script.sh"
os.system("HISAT2_INDEXES="+HISAT2out)

"""def STAR_Genome_Generate(dir, ref, gtf, N, overhang):
	os.system(' '.join([STAR, "--runMode", "genomeGenerate", "--runThreadN", N, "--genomeDir", dir, "--genomeFastaFiles", 
		' '.join([f for f in ref]), "--sjdbGTFfile", gtf, "--sjdbOverhang", overhang]))"""
	
def STAR_mapping(reads, ReadIsGzipped, N, dir):
	os.system(' '.join([STAR, "--runThreadN",N, "--genomeDir", dir, "--readFilesIn", 
		' '.join([read for read in reads]), "--alignIntronMin", "20", "--alignIntronMax", "500000", "--outFilterMismatchNmax", "10", 
"--outSAMtype", "BAM", "SortedByCoordinate", ''.join(["--readFilesCommand gunzip -c" for i in range(1) if ReadIsGzipped])])) 
	os.system(' '.join(["samtools index", "Aligned.sortedByCoord.out.bam"]))

"""def HISAT2_Index(ref, prefix):
	os.system(' '.join([HISAT2_BUILD, "-p", N, ref, prefix]))"""	

def HISAT2_mapping(reads, N, output, pairend):
	if pairend:
		os.system(' '.join([HISAT2, "--threads", N, "-q", "-x", "humanref", "-1", reads[0], "-2", reads[1], "-S", output]))
	else:
		os.system(' '.join([HISAT2, "--threads", N, "-q", "-x", "humanref", "-U", ' '.join([read for read in reads]), "-S", output]))
	#os.system(' '.join(["mv", output, HISAT2]))
	os.system(' '.join(["samtools view -Sb", output, "| samtools sort -o -", output+".sorted"]))
	os.system(' '.join(["samtools index", output+".sorted.bam"]))
 
def Variant_Calling(bam1, bam2, dir1, dir2, threads):
	for bam, dir in [(bam1, dir1), (bam2, dir2)]:
		os.system(' '.join(["python2.7",SRC+"multithread.py",REF, FREEBAYES, dir+"/"+bam, CHRO,threads, dir+"/"]))

def filter(output):
	
	os.system("export PERL5LIB="+VCFPERL)
	os.system("mkdir "+TEMP)	

	#Stage 1
	os.system(' '.join(["sh", SCRIPTS, VCFTOOLS, HISAT2out, STARout, TEMP]))

	#Stage 2 merging
	listFile = os.listdir(TEMP)
	os.chdir(TEMP)
	for f in listFile:
		os.system(' '.join(["bgzip -c", f, ">", f+".gz"]))
		os.system(' '.join(["tabix -p", "vcf", f+".gz"]))
	os.system(' '.join([VCF_MERGE, ' '.join([f+".gz" for f in listFile]), "| bgzip -c >", "human_variant.vcf.gz"]))

	#Stage 3
	os.system(' '.join([VCFTOOLS, "--gzvcf","human_variant.vcf.gz", "--exclude-positions"," ../lib/human_edit.txt", 
"--recode","--recode-INFO-all", "--out", "final"]))
	
	os.system("mv final.recode.vcf "+output)
	os.system("rm -fr "+TEMP)

if __name__ == '__main__':

	parser=argparse.ArgumentParser(description='Automatically generate SNPs variant for give RNA short reads')
	parser.add_argument('--reads', '-U', type=str,help='Input RNA unpaired reads', nargs='+')
	parser.add_argument('--outdir', '-o',type=str, help='Where the final result will be stored')
	parser.add_argument('-r1',metavar='1.fastq.gz', type=str, help='First reads')
	parser.add_argument('-r2',metavar='2.fastq.gz', type=str, help='Second reads')
	parser.add_argument('--ThreadsN', metavar='N', type=str, help='Number of threads')
	args=parser.parse_args()

	if not args.r1:
		reads = args.reads
	else:
		reads = [args.r1] + [args.r2]

	if '.gz' in reads[0]:
		iszipped = True
	else:
		iszipped = False

	for i in range(len(reads)):
		reads[i] = os.path.abspath(FILE+'/../'+reads[i])
	print reads
	
	#os.system("mkdir reference")
	#STAR_Genome_Generate("reference", [ref], GTF, threads, "100")
	
	os.chdir(STARout)
	STAR_mapping(reads, iszipped, args.ThreadsN, REFERENCE)

	os.chdir(HISAT2out)
	HISAT2_mapping(reads, args.ThreadsN,"HISAT2.Aligned", args.r1 != None)
	
	Variant_Calling("Aligned.sortedByCoord.out.bam", "HISAT2.Aligned", STARout, HISAT2out, args.ThreadsN)
	
	filter(args.outdir)

