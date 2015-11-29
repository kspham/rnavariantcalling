#!/usr/bin/env python

""" Pipeline calling variant from RNA-seq"""

import os,sys
import argparse
import subprocess


TEMP=os.environ["TEMP"]
STARout =os.environ["STARout"]
HISAT2out = os.environ["HISAT2out"]
STARref =os.environ["STARREF"]
HISAT2ref=os.environ["HISAT2REF"]
REF =os.environ["REF"]
SRC = os.environ["SRC"]

#How to execute a command
def exeCommand(sCommand):

###Get all output data
	outData, errData = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, close_fds=True).communicate()

###Get all response data
	for lineData in outData.splitlines():
		#if(self.RUNNING_DEBUG_FLAG == 1):
		outStringData = str(lineData)
		print("%s" % (outStringData))
###If there is error
	if((errData != None) and (len(errData) > 0)):
		print("Command has error:{0}".format(errData))

def shellEscape(s):
	return s.replace("(","\(").replace(")","\)")

#Mapping with STAR	
def STAR_mapping(reads, ReadIsGzipped, N, dir):
	exeCommand(shellEscape(' '.join(["STAR", "--runThreadN",N, "--genomeDir", dir, "--readFilesIn", 
		' '.join([read for read in reads]), "--alignIntronMin", "20", "--alignIntronMax", "500000", "--outFilterMismatchNmax", "10", 
"--outSAMtype", "BAM", "SortedByCoordinate", ''.join(["--readFilesCommand gunzip -c" for i in range(1) if ReadIsGzipped])])))
	exeCommand(shellEscape(' '.join(["$SAMBAMBA index -t 32", "Aligned.sortedByCoord.out.bam"])))


#Mapping with HISAT2
def HISAT2_mapping(reads, N, output, pairend):
	if pairend:
		exeCommand(shellEscape(' '.join(["hisat2", "--threads", N, "-q", "-x", "genome", "-1", reads[0], "-2", reads[1], "-S", 
output])))
	else:
		exeCommand(shellEscape(' '.join(["hisat2", "--threads", N, "-q", "-x", "genome", "-U", ' '.join([read for read in reads]), 
"-S", output])))

	exeCommand(shellEscape(' '.join(["$SAMBAMBA view -S -f bam -t 32", output, ">", output+".bam"])))
	exeCommand(shellEscape(' '.join(["$SAMBAMBA sort -t 32","-o", output+".sorted.bam", output+".bam"])))
	#exeCommand(shellEscape(' '.join(["$SAMBAMBA view -H", output+".sorted.bam", "| sed 's/SN:/SN:chr/g' /dev/stdin > header"])))
	exeCommand(shellEscape(' '.join(["samtools reheader $HEADER", output+".sorted.bam","> temp", "&& mv temp", output+".sorted.bam"])))
	exeCommand(shellEscape(' '.join(["$SAMBAMBA index -t 32", output+".sorted.bam"])))
 
def Variant_Calling(bam1, bam2, dir1, dir2, threads):
	for bam, dir in [(bam1, dir1), (bam2, dir2)]:
		exeCommand(shellEscape(' '.join(["multithread.py","$REF", "freebayes", dir+"/"+bam, "$CHRO",threads, 
dir+"/"])))

def filter(output):
	
	exeCommand(shellEscape("mkdir $TEMP"))	

	#Stage 1
	exeCommand(shellEscape(' '.join(["sh", "$FILTER", "vcftools", "$HISAT2out", "$STARout"])))

	#Stage 2 merging
	listFile = os.listdir(TEMP)
	os.chdir(TEMP)
	for f in listFile:
		exeCommand(shellEscape(' '.join(["bgzip -c", f, ">", f+".gz"])))
		exeCommand(shellEscape(' '.join(["tabix -p", "vcf", f+".gz"])))
	exeCommand(shellEscape(' '.join(["vcf-merge", ' '.join([f+".gz" for f in listFile]), "| bgzip -c >", "human_variant.vcf.gz"])))

	#Stage 3
	exeCommand(shellEscape(' '.join(["vcftools", "--gzvcf","human_variant.vcf.gz", "--exclude-positions","$LIB/human_edit.txt", 
"--recode","--recode-INFO-all", "--out", "final"])))
	
	exeCommand(shellEscape("mv final.recode.vcf "+output))
	exeCommand(shellEscape("rm -fr $TEMP"))

if __name__ == '__main__':

	parser=argparse.ArgumentParser(description='Automatically generate SNPs variant for give RNA short reads')
	parser.add_argument('--reads', '-U', type=str,help='Input RNA reads paired or unpaired', nargs='+', required=True)
	parser.add_argument('--outdir', '-o',type=str, help='Where the final result will be stored')
	parser.add_argument('--ThreadsN', metavar='N', type=str, help='Number of threads', default='4')
	args=parser.parse_args()

	reads = args.reads

	if '.gz' in reads[0]:
		iszipped = True
	else:
		iszipped = False

	for i in range(len(reads)):
		reads[i] = os.path.abspath(reads[i])

	exeCommand(shellEscape("mkdir "+args.outdir))
	output = os.path.abspath(args.outdir)

	os.chdir(STARout)
	STAR_mapping(reads, iszipped, args.ThreadsN, STARref)

	os.chdir(HISAT2out)
	HISAT2_mapping(reads, args.ThreadsN,"HISAT2.Aligned", len(reads)>1)

	Variant_Calling("Aligned.sortedByCoord.out.bam", "HISAT2.Aligned.sorted.bam", STARout, HISAT2out, args.ThreadsN)
	
	filter(output)

