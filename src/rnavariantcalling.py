#!/usr/bin/env python

""" Pipeline calling variant from RNA-seq"""

import os,sys
import argparse
import subprocess
import yaml
import uuid
import random
import hashlib

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
	os.chdir(STARout)
	exeCommand(shellEscape(' '.join([STAR, "--runThreadN",N, "--genomeDir", dir, "--readFilesIn", 
		' '.join([read for read in reads]), "--alignIntronMin", "20", "--alignIntronMax", "500000", "--outFilterMismatchNmax", "10", 
"--outSAMtype", "BAM", "SortedByCoordinate", ''.join(["--readFilesCommand gunzip -c" for i in range(1) if ReadIsGzipped])])))


#Mapping with HISAT2
def HISAT2_mapping(reads, N, output, pairend):
	os.chdir(HISAT2out)
	if pairend:
		exeCommand(shellEscape(' '.join([hisat2, "--threads", N, "-q", "-x", "genome", "-1", reads[0], "-2", reads[1], "-S", 
output])))
	else:
		exeCommand(shellEscape(' '.join([hisat2, "--threads", N, "-q", "-x", "genome", "-U", ' '.join([read for read in reads]), 
"-S", output])))

def Variant_Calling(bam, dir, threads):
	exeCommand(shellEscape(' '.join(["multithread.py",REF, freebayes, dir+"/"+bam, CHRO,threads, dir+"/"])))

def filter(output):	
	#Stage 1
	exeCommand(shellEscape(' '.join(["filter", vcftools, HISAT2out, STARout, TEMP])))
	#Stage 2
	exeCommand(shellEscape(' '.join(["filter2", vcftools, vcfconcat,TEMP, uname, output, editsite])))

def ParsingBAM(N):
	#Index Aligned.sortByCoord.out.bam
	os.chdir(STARout)
	exeCommand(shellEscape(' '.join([SAMBAMBA,"index -t",N, "Aligned.sortedByCoord.out.bam"])))
	#Reheader and sort, index HISAT2.Aligned
	os.chdir(HISAT2out)
	output="HISAT2.Aligned"
	exeCommand(shellEscape(' '.join([SAMBAMBA, "view -S -f bam -t", N, output, ">", output+".bam"])))
        exeCommand(shellEscape(' '.join(["samtools reheader",header, output+".bam","> temp", "&& mv temp", output+".bam"])))
        exeCommand(shellEscape(' '.join([SAMBAMBA,"sort -t",N,"-o", output+".sorted.bam", output+".bam"])))
        exeCommand(shellEscape(' '.join([SAMBAMBA, "index -t",N, output+".sorted.bam"])))

if __name__ == '__main__':

	parser=argparse.ArgumentParser(description='Automatically generate SNPs variant for give RNA short reads')
	parser.add_argument('--reads', '-U', type=str,help='Input RNA reads paired or unpaired', nargs='+', required=True)
	parser.add_argument('--outdir', '-o',type=str, help='Where the final result will be stored')
	parser.add_argument('--ThreadsN', metavar='N', type=str, help='Number of threads', default='4')
	parser.add_argument('--config', metavar='yamlFile', type=str, help='Config file as yaml format', required=True) 
	args=parser.parse_args()

	#Parse yaml file:
	with open(args.config, "r") as ymlfile:
        	cfg=yaml.load(ymlfile)
	header = cfg['lib']['header']
	TEMP= cfg['folder']['tmp']
	temporary=cfg['folder']['temporary']
	STARout = cfg['folder']['STARout']
	HISAT2out = cfg['folder']['HISAT2out']
	STARref = cfg['lib']['STARref']
	HISAT2ref= cfg['lib']['HISAT2ref']
	REF = cfg['lib']['REF']
	CHRO = cfg['lib']['chro']
	STAR= cfg['tools']['STAR']
	hisat2= cfg['tools']['hisat2']
	vcftools=cfg['tools']['vcftools']
	vcfconcat=cfg['tools']['vcfconcat']
	freebayes=cfg['tools']['freebayes']
	SAMBAMBA= cfg['tools']['sambamba']
	editsite=cfg['lib']['editsite']
	perl=cfg['lib']['PERL5LIB']
	
	os.environ['HISAT2_INDEXES']=HISAT2ref
	os.environ['PERL5LIB']=perl

        reads = args.reads
	print reads
	if '.gz' in reads[0]:
		iszipped = True
	else:
		iszipped = False

	for i in range(len(reads)):
		reads[i] = os.path.abspath(reads[i])

	if args.outdir:
		output = os.path.abspath(args.outdir)
	else:
		output = os.path.abspath(os.environ['PWD'])
	os.mkdir(output)

	# Generate UUID

	uname = hashlib.md5(reads[0]).hexdigest()

	#Main program

	TEMP+=uname
	STARout+=uname
	HISAT2out+=uname
	job=temporary+uname

	os.makedirs(TEMP)	
	os.makedirs(STARout)
	os.makedirs(HISAT2out)
	command = {}
	
	command[1] = [STAR_mapping,[reads, iszipped, args.ThreadsN, STARref]]

	command[2]=[HISAT2_mapping,[reads, args.ThreadsN,"HISAT2.Aligned", len(reads)>1]]

	command[3]=[ParsingBAM,[args.ThreadsN]]	

	command[4]=[Variant_Calling,["Aligned.sortedByCoord.out.bam", STARout, args.ThreadsN]]

	command[5]=[Variant_Calling,["HISAT2.Aligned.sorted.bam", HISAT2out, args.ThreadsN]]

	command[6]=[filter,[output]]

	# Initial steps
	stepsDone = {}
	for i in range(1,7):
		stepsDone[i] = "False"

	# Get done steps
	steps = open(job,"a")
	for line in steps:
		l = line.split()
		stepsDone[int(l[0])] = l[1]

	for step in range (1,7):
		if stepsDone[step] == "False":
			params = command[step][1]
			command[step][0](*params)
			steps.write(str(step)+"\t"+"True\n")
	
	print "Sucessfully!", "Job name: "+uname, "File name: "+uname+".recode.vcf"

	os.rmdir(TEMP)
	os.rmdir(STARout)
	os.rmdir(HISAT2out)
