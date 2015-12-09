#!/usr/bin/env python

""" Pipeline calling variant from RNA-seq"""

import os,sys
import argparse
import subprocess
import yaml
import uuid
import random
import hashlib
import logging
import logging.handlers

#MD5 function
def MD5string(strData):
	m=hashlib.md5()
	m.update(strData.encode('UTF-8'))
	return m.hexdigest()
	
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
	oLogger.debug("Done aligment with STAR")

#Mapping with HISAT2
def HISAT2_mapping(reads, N, output, pairend):
	os.chdir(HISAT2out)
	if pairend:
		exeCommand(shellEscape(' '.join([hisat2, "--threads", N, "-q", "-x", "genome", "-1", reads[0], "-2", reads[1], "-S", 
output])))
	else:
		exeCommand(shellEscape(' '.join([hisat2, "--threads", N, "-q", "-x", "genome", "-U", ' '.join([read for read in reads]), 
"-S", output])))
	oLogger.debug("Done aligment with HISAT2")

def Variant_Calling(bam, dir, threads):
	exeCommand(shellEscape(' '.join(["multithread.py",REF, freebayes, dir+"/"+bam, CHRO,threads, dir+"/"])))
	oLogger.debug("Done calling variant for:"+bam)

def filter1():	
	exeCommand(shellEscape(' '.join(["filter", vcftools, HISAT2out, STARout, TEMP])))
	oLogger.debug("Done filtering stage 1")

def filter2():
	exeCommand(shellEscape(' '.join(["filter2", vcftools, vcfconcat,TEMP, uname, output, editsite])))
	oLogger.debug("Done filtering stage 2")

def snpEff():
	os.chdir(output)
	exeCommand(' '.join(["java","-Xmx10g","-jar", SnpEff, "-v",
"GRCh37.75",uname+".recode.vcf",">",uname+"ann.vcf"]))
	oLogger.debug("Done annotation!")

def snpSift():
	os.chdir(output)
	exeCommand(' '.join(["sed 's/^chr//'",uname+"ann.vcf",">","tmp","&& mv tmp",uname+"ann.vcf"]))
	exeCommand(' '.join(["java","-jar", SnpSift, "annotate", "-id",
vcfdatabase,uname+"ann.vcf",">",uname+"annotated.vcf"]))
	oLogger.debug("Added rsID")

def ParsingBAM(N):
	#Index Aligned.sortByCoord.out.bam
	os.chdir(STARout)
	exeCommand(shellEscape(' '.join([SAMBAMBA,"index -t",N, "Aligned.sortedByCoord.out.bam"])))
	oLogger.debug("Indexed:Aligned.sortByCoord.out.bam")
	#Reheader and sort, index HISAT2.Aligned
	os.chdir(HISAT2out)
	output="HISAT2.Aligned"
	exeCommand(shellEscape(' '.join([SAMBAMBA, "view -S -f bam -t", N, output, ">", output+".bam"])))
	oLogger.debug("Convert %s to %s" %(output, output+".bam"))
        exeCommand(shellEscape(' '.join(["samtools reheader",header, output+".bam","> temp", "&& mv temp", output+".bam"])))
	oLogger.debug("Reheader BAM file")
        exeCommand(shellEscape(' '.join([SAMBAMBA,"sort -t",N,"-o", output+".sorted.bam", output+".bam"])))
	oLogger.debug("Sort %s" %(output+".bam"))
        exeCommand(shellEscape(' '.join([SAMBAMBA, "index -t",N, output+".sorted.bam"])))
	oLogger.debug("Indexed: %s" %(output+"sorted.bam.bai"))

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
	snpeff=cfg['tools']['snpeff']
	editsite=cfg['lib']['editsite']
	perl=cfg['lib']['PERL5LIB']
	java=cfg['tools']['java']
	SnpEff=cfg['tools']['snpeff']
	SnpSift=cfg['tools']['snpsift']
	vcfdatabase=cfg['lib']['vcfdatabase']

	os.environ['HISAT2_INDEXES']=HISAT2ref
	os.environ['PERL5LIB']=perl

        #Get reads abs path
        reads = args.reads
        if '.gz' in reads[0]:
                iszipped = True
        else:
             	iszipped = False

        for i in range(len(reads)):
                reads[i] = os.path.abspath(reads[i])

        # Generate job name
        uname = MD5string(reads[0])

	#Initial Logging

        oLogger=logging.getLogger("RNAvariantcalling")
        oLogger.setLevel(logging.DEBUG)
        oLoggerHandler=logging.handlers.RotatingFileHandler("./"+uname+".log",maxBytes=10485760, backupCount=10)
	oFormatter=logging.Formatter("%(levelname)s:[%(asctime)s]-[%(filename)s at line (%(lineno)d) of (%(funcName)s) function]-[%(message)s]")
	oLoggerHandler.setFormatter(oFormatter)
	oLogger.addHandler(oLoggerHandler)

	oLogger.debug("Job name as MD5 hash of file name:"+uname)
	oLogger.debug("Reads input: "+str(reads))

	# Get output directory
        if args.outdir:
                output = os.path.abspath(args.outdir)
        else:
                output = os.path.abspath(os.environ['PWD'])

	#Main program

	TEMP+=uname
	STARout+=uname
	HISAT2out+=uname
	job=temporary+uname
	final=output+"/"+uname+".recode.vcf"
	try:
		os.makedirs(TEMP)	
		os.makedirs(STARout)
		os.makedirs(HISAT2out)
	except OSError:
		pass
	try:
		os.makedirs(output)
	except OSError:
		pass
	oLogger.debug("Making TEMP STARout and HISAT2out directories"+"%s\n%s\n%s\n%s"%(TEMP,STARout,HISAT2out, output))

	command = {}
	
	command[1] = [STAR_mapping,[reads, iszipped, args.ThreadsN, STARref]]

	command[2]=[HISAT2_mapping,[reads, args.ThreadsN,"HISAT2.Aligned", len(reads)>1]]

	command[3]=[ParsingBAM,[args.ThreadsN]]	

	command[4]=[Variant_Calling,["Aligned.sortedByCoord.out.bam", STARout, args.ThreadsN]]

	command[5]=[Variant_Calling,["HISAT2.Aligned.sorted.bam", HISAT2out, args.ThreadsN]]

	command[6]=[filter1,[]]

	command[7]=[filter2,[]]

	command[8]=[snpEff,[]]
	
	command[9]=[snpSift,[]]
	oLogger.debug("Setup commands")

	# Initial steps
	stepsDone = {}
	for i in range(1,10):
		stepsDone[i] = "False"
	oLogger.debug("Setup tasks pool")

	# Get done steps
	try:
		steps = open(job,"r+")
		for line in steps:
			l = line.split()
			stepsDone[int(l[0])] = l[1]
	except IOError:
		steps = open(job,"w")
	
	oLogger.debug("Get job status" + str(stepsDone))

	for step in range (1,10):
		if stepsDone[step] == "False":
			params = command[step][1]
			command[step][0](*params)
			steps.write(str(step)+"\t"+"True\n")
	steps.close()
	oLogger.debug("Sucessfully! "+"Job name: "+uname+" File name: "+uname+".recode.vcf")
