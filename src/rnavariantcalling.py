#!/cluster/home/kspham/.local/bin/python

""" Pipeline calling variant from RNA-seq"""

import os, sys
import argparse
import subprocess
import yaml
import hashlib
import logging
import logging.handlers

# MD5 function
def MD5stringdata(strData):
    m = hashlib.md5()
    m.update(strData.encode('UTF-8'))
    return m.hexdigest()


#How to execute a command
def exeCommand(sCommand):
    ###Log command data
    oLogger.debug(sCommand)

    ###Get all output data
    outData, errData = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        close_fds=True).communicate()

    ###Get all response data
    for lineData in outData.splitlines():
        #outStringData = str(lineData)
        #print("%s" % (outStringData))
        continue

    ###If there is error
    if ((errData != None) and (len(errData) > 0)):
        oLogger.error("Command has error:{0}".format(errData))


def shellEscape(s):
    return s.replace("(", "\(").replace(")", "\)")


#Mapping with STAR
def STAR_mapping(reads, ReadIsGzipped, N, dir):
    oLogger.debug("Run STAR mapping with genome reference %s" %(dir))
    os.chdir(STARout)
    exeCommand(shellEscape(' '.join([STAR, "--runThreadN", N, "--genomeDir", dir, "--readFilesIn",
                                     ' '.join([read for read in reads]), "--alignIntronMin", "20", "--alignIntronMax",
                                     "500000", "--outFilterMismatchNmax", "10",
                                     "--outSAMtype", "BAM", "SortedByCoordinate",
                                     ''.join(["--readFilesCommand gunzip -c" for i in range(1) if ReadIsGzipped])])))
    oLogger.debug("Done aligment with STAR")


#Mapping with HISAT2
def HISAT2_mapping(reads, N, output, pairend, onlySTAR):
    if onlySTAR:
        pass
    else:
        oLogger.debug("Run HISAT2 mapping with genome reference %s" %(HISAT2ref))
        os.chdir(HISAT2out)
        if pairend:
            exeCommand(
                shellEscape(' '.join([hisat2, "--threads", N, "-q", "-x", "genome", "-1", reads[0], "-2", reads[1], "-S",
                                      output])))
        else:
            exeCommand(shellEscape(
                ' '.join([hisat2, "--threads", N, "-q", "-x", "genome", "-U", ' '.join([read for read in reads]),
                          "-S", output])))
        oLogger.debug("Done aligment with HISAT2")


def Variant_Calling(bam, dir, threads):
    ###Please put fasta_generate_regions.py in the current directory
    os.chdir(STARout)
    command = ' '.join(["freebayes-parallel",region,threads,
    "-f", REF, bam, ">",dir+"/"+uname+".vcf"])
    exeCommand(shellEscape(command))
    oLogger.debug(command)
    oLogger.debug("Done calling variant for:" + bam)


def filter1(onlySTAR):
    if not onlySTAR:
        exeCommand(shellEscape(' '.join(["filter", vcftools, HISAT2out, STARout, TEMP, output])))
    else:
        exeCommand(shellEscape(' '.join(["cp", STARout+"/"+uname+".vcf", TEMP+"/"+uname+".recode.vcf"])))
    oLogger.debug("Done filtering stage 1")


def filter2():
    exeCommand(shellEscape(' '.join(["filter2", vcftools, vcfconcat, TEMP, uname, output, editsite])))
    oLogger.debug("Done filtering stage 2")


def snpEff(ref, names={'hg19':'GRCh37.75', 'mm10':'GRCm38.82'}):

    os.chdir(output)
    exeCommand(' '.join([java, "-d64 -Xms4G -Xmx8G -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit", "-jar", SnpEff, "-v",
                         names[ref], uname + ".recode.vcf", ">", uname + "ann.vcf"]))
    oLogger.debug("Done annotation!")


def snpSift():
    os.chdir(output)
    exeCommand(' '.join(["sed 's/^chr//'", uname + "ann.vcf", ">", "tmp", "&& mv tmp", uname + "ann.vcf"]))
    exeCommand(' '.join(
        [java, "-d64 -Xms4G -Xmx8G -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit", "-jar", SnpSift, "annotate", "-id",
         vcfdatabase, uname + "ann.vcf", ">", uname + "annotated.vcf"]))
    exeCommand(' '.join(['bgzip -c', uname+"annotated.vcf", ">", uname+"annotated.vcf.gz"]))
    exeCommand(' '.join(['tabix -p', uname+"annotated.vcf.gz"]))
    oLogger.debug("Added rsID")


def ParsingBAM(N, onlySTAR):

    #Index Aligned.sortByCoord.out.bam
    os.chdir(STARout)
    exeCommand(shellEscape(' '.join([SAMBAMBA, "index -t", N, "Aligned.sortedByCoord.out.bam"])))
    oLogger.debug("Indexed:Aligned.sortByCoord.out.bam")

    #Convert and sort, index HISAT2.Aligned
    if not onlySTAR:
        os.chdir(HISAT2out)
        output = "HISAT2.Aligned"
        exeCommand(shellEscape(' '.join([SAMBAMBA, "view -S -f bam -t", N, output, ">", output + ".bam"])))
        oLogger.debug("Convert %s to %s" % (output, output + ".bam"))

        exeCommand(shellEscape(' '.join([SAMBAMBA, "sort -t", N, "-o", output + ".sorted.bam", output + ".bam"])))
        oLogger.debug("Sort %s" % (output + ".bam"))
        exeCommand(shellEscape(' '.join([SAMBAMBA, "index -t", N, output + ".sorted.bam"])))
        oLogger.debug("Indexed: %s" % (output + "sorted.bam.bai"))

def cleanBam(starDir, hisat2Dir):
    #exeCommand(shellEscape(' '.join(["rm", "-f", hisat2Dir + "/HISAT2.Aligned"])))
    #exeCommand(shellEscape(' '.join(["rm", "-f", hisat2Dir + "/HISAT2.Aligned.bam"])))
    exeCommand(shellEscape(' '.join(["cp", STARout+"/Aligned.sortedByCoord.out.bam*", output])))
    exeCommand(shellEscape(' '.join(["rm", "-rf", starDir])))
    #exeCommand(shellEscape(' '.join(["rm", "-rf", hisat2Dir + "/chr*.vcf"])))
    oLogger.debug("Clean:HISAT2.Aligned, HISAT2.Aligned.bam")

###Main FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatically generate SNPs variant for give RNA short reads')
    parser.add_argument('--reads', '-U', type=str, help='Input RNA reads paired or unpaired', nargs='+', required=True)
    parser.add_argument('--outdir', '-o', type=str, help='Where the final result will be stored')
    parser.add_argument('--ThreadsN', metavar='N', type=str, help='Number of threads', default='4')
    parser.add_argument('--config', metavar='yamlFile', type=str, help='Config file as yaml format', required=True)
    parser.add_argument('--set', metavar='Steps will be set Done', type=str, nargs='+')
    parser.add_argument('--unset', metavar='Steps will be set NOT done yet', type=str, nargs='+')
    parser.add_argument('--logdir', help='Logging directory', type=str)
    parser.add_argument('--species', '-s', type=str, help='hg19/mm10',required=True)
    parser.add_argument('--vcfdatabase', '-v', type=str, help='vcf database for annotation')
    parser.add_argument('--onlySTAR', help='only run with STAR_mapping, not HISAT2_mapping', dest='onlySTAR', action='store_true')
    parser.set_defaults(onlySTAR=True)
    args = parser.parse_args()


    #Parse yaml file:
    with open(args.config, "r") as ymlfile:
        cfg = yaml.load(ymlfile)
        ymlfile.close()

    #headers = cfg['lib']['headers']
    TEMP = cfg['folder']['tmp']
    temporary = cfg['folder']['temporary']
    STARout = cfg['folder']['STARout']
    HISAT2out = cfg['folder']['HISAT2out']

    #hg19
    if args.species == 'hg19':
        STARref = cfg['lib']['hg19STARref']
        HISAT2ref = cfg['lib']['hg19HISAT2ref']
        REF = cfg['lib']['hg19REF']
        CHRO = cfg['lib']['hg19chro']
        editsite = cfg['lib']['hg19editsite']
        #vcfdatabase = cfg['lib']['hg19vcfdatabase']
        region=cfg['lib']['hg19region']
    else:
        #mm10
        STARref = cfg['lib']['mm10STARref']
        HISAT2ref = cfg['lib']['mm10HISAT2ref']
        REF = cfg['lib']['mm10REF']
        CHRO = cfg['lib']['mm10chro']
        editsite = cfg['lib']['mm10editsite']
        region=cfg['lib']['mm10region']
        #vcfdatabase = cfg['lib']['mm10vcfdatabase']

    STAR = cfg['tools']['STAR']
    hisat2 = cfg['tools']['hisat2']
    vcftools = cfg['tools']['vcftools']
    vcfconcat = cfg['tools']['vcfconcat']
    freebayes = cfg['tools']['freebayes']
    SAMBAMBA = cfg['tools']['sambamba']
    snpeff = cfg['tools']['snpeff']
    perl = cfg['lib']['PERL5LIB']
    java=cfg['tools']['java']
    SnpEff = cfg['tools']['snpeff']
    SnpSift = cfg['tools']['snpsift']
    if args.vcfdatabase:
        vcfdatabase = args.vcfdatabase
    else:
        vcfdatabase = cfg['lib'][args.species+'vcfdatabase']

    os.environ['HISAT2_INDEXES'] = HISAT2ref
    os.environ['PERL5LIB'] = perl

    #Get reads abs path
    reads = args.reads
    if '.gz' in reads[0]:
        iszipped = True
    else:
        iszipped = False

    for i in range(len(reads)):
        reads[i] = os.path.abspath(reads[i])

    # Generate job name
    uname = MD5stringdata(reads[0])

    #Initial Logging
    if args.logdir:
        logdir = os.path.abspath(args.logdir)
    else:
        logdir = "."
    oLogger = logging.getLogger("RNAvariantcalling")
    oLogger.setLevel(logging.DEBUG)
    oLoggerHandler = logging.handlers.RotatingFileHandler(logdir + "/" + uname + ".log", maxBytes=10485760, backupCount=10)
    oLoggerStreamHandler = logging.StreamHandler()
    oLoggerStreamHandler.setLevel(logging.DEBUG)

    oFormatter = logging.Formatter(
        "%(levelname)s:[%(asctime)s]-[%(filename)s at line (%(lineno)d) of (%(funcName)s) function]-[%(message)s]")
    oLoggerHandler.setFormatter(oFormatter)
    oLoggerStreamHandler.setFormatter(oFormatter)
    oLogger.addHandler(oLoggerHandler)
    oLogger.addHandler(oLoggerStreamHandler)

    oLogger.debug("Job name as MD5 hash of file name:" + uname)
    oLogger.debug("Reads input: " + str(reads))

    # Get output directory
    if args.outdir:
        output = os.path.abspath(args.outdir)
    else:
        output = os.path.abspath(os.environ['PWD'])

    #Main program

    TEMP += uname
    STARout += uname
    HISAT2out += uname
    job = temporary + uname
    final = output + "/" + uname + ".recode.vcf"
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
    oLogger.debug("Making TEMP STARout and HISAT2out directories" + "%s\n%s\n%s\n%s" % (TEMP, STARout, HISAT2out, output))

    command = {}

    command[1] = [STAR_mapping, [reads, iszipped, args.ThreadsN, STARref]]

    command[2] = [HISAT2_mapping, [reads, args.ThreadsN, "HISAT2.Aligned", len(reads) > 1, args.onlySTAR ]]

    command[3] = [ParsingBAM, [args.ThreadsN, args.onlySTAR,]]

    command[4] = [Variant_Calling, ["Aligned.sortedByCoord.out.bam", STARout, args.ThreadsN]]

    command[5] = [Variant_Calling, ["HISAT2.Aligned.sorted.bam", HISAT2out, args.ThreadsN]]

    command[6] = [filter1, [args.onlySTAR]]

    command[7] = [filter2, []]

    command[8] = [snpEff, [args.species]]

    command[9] = [snpSift, []]
    command[10] = [cleanBam, [STARout, HISAT2out]]
    oLogger.debug("Setup commands")

    # Initial steps
    stepsDone = {}
    for i in range(1, len(command)+1):
        stepsDone[i] = "False"
    oLogger.debug("Setup tasks pool")

    # Get done steps
    try:
        steps = open(job, "r+")
        for line in steps:
            l = line.split()
            stepsDone[int(l[0])] = l[1]
    except IOError:
        steps = open(job, "w")

    # Get done steps from user input
    if args.set:
        for step in args.set:
            stepsDone[int(step)] = "True"

    # Get NOT done yet steps from user input
    if args.unset:
        for step in args.unset:
            stepsDone[int(step)] = "False"

    # Set command 5 is False due to args.onlySTAR
    if args.onlySTAR:
        stepsDone[5] = "True"

    oLogger.debug("Get job status" + str(stepsDone))

    for step in range(1, len(command)+1):
        if stepsDone[step] == "False":
            params = command[step][1]
            command[step][0](*params)
            steps.write(str(step) + "\t" + "True\n")
    steps.close()
    oLogger.debug("Sucessfully! " + "Job name: " + uname + " File name: " + uname + ".recode.vcf")
