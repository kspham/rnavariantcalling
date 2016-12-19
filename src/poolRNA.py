#!/usr/bin.env python

import sys
import os
import logging
import subprocess

def LogInstance(log_file, debug=True):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    logger = logging.getLogger()
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if debug:
	console_log.setLevel(logging.DEBUG)
    else:
	console_log.setLevel(logging.INFO)
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)
    return logger

def exeCommand(sCommand):
    ###Log command data
    #oLogger.debug(sCommand)

    ###Get all output data
    outData, errData = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        close_fds=True).communicate()

def _read_samples(samples_BAMs):
    """
    read uuid, BAM tupple
    :param samples_uuid:
    :param samples_BAMs:
    :return:
    """
    bams = open(samples_BAMs).read().strip().split("\n")
    uuids = [b.split("/")[10] for b in bams]
    return zip(uuids, bams)

def _init_batch_commnads(sample_pairs, binpool, bin, reg, ref, threads, outdir):
    """
    initialize a bunch of variant calling commands
    :param sample_pairs:
    :return:
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    oLogger = LogInstance("./%s/batch_log.txt" %outdir)
    for outFile, bamFile in sample_pairs:
        vcfFile = "%s/%s.vcf" %(outdir, outFile)
        if os.path.exists(vcfFile) and os.stat(vcfFile).st_size!=0:
            oLogger.info("Found %s" %vcfFile)
            continue
        cmd = "python %s --path %s --thread %s --num 50 --reg %s --ref %s --min 2 --bam %s.bam --out ./%s/%s.vcf" %(binpool, bin, threads, reg, ref, outFile, outdir, outFile)
        scpcmd = "scp -r -i /cluster/home/kspham/.ssh/id_rsa_bioturing www-data@same.ucsd.edu:%s %s.bam" %(bamFile, outFile)
        oLogger.info("Downloading bam file %s" %bamFile)
        scpcmd2 = "scp -r -i /cluster/home/kspham/.ssh/id_rsa_bioturing www-data@same.ucsd.edu:%s.bai %s.bam.bai" %(bamFile, outFile)
        if not os.path.exists(bamFile):
            exeCommand(scpcmd)
            exeCommand(scpcmd2)
        oLogger.info("Started %s" %cmd)
        exeCommand(cmd)
        oLogger.info("Finished %s" %outFile)

if __name__ =="__main__":
   samples_BAMs, binpool, bin, reg, ref, threads , outdir= sys.argv[1:]
   pairs = _read_samples(samples_BAMs)
   _init_batch_commnads(pairs, binpool, bin, reg, ref, threads, outdir)
