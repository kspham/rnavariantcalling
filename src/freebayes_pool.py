#!/usr/bin/env python
"""
Created on Sep 16, 2015
author: HoaPT - BioTuring
FreeBayes Thread POOL
"""
__author__    = 'Linuxpham <thaihoabo@gmail.com>'
__contact__   = 'thaihoabo@gmail.com'
__version__   = "1.0"
__date__      = '2016-07-17'
import subprocess
import threadpool
import optparse
import sys
import time
import os

###Parse command line options
parserInstance = optparse.OptionParser(usage="usage: python %prog [options]", version="FreeBayes Pool 1.0 - BioTuring")
parserInstance.add_option("-p", "--path", default="freebayes", help="FreeBayes bin path")
parserInstance.add_option("-r", "--reg", default="chr1:0-1", help="Region ignore chrM file path")
parserInstance.add_option("-m", "--regm", default="chrM:577-647", help="Region chrM file path")
parserInstance.add_option("-c", "--min", default="5", help="Min alternate count")
parserInstance.add_option("-b", "--bam", default="", help="Bam file")
parserInstance.add_option("-f", "--ref", default="", help="Fasta reference file path")
parserInstance.add_option("-t", "--thread", default="4", help="Thread number")
parserInstance.add_option("-d", "--debug", default="0", help="Debug mode")
parserInstance.add_option("-o", "--out", default="/tmp/out.vcf", help="Output VCF file")
parserInstance.add_option("-n", "--num", default="50", help="Region chunk size")
parserInstance.parse_args()

###VCF output file
debugMode = int(parserInstance.values.debug)
freebayesBinPath = str(parserInstance.values.path)
regionFilePath = str(parserInstance.values.reg)
regionChrMFilePath = str(parserInstance.values.regm)
minAlterbateCount = int(parserInstance.values.min)
bamFilePath = str(parserInstance.values.bam)
refFilePath = str(parserInstance.values.ref)
outFilePath = str(parserInstance.values.out)
threadCount = int(parserInstance.values.thread)
numCount = int(parserInstance.values.num)

###Open VCF output file
outFileTmpPath = "%s.tmp" % (outFilePath)
if os.path.exists(outFileTmpPath):
    os.remove(outFileTmpPath)
outTmpFile = open(outFileTmpPath, 'w')

###Execute the command
def executeCommand(sCommand):
    global outTmpFile
    if (debugMode > 0):
        print("Command: %s" % (sCommand))
    subprocess.Popen(sCommand, shell=True, stdout=outTmpFile, stderr=subprocess.STDERR, close_fds=True).communicate()

###Create freebayes command
def createCommand(prefixCommand, arrRegionTMP, ignoreHeader):
    sCommand = ""
    if(len(arrRegionTMP) > 0):
        sCommand = prefixCommand + " --region %s" % (','.join(list(arrRegionTMP)))
        if (ignoreHeader == True):
            sCommand = sCommand + " --no-header"
    return sCommand

###Logic business
def main():
    global outTmpFile
    if ((len(regionFilePath) == 0) or (len(refFilePath) == 0) or (len(outFilePath) == 0)):
        parserInstance.print_help()
        sys.exit(1)

    iStartTime = time.time()
    arrListParam = []
    prefixCommand = "%s -f %s -C %d %s" % (freebayesBinPath, refFilePath, minAlterbateCount, bamFilePath)

    if (debugMode > 0):
        print("Create VCF HEADER")

    arrRegionTMP = set()
    arrFirstRegionTMP = set()
    arrFirstRegionTMP.add('chr1:0-1')
    sCommand = createCommand(prefixCommand, arrFirstRegionTMP, False)
    if len(sCommand) > 0:
        arrListParam.append(sCommand)

    if (debugMode > 0):
        print("Prepare chrM region for freebayes")

    ###USE for CHR-M
    oRegionChrMFile = open(regionChrMFilePath)
    for line in oRegionChrMFile:
        line = line.strip()
        if (len(arrRegionTMP) == 13):
            sCommand = createCommand(prefixCommand, arrRegionTMP, True)
            if(len(sCommand) > 0):
                arrListParam.append(sCommand)
            arrRegionTMP.clear()
        else:
            if (len(line) > 0):
                arrRegionTMP.add(line)
    oRegionChrMFile.close()

    if (len(arrRegionTMP) > 0):
        sCommand = createCommand(prefixCommand, arrRegionTMP, True)
        if (len(sCommand) > 0):
            arrListParam.append(sCommand)
        arrRegionTMP.clear()

    if (debugMode > 0):
        print("Prepare other region for freebayes")

    ###USE for other CHRM
    oRegionFile = open(regionFilePath)
    for line in oRegionFile:
        line = line.strip()
        if (len(arrRegionTMP) == numCount):
            sCommand = createCommand(prefixCommand, arrRegionTMP, True)
            if (len(sCommand) > 0):
                arrListParam.append(sCommand)
            arrRegionTMP.clear()
        else:
            if(len(line) > 0):
                arrRegionTMP.add(line)
    oRegionFile.close()

    if(len(arrRegionTMP) > 0):
        sCommand = createCommand(prefixCommand, arrRegionTMP, True)
        if (len(sCommand) > 0):
            arrListParam.append(sCommand)
        arrRegionTMP.clear()

    if (debugMode > 0):
        print("Execute all %d freebayes command in POOL" % (len(arrListParam)))

    if (len(arrListParam) > 0):
        oFastQCPool = threadpool.ThreadPool(threadCount)
        arrRequest = threadpool.makeRequests(executeCommand, arrListParam)
        [oFastQCPool.putRequest(req) for req in arrRequest]
        oFastQCPool.wait()

    ###Close file to SYNC DISK
    sys.stdout.flush()
    sleep(2)
    outTmpFile.close()

    ###START VALID VCF file########
    ###Check output file data
    if os.path.exists(outFilePath):
        os.remove(outFilePath)
    outFile = open(outFilePath, 'w')
    ###Write the FIRST header
    with open(outFileTmpPath, 'r') as oReadFile:
        for lineData in oReadFile:
            lineData = lineData.strip()
            if(len(lineData) == 0):
                continue
            if lineData[0:2] != "##":
                break
            outFile.write(lineData)
            outFile.write("\n")
        oReadFile.close()
    ###Write the valid VCF line
    with open(outFileTmpPath, 'r') as oReadFile:
        for lineData in oReadFile:
            lineData = lineData.strip()
            if (len(lineData) == 0):
                continue
            if lineData[0:2] != "chr":
                continue
            outFile.write(lineData)
            outFile.write("\n")
        oReadFile.close()
    outFile.close()
    ###END VALID VCF file########

    ###Debug information
    iEndTime = time.time()
    print("Finish freebayes parallel in %d seconds" % (int(iEndTime) - int(iStartTime)))

if __name__ == "__main__":
    sys.exit(main())
