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
ignoreChrM = True

###Execute the command
def executeCommand(sCommand):
    global outFilePath
    arrCommand = sCommand.split("|||")
    outFileTmpPath = "%s.%s.tmp" % (outFilePath, arrCommand[1])
    if os.path.exists(outFileTmpPath):
        os.remove(outFileTmpPath)
    outTmpFile = open(outFileTmpPath, 'w')
    if (debugMode > 0):
        print("Command: %s" % (arrCommand[0]))
    subprocess.Popen(arrCommand[0], shell=True, stdout=outTmpFile, stderr=subprocess.STDOUT, close_fds=True).communicate()
    sys.stdout.flush()
    outTmpFile.close()

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
    global ignoreChrM
    if ((len(regionFilePath) == 0) or (len(refFilePath) == 0) or (len(outFilePath) == 0)):
        parserInstance.print_help()
        sys.exit(1)

    iStartTime = time.time()
    arrListParam = []
    jFileIndex = 0
    arrRegionTMP = set()
    arrFirstRegionTMP = set()
    prefixCommand = "%s -f %s -C %d %s" % (freebayesBinPath, refFilePath, minAlterbateCount, bamFilePath)

    ###USE for CHR1-0-1
    if (debugMode > 0):
        print("Create VCF HEADER")
    arrFirstRegionTMP.add('chr1:0-1')
    sCommand = createCommand(prefixCommand, arrFirstRegionTMP, False)
    if len(sCommand) > 0:
        sFullCommand = "%s|||%s" % (sCommand, str(jFileIndex))
        jFileIndex = jFileIndex + 1
        arrListParam.append(sFullCommand)

    ###USE for CHR-M
    if ignoreChrM == False:
        if (debugMode > 0):
            print("Prepare chrM region for freebayes")

        oRegionChrMFile = open(regionChrMFilePath)
        for line in oRegionChrMFile:
            line = line.strip()
            if (len(arrRegionTMP) == 13):
                sCommand = createCommand(prefixCommand, arrRegionTMP, True)
                if(len(sCommand) > 0):
                    sFullCommand = "%s|||%s" % (sCommand, str(jFileIndex))
                    jFileIndex = jFileIndex + 1
                    arrListParam.append(sFullCommand)
                arrRegionTMP.clear()
            else:
                if (len(line) > 0):
                    arrRegionTMP.add(line)
        oRegionChrMFile.close()

        if (len(arrRegionTMP) > 0):
            sCommand = createCommand(prefixCommand, arrRegionTMP, True)
            if (len(sCommand) > 0):
                sFullCommand = "%s|||%s" % (sCommand, str(jFileIndex))
                jFileIndex = jFileIndex + 1
                arrListParam.append(sFullCommand)
            arrRegionTMP.clear()

    ###USE for other CHRM
    if (debugMode > 0):
        print("Prepare other region for freebayes")
    oRegionFile = open(regionFilePath)
    for line in oRegionFile:
        line = line.strip()
        if (len(arrRegionTMP) == numCount):
            sCommand = createCommand(prefixCommand, arrRegionTMP, True)
            if (len(sCommand) > 0):
                sFullCommand = "%s|||%s" % (sCommand, str(jFileIndex))
                jFileIndex = jFileIndex + 1
                arrListParam.append(sFullCommand)
            arrRegionTMP.clear()
        else:
            if(len(line) > 0):
                arrRegionTMP.add(line)
    oRegionFile.close()

    if(len(arrRegionTMP) > 0):
        sCommand = createCommand(prefixCommand, arrRegionTMP, True)
        if (len(sCommand) > 0):
            sFullCommand = "%s|||%s" % (sCommand, str(jFileIndex))
            jFileIndex = jFileIndex + 1
            arrListParam.append(sFullCommand)
        arrRegionTMP.clear()

    if (debugMode > 0):
        print("Execute all %d freebayes command in POOL" % (len(arrListParam)))

    if (len(arrListParam) > 0):
        oFastQCPool = threadpool.ThreadPool(threadCount)
        arrRequest = threadpool.makeRequests(executeCommand, arrListParam)
        [oFastQCPool.putRequest(req) for req in arrRequest]
        oFastQCPool.wait()

    ###START VALID VCF file########
    if os.path.exists(outFilePath):
        os.remove(outFilePath)
    outFile = open(outFilePath, 'w')
    for iFileIndex in range(0, jFileIndex):
        outFileTmpPath = "%s.%s.tmp" % (outFilePath, str(iFileIndex))
        with open(outFileTmpPath, 'r') as oReadFile:
            for lineData in oReadFile:
                lineData = lineData.strip()
                if (len(lineData) == 0):
                    continue
                elif lineData[0:1] == "#":
                    outFile.write(lineData)
                    outFile.write("\n")
                elif lineData[0:3] == "chr":
                    arrLineData = lineData.split("\t")
                    if len(arrLineData) > 6:
                        arrChrID = arrLineData[0].split("chr")
                        outFile.write(arrChrID[1])
                        outFile.write("\t")
                        outFile.write("\t".join(arrLineData[1:]))
                        outFile.write("\n")
            oReadFile.close()
            os.remove(outFileTmpPath)
    outFile.close()
    ###END VALID VCF file########

    ###Debug information
    iEndTime = time.time()
    if (debugMode > 0):
        print("Finish freebayes parallel in %d seconds (%d files)" % ((int(iEndTime) - int(iStartTime)), jFileIndex))

if __name__ == "__main__":
    sys.exit(main())
