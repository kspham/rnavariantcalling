#!/cluster/home/kspham/.local/bin/python

""" Pipeline calling variant from RNA-seq"""
import subprocess
import threadpool
import optparse
import sys
import time
import os

###Logic business
def main():
    iStartTime = time.time()
    vcfFile = "/mnt/hdd2/hoapt/test_freebayes/seqlib.4.vcf.bk"
    vcfOutFile = "/mnt/hdd2/hoapt/test_freebayes/seqlib.4.vcf"

    ###START VALID VCF file########
    outFile = open(vcfOutFile, 'w')

    ###Write the FIRST header
    with open(vcfFile, 'r') as oReadFile:
        for lineData in oReadFile:
            lineData = lineData.strip()
            if(len(lineData) == 0):
                continue
            if lineData[0:1] != "#":
                break
            outFile.write(lineData)
            outFile.write("\n")
        oReadFile.close()

    ###Write the valid VCF line
    with open(vcfFile, 'r') as oReadFile:
        for lineData in oReadFile:
            lineData = lineData.strip()
            if (len(lineData) == 0):
                continue
            if lineData[0:3] != "chr":
                continue
            arrData = lineData.split("\t")
            if len(arrData) > 6:
                arrChrID = arrData[0].split("chr")
                outFile.write(arrChrID[1])
                outFile.write("\t")
                outFile.write("\t".join(arrData[1:]))
                outFile.write("\n")
        oReadFile.close()
    outFile.close()
    ###END VALID VCF file########

    ###Debug information
    iEndTime = time.time()
    print("Finish in %d seconds" % (int(iEndTime) - int(iStartTime)))

if __name__ == "__main__":
    sys.exit(main())