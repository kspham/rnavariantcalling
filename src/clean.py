__author__ = 'linuxpham'
import os

#How to execute a command
def exeCommand(sCommand):
    ###Get all output data
    outData, errData = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        close_fds=True).communicate()

    ###Get all response data
    for lineData in outData.splitlines():
        #if(self.RUNNING_DEBUG_FLAG == 1):
        outStringData = str(lineData)
        print("%s" % (outStringData))
    ###If there is error
    if ((errData != None) and (len(errData) > 0)):
        print("Command has error:{0}".format(errData))


def shellEscape(s):
    return s.replace("(", "\(").replace(")", "\)")

if __name__ == '__main__':
    ROOT_PATH="/home/hoabo/smartpipe/rnavariantcalling/HISAT2"
    for path, subdirs, files in os.walk(ROOT_PATH):
        for subdir in subdirs:
            if(len(subdir)> 1):
                SUB_PATH1 = "%s/%s/HISAT2.Aligned" % (ROOT_PATH, subdir)
                SUB_PATH2 = "%s/%s/HISAT2.Aligned.bam" % (ROOT_PATH, subdir)
                exeCommand(shellEscape(' '.join(["rm", "-f", SUB_PATH1])))
                exeCommand(shellEscape(' '.join(["rm", "-f", SUB_PATH2])))
                print("Done:%s" % (SUB_PATH1))
                print("Done:%s" % (SUB_PATH2))