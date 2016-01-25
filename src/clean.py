__author__ = 'linuxpham'
import os

if __name__ == '__main__':
    ROOT_PATH="/home/hoabo/smartpipe/rnavariantcalling/HISAT2"
    for path, subdirs, files in os.walk(ROOT_PATH):
        for subdir in subdirs:
            print(subdir)