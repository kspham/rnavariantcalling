#!/usr/bin/env python

import subprocess
import os
import sys
import numpy as np

from tempfile import *
from threading import Thread
from Queue import Queue

BAYES, REF, REG, NUM, BAM, OUT = sys.argv[1:]

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

oLogger = LogInstance("variant_log.txt")

def exeCommand(sCommand):
    process = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    errcode = process.returncode
    sys.stdout.write(out)
    if err:
        oLogger.debug("ERROR: " + err)
        oLogger.debug("COMMAND: " + sCommand)

def shellEscape(s):
    return s.replace("(", "\(").replace(")", "\)")

class Worker(Thread):
    """Thread executing tasks from a given tasks queue"""

    def __init__(self, tasks):
        Thread.__init__(self)
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.tasks.get()
            try: func(*args, **kargs)
            except Exception, e: print e
            self.tasks.task_done()

class ThreadPool:
    """Pool of threads consuming tasks from a queue"""

    def __init__(self, num_threads):
        self.tasks = Queue(num_threads)
        for _ in range(num_threads): Worker(self.tasks)

    def add_task(self, func, *args, **kargs):
        """Add a task to the queue"""
        self.tasks.put((func, args, kargs))

    def wait_completion(self):
        """Wait for completion of all the tasks in the queue"""
        self.tasks.join()

def call(f_tmp):
    exeCommand(shellEscape(' '.join([BAYES, "-f", REF, "-C 2","-t", f_tmp ,BAM, ">", f_tmp+".vcf"])))
    oLogger.debug("DONE: " + f_tmp.name)


if __name__ == '__main__':
    # 1) Init a Thread pool with the desired number of threads
    pool = ThreadPool(int(THREADS))

    regions = np.array(open(REG).read().strip().split("\n"), dtype=None)
    regions = regions[np.random.permutation(len(regions))]

    temps = []
    # 2) Add the task to the queue
    for i in range(0, len(regions), int(NUM)):
        f_tmp = NamedTemporaryFile(delete=False)
        f_tmp.write("\n".join(regions[i:i+int(NUM)]))
        f_tmp.close()
        temps.append(f_tmp.name)
        pool.add_task(call,f_tmp)
    f_tmp = NamedTemporaryFile(delete=False)
    f_tmp.write("\n".join(regions[i:]))
    f_tmp.close()
    temps.append(f_tmp.name)
    pool.add_task(call,f_tmp)

    # 3) Wait for completion
    pool.wait_completion()

    f_tmps = open("tmp_vcfs.txt", "w")
    f_tmps.write("\n".join([f+".vcf" for f in temps]))
    exeCommand(shellEscape(' '.join(["vcf-concat", "-f", f_tmps.name, ">", OUT])))
    for f in temps:
        os.unlink(f.name)
        os.unlink(f.name+".vcf")
