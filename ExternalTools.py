#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess
import os

def indexGenome(filepath):
    print "Indexing %s..."%filepath
    commandLine = "samtools faidx %s" % filepath
    print "Running %s"%commandLine
    subprocess.check_call(commandLine.split())
    print "Indexing %s - Done!" % filepath

def gzipFile(file):
    print "Gzipping %s..."%file
    if os.path.exists(file+".gz"):
        os.remove(file+".gz")
    commandLine = "gzip -k %s" % file
    print "Running %s"%commandLine
    subprocess.check_call(commandLine.split())
    print "Gzipping %s - Done!" % file

def processBam(bam, outputBam, threads):
    print "Sorting and indexing %s..." % bam
    scriptDir = os.path.dirname(os.path.realpath(__file__))
    commandLine = "bash %s/processBam.sh %s %s %d"%(scriptDir, bam, outputBam, threads)
    print "Running %s"%commandLine
    subprocess.check_call(commandLine.split())
    print "Sorting and indexing %s - Done!" % bam


def runAlignQC(tool, bam, genome, gtf, outputFolder, threads):
    print "Running AlignQC for %s..." % bam


    commandLine = "alignqc analyze %s -g %s -t %s --output_folder %s/alignqc_out_on_%s --threads %d" % \
                            (bam, genome, gtf, outputFolder, tool, threads)
    print "Running %s" % commandLine
    subprocess.check_call(commandLine.split())
    print "Running AlignQC for %s - Done!" % bam
