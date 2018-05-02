#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess
import os


scriptDir = os.path.dirname(os.path.realpath(__file__))

def executeCommandLine(commandLine):
    print "Running %s" % commandLine
    subprocess.check_call(commandLine.split())

def indexGenome(filepath):
    print "Indexing %s..."%filepath
    executeCommandLine("samtools faidx %s" % filepath)
    print "Indexing %s - Done!" % filepath

def gzipFile(file):
    print "Gzipping %s..."%file
    if os.path.exists(file+".gz"):
        os.remove(file+".gz")
    executeCommandLine("gzip -k %s" % file)
    print "Gzipping %s - Done!" % file

def processBam(bam, outputBam, threads):
    print "Sorting and indexing %s..." % bam
    executeCommandLine("samtools sort -@ %d %s -o %s"%(threads, bam, outputBam))
    executeCommandLine("samtools index %s" % (outputBam))
    print "Sorting and indexing %s - Done!" % bam


def runAlignQC(tool, bam, genome, gtf, outputFolder, threads):
    print "Running AlignQC for %s..." % bam
    executeCommandLine("python %s/AlignQC/alignqc/alignqc.py analyze %s -g %s -t %s --output_folder %s/alignqc_out_on_%s --threads %d" % \
                            (scriptDir, bam, genome, gtf, outputFolder, tool, threads))
    print "Running AlignQC for %s - Done!" % bam

def getBamInCoordinates(inputBam, coordinate, outputBam):
    executeCommandLine('samtools view -b %s %s -o %s'%(inputBam, coordinate, outputBam))
    executeCommandLine('samtools index %s' % outputBam)
