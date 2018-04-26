#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess
import os
import gzip


"""
Input:
python run_LR_EC_analyser.py -g <genome> -gtf <transcriptome> raw_reads.bam LoRDEC.bam LoRMA.bam PBCR.bam ...

Data structures (first version gene only):
gene2ECTool2MappingInfo
gene2coordinates

MappingInfo is just a class with relevant mapping info to compare 2 EC tools
For now, I think it should have:
    1. # of reads mapped to transcript/gene
    2. overlap size between read and transcript (all the values / average)?
    

0/ Process the gtf and get the gene coordinates
    
1/ For each .bam file:
    1/ Filter the bam to keep only the best hit
    2/ Run AlignQC on this bam
    3/ Process output/data/annotbest.txt and populate gene2ECTool2MappingInfo and transcript2ECTool2MappingInfo
    4/ Create a discrepancy measure (I think now it should be like the largest difference between the # of reads in the raw_reads and one of the tools)
    5/ Create an HTML page sorted by this discrepancy measure and show all the infos we have gathered
    6/ On clicking in one of the genes, we give to igv viewer the genome, transcriptome, all bams, and the gene coordinate



To view results, execute
python run_LR_EC_analyser --view <path_to_results>
    This will start a ftp server on <path_to_results>
    Then will open a browser to view <path_to_results>

"""


"""
annotbest.txt is the AlignQC file that contains the best mapping of the ANNOTATED READS (reads that we were able to align and that AlignQC was able to assign to a transcript) 
annotbest.txt format:
columns:
0: id of the line
1: read name
2: gene name
3: transcript name
4: partial or full (if it mapped partially or fully in the genome according to AlignQC)
5: # ?
6: # most consecutive exons in read
7: # exons in read
8: # exons in the transcript
9: overlap size between read and transcript
10: read length
11: transcript length
12: alignment coordinates
13: transcript coordinates
5	m150117_043342_42142_c100769800150000001823165407071580_s1_p0/144819/ccs	ENSG00000274276.4	ENST00000624934.3	partial	15	11	16	18	1747	3884	1992	chr21:6446736-6467516	chr21:6445433-6467532	2454
"""

scriptDir = os.path.dirname(os.path.realpath(__file__))

class Coordinate:
    def __init__ (self, chromosome, begin, end, strand):
        self.chromosome=chromosome
        self.begin=begin
        self.end=end
        self.strand=strand

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)


def parseGTFToGetGenes(gtf):
    print "Parsing %s..."%gtf
    genes=[]
    with gzip.open(gtf) as gtfFile:
        for line in gtfFile:
            lineSplit = line.rstrip().split()
            if lineSplit[2]=="gene":
                genes.append(Coordinate(lineSplit[0], int(lineSplit[3]), int(lineSplit[4]), lineSplit[6]))
    print "Parsing %s - Done!" % gtf
    return genes

def getOnlyBestHitsFromBam(bam):
    print "Getting only the best hits from %s..." % bam
    commandLine = "bash %s/getOnlyBestHitsFromBam.sh %s"%(scriptDir, bam)
    print "Running %s"%commandLine
    subprocess.check_call(commandLine.split())
    print "Getting only the best hits from %s - Done!" % bam

def runAlignQC(bam, genome, gtf, threads):
    print "Running AlignQC for %s..." % bam
    commandLine = "alignqc analyze %s -g %s -t %s --output_folder %s --threads %d" % \
                            (bam, genome, gtf, "alignqc_out_on_"+os.path.basename(bam), threads)
    print "Running %s" % commandLine
    subprocess.check_call(commandLine.split())
    print "Running AlignQC for %s - Done!" % bam

def main():
    parser = argparse.ArgumentParser(description='Long read error corrector analyser.')
    parser.add_argument('bams', metavar='<file.bam>', type=str, nargs='+',
                        help='BAM files of the Fastas output by the correctors')
    parser.add_argument("--genome", dest="genome", help="The genome in gzipped fasta file (this file needs to be a .gz)")
    parser.add_argument("--gtf", dest="gtf", help="The transcriptome as gzipped GTF file (this file needs to be a .gz)")
    parser.add_argument("--raw", dest="rawBam", help="The BAM file of the raw reads (i.e. the uncorrected long reads file)")
    parser.add_argument("-t", dest="threads", type=int, help="Number of threads to use")
    args=parser.parse_args()

    #get genes as coordinates
    genes = parseGTFToGetGenes(args.gtf)

    #get only the best hits from the BAM file
    for bam in args.bams:
        getOnlyBestHitsFromBam(bam)

    #update the bam files
    bestHitsBams = [bamfile+".best_hits.bam" for bamfile in args.bams]

    #Run AlignQC on the bestHitsBams
    for bam in bestHitsBams:
        runAlignQC(bam, args.genome, args.gtf, args.threads)






if __name__=="__main__":
    main()