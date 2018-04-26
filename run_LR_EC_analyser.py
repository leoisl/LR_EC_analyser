#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import sys

"""
Input:
python run_LR_EC_analyser.py -g <genome> -t <transcriptome> raw_reads.fa LoRDEC.fa LoRMA.fa PBCR.fa ...

Data structures:
gene2ECTool2MappingInfo
transcript2ECTool2MappingInfo
MappingInfo is just a class with relevant mapping info to compare 2 EC tools
For now, I think it should have:
    1. # of reads mapped to transcript/gene
    2. overlap size between read and transcript (all the values / average)?

1/ For each .fa file:
    1/ Map it to the genome+transcriptome using GMAP
    2/ Filter the bam to keep only the best hit
    3/ Run AlignQC on this bam
    4/ Process output/data/annotbest.txt and populate gene2ECTool2MappingInfo and transcript2ECTool2MappingInfo
    5/ Create a discrepancy measure (I think now it should be like the largest difference between the # of reads of two tools)
    6/ Create an HTML page sorted by this discrepancy measure and show all the infos we have gathered
"""

"""
best.sorted.gpd is the AlignQC file containing the best mapping of the reads
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

def readBam(bamFilename):
    inputBam = pysam.AlignmentFile(bamFilename, "rb")

    for bamLine in inputBam:

    inputBam.close()

def main():


if __name__=="__main__":
    main()