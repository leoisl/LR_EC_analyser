#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Feature import *

def parseGTFToGetGenes(gtf, tools):
    print "Parsing %s..."%gtf
    geneID2gene={}
    with open(gtf) as gtfFile:
        for line in gtfFile:
            if line[0]!="#":
                lineSplit = line.rstrip().split()
                feature = lineSplit[2]
                chromosome = lineSplit[0]
                begin = int(lineSplit[3])
                end = int(lineSplit[4])
                strand = lineSplit[6]

                if feature=="gene":
                    geneId = lineSplit[lineSplit.index("gene_id") + 1][1:-2]
                    geneID2gene[geneId] = Gene(geneId, chromosome, begin, end, strand, tools)
                elif feature=="transcript":
                    transcriptId = lineSplit[lineSplit.index("transcript_id") + 1][1:-2]
                    geneID2gene[geneId].transcriptId2Transcript[transcriptId]=Feature(transcriptId, chromosome, begin, end, strand, tools)
    print "Parsing %s - Done!" % gtf
    return geneID2gene