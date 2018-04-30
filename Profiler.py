#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import urllib

class FeatureProfiler:
    """
    Represents all features and their profiles
    """
    def __init__(self, geneID2gene, tools, tool2Bam):
        self.geneId2Gene = geneID2gene
        self.tools=tools
        self.tool2Bam = tool2Bam

    def populateFromAnnotbest(self, tool, outputFolder):
        """
        Reads the annotbest.txt.gz file from AlignQC and populates this profiler
        annotbest.txt is the AlignQC file that contains the best mapping of the ANNOTATED READS (reads that we were able to align and that AlignQC was able to assign to a transcript)

        annotbest.txt format:
        columns:
        0: id of the line (?)
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

        Example:
        5	m150117_043342_42142_c100769800150000001823165407071580_s1_p0/144819/ccs	ENSG00000274276.4	ENST00000624934.3	partial	15	11	16	18	1747	3884	1992	chr21:6446736-6467516	chr21:6445433-6467532	2454
        """
        dataFolder = outputFolder+"/alignqc_out_on_%s/data"%tool
        with gzip.open(dataFolder+"/annotbest.txt.gz") as file:
            for line in file:
                lineSplit = line.rstrip().split()
                geneId = lineSplit[2]
                transcriptId = lineSplit[3]

                self.geneId2Gene[geneId].profile.tool2NbOfMappedReads[tool] += 1
                self.geneId2Gene[geneId].transcriptId2Transcript[transcriptId].profile.tool2NbOfMappedReads[tool] += 1


    def buildIGVInfo(self, gene, genome, gtf):
        igvInfo = {}
        igvInfo["fastaURL"] = genome
        igvInfo["locus"] = gene.getLocusInIGVJSFormat()[1:-1]
        igvInfo["annotationURL"] = gtf
        for i, tool in enumerate(self.tools):
            igvInfo["tool_%d"%i]=tool
            igvInfo["bam_%d"%i]=self.tool2Bam[tool]

        listOfKeyValues=["%s=%s"%(key,value) for key,value in igvInfo.items()]
        return ["'" + urllib.quote("&".join(listOfKeyValues), safe='') + "'"]

    def toJSArrayForHOT(self, genome, gtf):
        """
        Transforms this object in a format to JSArray to be put in a Hands-on table
        :return: a string like:
        [ [list with first gene infos], [list with second gene infos] ... ]
        """
        stringList=[]
        for gene in self.geneId2Gene:
            if self.geneId2Gene[gene].profile.isExpressedInAnyTool():
                stringList.append("[" + ",".join(self.geneId2Gene[gene].getASArrayForHOT() +
                                                 self.buildIGVInfo(self.geneId2Gene[gene], genome, gtf)) + "]")
        return "[" + ",".join(stringList) + "]"


class StatProfiler:
    @staticmethod
    def readFileComposedOfPairStringIntToDict(filename):
        stats = {}
        with open(filename) as file:
            for line in file:
                lineSplit = line.split()
                stats[lineSplit[0]] = int(lineSplit[1])
        return stats

    @staticmethod
    def processRarefractionFile(filename, totalReads):
        """
        Goes through the rarefraction file and get the values for rarefraction in the last line
        :param filename:
        :return:
        """
        with open(filename) as file:
            lastLineSplitted = file.readlines()[-1].split()

        nbOfReads = int(lastLineSplitted[0])
        rarefractionLow = int(float(lastLineSplitted[1]))
        rarefractionMed = int(float(lastLineSplitted[2]))
        rarefractionHigh = int(float(lastLineSplitted[3]))

        if rarefractionLow == rarefractionMed and rarefractionLow == rarefractionHigh and nbOfReads == totalReads:
            return rarefractionLow
        else:
            raise Exception("Rarefraction processing error!")

    def parseAlignQCOutput(self, tool, outputFolder):
        dataFolder = outputFolder+"/alignqc_out_on_%s/data" % tool
        statsDic = {}

        '''
        Added the following fields:
        TOTAL_READS	1282
        UNALIGNED_READS	200
        ALIGNED_READS	1082
        SINGLE_ALIGN_READS	996
        GAPPED_ALIGN_READS	73
        CHIMERA_ALIGN_READS	13
        TRANSCHIMERA_ALIGN_READS	0
        SELFCHIMERA_ALIGN_READS	13
        TOTAL_BASES	4561236
        UNALIGNED_BASES	2257927
        ALIGNED_BASES	2303309
        SINGLE_ALIGN_BASES	2222158
        GAPPED_ALIGN_BASES	58561
        CHIMERA_ALIGN_BASES	22590
        TRANSCHIMERA_ALIGN_BASES	0
        SELFCHIMERA_ALIGN_BASES	22590
        '''
        statsDic.update(StatProfiler.readFileComposedOfPairStringIntToDict(dataFolder + "/alignment_stats.txt"))

        '''
        Added the following fields:
        ALIGNMENT_COUNT	501
        ALIGNMENT_BASES	1097266
        ANY_ERROR	118039
        MISMATCHES	40866
        ANY_DELETION	28083
        COMPLETE_DELETION	13374
        HOMOPOLYMER_DELETION	14709
        ANY_INSERTION	49090
        COMPLETE_INSERTION	28668
        HOMOPOLYMER_INSERTION	20422
        '''
        statsDic.update(StatProfiler.readFileComposedOfPairStringIntToDict(dataFolder + "/error_stats.txt"))

        statsDic["GENES_DETECTED_ANY_MATCH"] = StatProfiler.processRarefractionFile(dataFolder + "/gene_rarefraction.txt",
                                                                       statsDic["TOTAL_READS"])
        statsDic["GENES_DETECTED_FULL_MATCH"] = StatProfiler.processRarefractionFile(dataFolder + "/gene_full_rarefraction.txt",
                                                                        statsDic["TOTAL_READS"])

        return statsDic

    def __init__(self, tools, outputFolder):
        self.tools = tools
        self.tool2Stats = {tool: self.parseAlignQCOutput(tool, outputFolder) for tool in tools}

    def toJSArrayForHOT(self):
        jsArray=[]
        for feature in self.tool2Stats[self.tools[0]]:
            line=["'%s'"%feature]
            for tool in self.tools:
                line.append(str(self.tool2Stats[tool][feature]))
            jsArray.append("[" + ",".join(line) + "]")
        return "[" + ",".join(jsArray) + "]"