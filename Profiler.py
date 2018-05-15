#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import urllib
from ExternalTools import *

class FeatureProfiler:
    """
    Represents all features and their profiles
    """
    def __init__(self, geneID2gene, tools, tool2Bam, genome, gtf, outputFolder):
        self.geneId2Gene = geneID2gene
        self.tools=tools
        self.tool2Bam = tool2Bam
        self.genome = genome
        self.gtf = gtf
        self.outputFolder = outputFolder

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


    def __buildIGVInfo(self, feature):
        igvInfo = {}
        igvInfo["fastaURL"] = self.genome
        igvInfo["locus"] = feature.getLocusInIGVJSFormat()[1:-1]
        igvInfo["annotationURL"] = self.gtf
        igvInfo["visibilityWindow"] = feature.getFeatureLength()
        for i, tool in enumerate(self.tools):
            igvInfo["tool_%d"%i] = tool

            #single bam version
            igvInfo["bam_%d"%i] = self.tool2Bam[tool]

            '''
            Version where we split the big bam into several bams
            outputFile = self.outputFolder+"/bams/%s_%s.bam"%(tool, feature.id)

            if not self.skip_splitting_bam:
                getBamInCoordinates(self.tool2Bam[tool], igvInfo["locus"], outputFile)

            igvInfo["bam_%d" % i] = "bams/%s_%s.bam"%(tool, feature.id)
            '''




        listOfKeyValues=["%s=%s"%(key,value) for key,value in igvInfo.items()]
        return ["'" + urllib.quote("&".join(listOfKeyValues), safe='') + "'"]

    def __getFeatureASArrayForHOT(self, feature):
        """
        Basically this get the feature as Array for HOT and add IGV info just after the gene in featureAsArrayForHOT
        :param featureAsArrayForHOT:
        :return: a javascript array
        """
        featureAsArrayForHOT = feature.getASArrayForHOT()
        featureAsArrayForHOT = featureAsArrayForHOT [:2] + self.__buildIGVInfo(feature) + featureAsArrayForHOT[2:]
        return "[%s]"%(",".join(featureAsArrayForHOT))

    def geneProfileToJSArrayForHOT(self):
        """
        Transforms this object in a format to JSArray to be put in a Hands-on table
        :return: a string like:
        [ [list with first gene infos], [list with second gene infos] ... ]
        """
        stringList=[]
        for gene in self.geneId2Gene.values():
            if gene.profile.isExpressedInAnyTool():
                stringList.append(self.__getFeatureASArrayForHOT(gene))
        return "[" + ",".join(stringList) + "]"

    def transcriptProfileToJSArrayForHOT(self):
        """
        Transforms this object in a format to JSArray to be put in a Hands-on table
        :return: a string like:
        [ [list with first gene infos], [list with second gene infos] ... ]
        """
        stringList=[]
        for gene in self.geneId2Gene.values():
            if gene.profile.isExpressedInAnyTool():
                for transcript in gene.transcriptId2Transcript.values():
                    if transcript.profile.isExpressedInAnyTool():
                        stringList.append(self.__getFeatureASArrayForHOT(transcript))
        return "[" + ",".join(stringList) + "]"




class StatProfiler:
    def __init__(self, tools, outputFolder):
        self.tools = tools
        self.tool2Stats = {tool: self.parseAlignQCOutput(tool, outputFolder) for tool in tools}
        self.readStatsFeatures = ["TOTAL_READS", "UNALIGNED_READS", "ALIGNED_READS", "MEAN_LENGTH", "SINGLE_ALIGN_READS", "GAPPED_ALIGN_READS", "CHIMERA_ALIGN_READS", "TRANSCHIMERA_ALIGN_READS", "SELFCHIMERA_ALIGN_READS"]
        self.baseStatsFeatures = ["TOTAL_BASES", "UNALIGNED_BASES", "ALIGNED_BASES", "SINGLE_ALIGN_BASES", "GAPPED_ALIGN_BASES", "CHIMERA_ALIGN_BASES", "TRANSCHIMERA_ALIGN_BASES", "SELFCHIMERA_ALIGN_BASES"]
        self.errorStatsFeatures = ["ANY_ERROR", "MISMATCHES", "ANY_DELETION", "ANY_INSERTION", "COMPLETE_DELETION", "HOMOPOLYMER_DELETION", "COMPLETE_INSERTION", "HOMOPOLYMER_INSERTION"]
        self.allFeatures = self.readStatsFeatures + self.baseStatsFeatures + self.errorStatsFeatures

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

    @staticmethod
    def getReadsMeanLength (lengthsFilename):
        lengths=[]
        with gzip.open(lengthsFilename) as file:
            for line in file:
                lineSplit = line.rstrip().split()
                lengths.append(int(lineSplit[4]))
        return float(sum(lengths))/len(lengths)

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

        statsDic["MEAN_LENGTH"] = StatProfiler.getReadsMeanLength(dataFolder+"/lengths.txt.gz")

        return statsDic


    def __toJSArrayForHOT(self, features):
        jsArray=[]
        for feature in features:
            line=["'%s'"%feature]
            for tool in self.tools:
                line.append(str(self.tool2Stats[tool][feature]))
            jsArray.append("[" + ",".join(line) + "]")
        return "[" + ",".join(jsArray) + "]"

    def getReadStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.readStatsFeatures)

    def getBaseStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.baseStatsFeatures)

    def getErrorStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.errorStatsFeatures)
