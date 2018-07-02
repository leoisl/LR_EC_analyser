#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import urllib
from Category import *

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
        0: read line number (the readline of this annotation in file best.sorted.bed.gz or best.sorted.gpd.gz)
        1: read name
        2: gene name
        3: transcript name
        4: match type (partial or full (if it mapped partially or fully in the genome according to AlignQC))
        5: number of matching exons
        6: highest number of consecutive_exons
        7: number of exons in read
        8: number of exons in reference transcript
        9: number of bp overlapping between read and transcript
        10: read length
        11: transcript length
        12: read range
        13: transcript range
        14: reference line number

        Example:
        5	m150117_043342_42142_c100769800150000001823165407071580_s1_p0/144819/ccs	ENSG00000274276.4	ENST00000624934.3	partial	15	11	16	18	1747	3884	1992	chr21:6446736-6467516	chr21:6445433-6467532	2454


        From seqtools/cli/utilities/gpd_annotate.py:
           1. read line number
		   2. read name
		   3. gene name
		   4. transcript name
		   5. match type
		   6. number of matching exons
		   7. highst number of consecutive_exons
		   8. number of exons in read
		   9. number of exons in reference transcript
		   10. number of bp overlapping
		   11. read lengthread_length
		   12. transcript length
		   13. read range
		   14. transcript range
		   15. reference line number
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
    def __init__(self, tools, outputFolder, start=0, end=4500, step=500):
        self.tools = tools
        self.outputFolder = outputFolder
        self.start = start
        self.end = end
        self.step = step

        self.readStatsFeatures = ["TOTAL_READS", "UNALIGNED_READS", "ALIGNED_READS", "MEAN_LENGTH", "SINGLE_ALIGN_READS", "GAPPED_ALIGN_READS", "CHIMERA_ALIGN_READS", "TRANSCHIMERA_ALIGN_READS", "SELFCHIMERA_ALIGN_READS"]
        self.baseStatsFeatures = ["TOTAL_BASES", "UNALIGNED_BASES", "ALIGNED_BASES", "SINGLE_ALIGN_BASES", "GAPPED_ALIGN_BASES", "CHIMERA_ALIGN_BASES", "TRANSCHIMERA_ALIGN_BASES", "SELFCHIMERA_ALIGN_BASES"]
        self.errorStatsFeatures = ["ANY_ERROR", "MISMATCHES", "ANY_DELETION", "ANY_INSERTION", "COMPLETE_DELETION", "HOMOPOLYMER_DELETION", "COMPLETE_INSERTION", "HOMOPOLYMER_INSERTION"]
        self.annotationStatsFeatures = ["GENES_DETECTED_ANY_MATCH", "GENES_DETECTED_FULL_MATCH", "TRANSCRIPTS_DETECTED_ANY_MATCH", "TRANSCRIPTS_DETECTED_FULL_MATCH"]
        self.allFeatures = self.readStatsFeatures + self.baseStatsFeatures + self.errorStatsFeatures + self.annotationStatsFeatures

        # this stores as key a feature that should be shown as % instead of raw numbers, and which nb to use as %
        self.feature2UpperLimitToBeUsedInPercentage={
            "UNALIGNED_READS": "TOTAL_READS",
            "ALIGNED_READS": "TOTAL_READS",
            "SINGLE_ALIGN_READS": "ALIGNED_READS",
            "GAPPED_ALIGN_READS": "ALIGNED_READS",
            "CHIMERA_ALIGN_READS": "ALIGNED_READS",
            "TRANSCHIMERA_ALIGN_READS": "ALIGNED_READS",
            "SELFCHIMERA_ALIGN_READS": "ALIGNED_READS",
            "UNALIGNED_BASES": "TOTAL_BASES",
            "ALIGNED_BASES": "TOTAL_BASES",
            "SINGLE_ALIGN_BASES": "ALIGNED_BASES",
            "GAPPED_ALIGN_BASES": "ALIGNED_BASES",
            "CHIMERA_ALIGN_BASES": "ALIGNED_BASES",
            "TRANSCHIMERA_ALIGN_BASES": "ALIGNED_BASES",
            "SELFCHIMERA_ALIGN_BASES": "ALIGNED_BASES",
            "ANY_ERROR": "ALIGNMENT_BASES",
            "MISMATCHES": "ALIGNMENT_BASES",
            "ANY_DELETION": "ALIGNMENT_BASES",
            "ANY_INSERTION": "ALIGNMENT_BASES",
            "COMPLETE_DELETION": "ALIGNMENT_BASES",
            "HOMOPOLYMER_DELETION": "ALIGNMENT_BASES",
            "COMPLETE_INSERTION": "ALIGNMENT_BASES",
            "HOMOPOLYMER_INSERTION": "ALIGNMENT_BASES"
        }

        self.feature2Translations={
            "ANY_ERROR": "ERROR_RATE",
            "ANY_DELETION": "DELETION",
            "ANY_INSERTION": "INSERTION",
            "COMPLETE_DELETION": "NON_HOMOPOLYMER_DELETION",
            "COMPLETE_INSERTION": "NON_HOMOPOLYMER_INSERTION"
        }

        self.parseAlignQCOutputForAllTools()


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


    def __processLengthsFile (self, tool, dataFolder):
        '''
        Fills:
            self.tool2Stats[tool]["MEAN_LENGTH"]
            self.tool2Stats[tool]["ALIGNED_SIZE_BINS"]
            self.tool2Stats[tool]["UNALIGNED_SIZE_BINS"]


        /data/lengths.txt.gz file:
            -nb of lines = nb of reads
            -About column 1:
                Aligned reads (the sum of):
                    Single-align reads = line[1]=="original"
                    Gapped-align reads = line[1]=="gapped"
                    Chimeric reads (sum of):
                        Self-chimera: line[1]=="self-chimera" or line[1]=="self-chimera-atypical"
                        Trans-chimera: line[1]=="chimera"
                Unaligned reads:
                    line[1]=="unaligned"
            -Column 4 is read length
        '''
        alignedSizeBins = NumberCategory(self.start, self.end, self.step)
        unalignedSizeBins = NumberCategory(self.start, self.end, self.step)
        totalSizeBins = NumberCategory(self.start, self.end, self.step)
        alignedReadsClassifications=["original", "gapped", "self-chimera", "self-chimera-atypical", "chimera"]
        unalignedReadsClassifications = ["unaligned"]


        #process the lengths file
        lengths=[]
        with gzip.open(dataFolder + "/lengths.txt.gz") as file:
            for line in file:
                lineSplit = line.rstrip().split()

                #get the fields we are interested
                classification=lineSplit[1]
                length=int(lineSplit[4])

                #add the read to the appropriate size bin
                if classification in alignedReadsClassifications:
                    alignedSizeBins.addDataPoint(length, None)
                elif classification in unalignedReadsClassifications:
                    unalignedSizeBins.addDataPoint(length, None)
                else:
                    raise Exception("Unknown classification %s in __processLengthsFile()" % classification)
                totalSizeBins.addDataPoint(length, None)


                #add length to lengths
                lengths.append(length)

        #fills the object
        self.tool2Stats[tool]["MEAN_LENGTH"] = float(sum(lengths))/len(lengths)
        self.tool2Stats[tool]["TOTAL_SIZE_BINS"] = totalSizeBins
        self.tool2Stats[tool]["ALIGNED_SIZE_BINS"] = alignedSizeBins
        self.tool2Stats[tool]["UNALIGNED_SIZE_BINS"] = unalignedSizeBins



    def __parseAlignQCOutput(self, tool):
        dataFolder = self.outputFolder+"/alignqc_out_on_%s/data" % tool
        self.tool2Stats[tool] = {}


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
        self.tool2Stats[tool].update(StatProfiler.readFileComposedOfPairStringIntToDict(dataFolder + "/alignment_stats.txt"))

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
        self.tool2Stats[tool].update(StatProfiler.readFileComposedOfPairStringIntToDict(dataFolder + "/error_stats.txt"))

        self.tool2Stats[tool]["GENES_DETECTED_ANY_MATCH"] = StatProfiler.processRarefractionFile(dataFolder + "/gene_rarefraction.txt",
                                                                                                 self.tool2Stats[tool]["TOTAL_READS"])
        self.tool2Stats[tool]["GENES_DETECTED_FULL_MATCH"] = StatProfiler.processRarefractionFile(dataFolder + "/gene_full_rarefraction.txt",
                                                                                                  self.tool2Stats[tool]["TOTAL_READS"])
        self.tool2Stats[tool]["TRANSCRIPTS_DETECTED_ANY_MATCH"] = StatProfiler.processRarefractionFile(dataFolder + "/transcript_rarefraction.txt",
                                                                                                 self.tool2Stats[tool]["TOTAL_READS"])
        self.tool2Stats[tool]["TRANSCRIPTS_DETECTED_FULL_MATCH"] = StatProfiler.processRarefractionFile(dataFolder + "/transcript_full_rarefraction.txt",
                                                                                                  self.tool2Stats[tool]["TOTAL_READS"])
        '''
        Added the following fields:
        self.tool2Stats[tool]["MEAN_LENGTH"]
        self.tool2Stats[tool]["ALIGNED_SIZE_BINS"]
        self.tool2Stats[tool]["UNALIGNED_SIZE_BINS"]
        '''
        self.__processLengthsFile(tool, dataFolder)



    def parseAlignQCOutputForAllTools(self):
        self.tool2Stats = {}
        #populate self.tool2Stats
        for tool in self.tools:
            self.__parseAlignQCOutput(tool)

    def getStatsForToolAndMetric(self, tool, metric):
        if metric not in self.feature2UpperLimitToBeUsedInPercentage:
            return self.tool2Stats[tool][metric]
        else:
            return float(self.tool2Stats[tool][metric]) / float(self.tool2Stats[tool][self.feature2UpperLimitToBeUsedInPercentage[metric]]) * 100


    def getFeatureName(self, feature):
        """
        Returns a nicer name for the feature
        """
        return feature if feature not in self.feature2Translations else self.feature2Translations[feature]


    def getNiceDescriptionForFeature(self, feature):
        if feature not in self.feature2UpperLimitToBeUsedInPercentage:
            return self.getFeatureName(feature)
        else:
            return "%s in %% over %s"%(self.getFeatureName(feature), self.getFeatureName(self.feature2UpperLimitToBeUsedInPercentage[feature]))

    def __toJSArrayForHOT(self, features):
        jsArray=[]
        for feature in features:
            if feature not in self.feature2UpperLimitToBeUsedInPercentage:
                line=["'%s'" % self.getFeatureName(feature)]
            else:
                line = ["'%s (%%)'" % self.getFeatureName(feature)]

            for tool in self.tools:
                if feature not in self.feature2UpperLimitToBeUsedInPercentage:
                    line.append(str(self.tool2Stats[tool][feature]))
                else:
                    line.append("%.3f" % (float(self.tool2Stats[tool][feature]) / float(self.tool2Stats[tool][self.feature2UpperLimitToBeUsedInPercentage[feature]]) * 100))
            jsArray.append("[" + ",".join(line) + "]")
        return "[" + ",".join(jsArray) + "]"

    def getReadStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.readStatsFeatures)

    def getBaseStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.baseStatsFeatures)

    def getErrorStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.errorStatsFeatures)

    def getAnnotationStatsAsJSArrayForHOT(self):
        return self.__toJSArrayForHOT(self.annotationStatsFeatures)
