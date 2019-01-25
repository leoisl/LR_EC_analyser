#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import urllib
import Utils
import copy
from Category import *

class FeatureProfiler:
    """
    Represents all features and their profiles
    Used in the Gene read alignment viewer and Transcript read alignment viewer
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
        0: read line number (the readline of this annotation in file best.sorted.bed.gz or best.sorted.gpd.gz - I guess it is in increasing consecutive order)
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

        self.readStatsFeatures = ["TOTAL_READS", "MEAN_LENGTH", "ALIGNED_READS", "UNALIGNED_READS", "SINGLE_ALIGN_READS", "GAPPED_ALIGN_READS", "CHIMERA_ALIGN_READS", "TRANSCHIMERA_ALIGN_READS", "SELFCHIMERA_ALIGN_READS"]
        self.baseStatsFeatures = ["TOTAL_BASES", "ALIGNED_BASES", "UNALIGNED_BASES", "SINGLE_ALIGN_BASES", "GAPPED_ALIGN_BASES", "CHIMERA_ALIGN_BASES", "TRANSCHIMERA_ALIGN_BASES", "SELFCHIMERA_ALIGN_BASES"]
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
        self.tool2Stats[tool].update(Utils.FileUtils.readFileComposedOfPairStringIntToDict(dataFolder + "/alignment_stats.txt"))

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
        self.tool2Stats[tool].update(Utils.FileUtils.readFileComposedOfPairStringIntToDict(dataFolder + "/error_stats.txt"))

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


class SplicingSitesProfiler:
    """
    Represents the splicing sites data read from <alignqc_folder/data/junvar.txt>
    """
    def __init__(self, tools, outputFolder):
        self.__tools=tools
        self.__outputFolder = outputFolder

        #represents a map tool -> distance to the nearest SS -> count
        self.__tool2SpliceSiteDistance2Count = {tool : Utils.FileUtils.readFileComposedOfPairIntsToDict(outputFolder + "/alignqc_out_on_%s/data/junvar.txt" % tool) for tool in tools}

        '''
        TODO: discover why AlignQC provides only from -39 to +39 distance - what does it do with longer distances? Understand the code and the algorithm
        TODO: here we fix the dictionary so it will always contains from -39 to +39
        '''
        for tool in tools:
            for i in xrange(-39, 40):
                if i not in self.__tool2SpliceSiteDistance2Count[tool]:
                    self.__tool2SpliceSiteDistance2Count[tool][i]=0

        #represent the total number of splice sites per tool
        self.__tool2NbOfSpliceSites = {tool: sum(self.__tool2SpliceSiteDistance2Count[tool].values()) for tool in tools}

    def __fixForPercentage(self, tool2Values):
        """
        tool2Values: a dict of tool 2 values
        Divides all values by self.__tool2NbOfSpliceSites[tool]
        :return: Fixed tool2Values
        """
        for tool in tool2Values:
            tool2Values[tool] = float(tool2Values[tool])/float(self.__tool2NbOfSpliceSites[tool])*100
        return tool2Values

    def __getTool2NbOfSpliceSites(self, condition, inPercentage):
        """
        :param condition: the condition that the splice site must satisfy
        :param inPercentage: %-wise?
        :return: tool2NbOfSpliceSites
        """
        #build tool2NbOfSpliceSites
        tool2NbOfSpliceSites={}
        for tool in self.__tools:
            SpliceSiteDistance2Count = self.__tool2SpliceSiteDistance2Count[tool]
            tool2NbOfSpliceSites[tool] = sum([SpliceSiteDistance2Count[ssDistance] for ssDistance in SpliceSiteDistance2Count if condition(ssDistance)])

        # fix for the percentage
        if inPercentage:
            tool2NbOfSpliceSites = self.__fixForPercentage(tool2NbOfSpliceSites)

        return tool2NbOfSpliceSites

    def getTool2NbOfCorrectSpliceSites(self, inPercentage=False):
        """
        Correct SS is where the distance is 0
        :param inPercentage: %-wise?
        """
        def correctSSCondition(ssDistance):
            return ssDistance == 0
        return self.__getTool2NbOfSpliceSites(correctSSCondition, inPercentage)

    def getTool2NbOfIncorrectSpliceSites(self, inPercentage=False):
        """
        Incorrect SS is where the distance is != 0
        :param inPercentage: %-wise?
        """
        def incorrectSSCondition(ssDistance):
            return ssDistance != 0
        return self.__getTool2NbOfSpliceSites(incorrectSSCondition, inPercentage)

    def getTool2NbOfIncorrectNearSpliceSites(self, inPercentage=False):
        """
        Incorrect Near SS is where the distance is != 0, but at most 2
        :param inPercentage: %-wise?
        """
        def incorrectNearSSCondition(ssDistance):
            return ssDistance != 0 and ssDistance >= -2 and ssDistance <= 2
        return self.__getTool2NbOfSpliceSites(incorrectNearSSCondition, inPercentage)


    def getTool2NbOfIncorrectMultipleOf3SpliceSites(self, inPercentage=False):
        """
        Incorrect Multiple Of 3 SS is where the distance is != 0, but multiple of 3
        :param inPercentage: %-wise?
        """
        def incorrectMultipleOf3SSCondition(ssDistance):
            return ssDistance != 0 and ssDistance%3 == 0
        return self.__getTool2NbOfSpliceSites(incorrectMultipleOf3SSCondition, inPercentage)

    def getTool2NbOfIncorrectFarSpliceSites(self, inPercentage=False):
        """
        Incorrect Far SS is where the distance is != 0, larger than 2, and not multiple of 3
        :param inPercentage: %-wise?
        """
        def incorrectFarSSCondition(ssDistance):
            return (ssDistance < -2 or ssDistance > 2) and ssDistance%3 != 0
        return self.__getTool2NbOfSpliceSites(incorrectFarSSCondition, inPercentage)

    def getTool2SpliceSiteDistance2Count(self, inPercentage=False):
        """
        Returns self.__tool2SpliceSiteDistance2Count
        :param inPercentage: %-wise?
        :return: self.__tool2SpliceSiteDistance2Count
        """
        if not inPercentage:
            tool2SpliceSiteDistance2CountCopy = copy.deepcopy(self.__tool2SpliceSiteDistance2Count)
            #iterating the original, not the copy, because of modifications in parallel to iteration
            for tool in self.__tool2SpliceSiteDistance2Count:
                for ssDistance, count in self.__tool2SpliceSiteDistance2Count[tool].iteritems():
                    if ssDistance == 0:
                        del tool2SpliceSiteDistance2CountCopy[tool][ssDistance]
        else:
            tool2SpliceSiteDistance2CountCopy = copy.deepcopy(self.__tool2SpliceSiteDistance2Count)
            # iterating the original, not the copy, because of modifications in parallel to iteration
            for tool in self.__tool2SpliceSiteDistance2Count:
                for ssDistance, count in self.__tool2SpliceSiteDistance2Count[tool].iteritems():
                    if ssDistance == 0:
                        del tool2SpliceSiteDistance2CountCopy[tool][ssDistance]
                    else:
                        tool2SpliceSiteDistance2CountCopy[tool][ssDistance] = \
                            float(tool2SpliceSiteDistance2CountCopy[tool][ssDistance])/float(self.__tool2NbOfSpliceSites[tool])*100
        return tool2SpliceSiteDistance2CountCopy


class ReadSetProfiler:
    """
    Represents the read set for each tool with some additional characteristics obtained from annotbest.txt.gz
    """
    def __init__(self, tools, outputFolder):
        self.tools = tools

        #create the matchType2Tool2Count structure
        #only for multi-exonic transcripts
        self.tool2MatchType={tool: TextCategory(["full", "partial"]) for tool in tools}
        # create the tool2NbOfMatchingExons and tool2HighestNbOfConsecutiveExons structures
        self.tool2NbOfMatchingExons = {tool:NumberCategory(1, 16, 1) for tool in tools}
        self.tool2NbOfMatchingExonsCompacted = {tool: NumberCategory(1, 10, 4) for tool in tools}
        self.tool2HighestNbOfConsecutiveExons = {tool: NumberCategory(1, 16, 1) for tool in tools}

        #populate
        for tool in tools:
            self.__populateFromAnnotbest(tool, outputFolder)


    def __populateFromAnnotbest(self, tool, outputFolder):
        """
        Reads the annotbest.txt.gz file from AlignQC and populates this profiler
        annotbest.txt is the AlignQC file that contains the best mapping of the ANNOTATED READS (reads that we were able to align and that AlignQC was able to assign to a transcript)

        annotbest.txt format (the only interesting columns, starting at 0):
        4: match type ("partial" or "full" (if it mapped partially or fully to the best transcript, according to AlignQC))
            -use it to get # and % of reads that cover a transcript fully or partially;
            -consider only *multi-exonic transcripts*
                -you can get which transcripts are multi-exonic by looking at 8: number of exons in reference transcript
            -single-exonic transcript seems hard to evaluate if a match is full or partial (it is always partial - example on the raw dataset):
                leandro@ngs-provisoire:/data2/ASTER/error_correction/LR_EC_analyser/temp$ awk '$9==1{print $5}' annotbest.txt | sort | uniq -c
                24008 partial
            -for multi-exonic, it seems far better:
				leandro@ngs-provisoire:/data2/ASTER/error_correction/LR_EC_analyser/temp$ awk '$9>1{print $5}' annotbest.txt | sort | uniq -c
                121095 full
                309595 partial
        5: number of matching exons
            -# of exons matching to the best transcript (the best transcript is the reference transcript from the best alignment)
            -use for stats - # of identified exons as a distribution (also as %, according to the total # of reads)
        6: highest number of consecutive_exons
            -use for stats - long reads keep or destroy exon connectivity? - distribution of highest number of consecutive exons
        8: number of exons in reference transcript
            -use this to identify muti-exonic transcripts


        Worth a mention:
        7: number of exons in read
            -DO NOT USE THIS
            -# of identified exons in a read, not necessarily stemming from the best transcript mapping
            -TODO: check if it is really this with a couple of examples


        Example:
        5	m150117_043342_42142_c100769800150000001823165407071580_s1_p0/144819/ccs	ENSG00000274276.4	ENST00000624934.3	partial	15	11	16	18	1747	3884	1992	chr21:6446736-6467516	chr21:6445433-6467532	2454
        """
        dataFolder = outputFolder+"/alignqc_out_on_%s/data"%tool
        with gzip.open(dataFolder+"/annotbest.txt.gz") as file:
            for line in file:
                lineSplit = line.rstrip().split()
                matchType = lineSplit[4]
                nbOfMatchingExons = int(lineSplit[5])
                highestNbOfConsecutiveExons = int(lineSplit[6])
                nbOfExonsInRefTranscript = int(lineSplit[8])

                # only for multi-exonic transcripts
                if nbOfExonsInRefTranscript >= 2:
                    self.tool2MatchType[tool].addDataPointAndIAlreadyKnowTheCategory(matchType)

                #all transcripts
                self.tool2NbOfMatchingExons[tool].addDataPoint(nbOfMatchingExons)
                self.tool2NbOfMatchingExonsCompacted[tool].addDataPoint(nbOfMatchingExons)
                self.tool2HighestNbOfConsecutiveExons[tool].addDataPoint(highestNbOfConsecutiveExons)