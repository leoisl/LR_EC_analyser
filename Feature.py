#!/usr/bin/env python
# -*- coding: utf-8 -*-

#TODO: we should add more data here...
class Profile:
    """
    Represents the profile of a Feature for all tools - nb of reads mapped to it, quality of the reads, etc...
    """
    def __init__(self, tools):
        # TODO: we should add more data here...
        self.toolsOrder = tools #will remember the tools' order
        self.tool2NbOfMappedReads = {tool: 0 for tool in tools}

    def computeLargestDiscrepancy(self):
        """
        :return: the largest value between raw_reads and a tool
        """
        nbOfMappedRawReads = self.tool2NbOfMappedReads["raw.bam"]
        return max([abs(nbOfMappedRawReads-self.tool2NbOfMappedReads[tool]) for tool in self.tool2NbOfMappedReads] )

    def getProfileAsArrayForHot(self):
        return [str(self.computeLargestDiscrepancy())] + [str(self.tool2NbOfMappedReads[tool]) for tool in self.toolsOrder]

    def isExpressedInTool(self, tool):
        return self.tool2NbOfMappedReads[tool]>0

    def isExpressedInAnyTool(self):
        return sum(self.tool2NbOfMappedReads.values())>0



class Feature:
    """
    Represents a Feature (gene, transcript, exon, intron, etc)
    """
    def __init__ (self, id, chromosome, begin, end, strand, tools, parent=None):
        self.id=id
        self.chromosome=chromosome
        self.begin=begin
        self.end=end
        self.strand=strand
        self.profile = Profile(tools)
        self.parent=parent

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def getLocusInIGVJSFormat(self):
        return "'%s:%d-%d'"%(self.chromosome, self.begin, self.end)

    def getFeatureLength(self):
        return self.end-self.begin+1

    def getASArrayForHOT(self):
        #TODO: is it better to return just the gene name (this line commented out returns gene name, chr and strand)
        #return ["'%s'"%self.id, self.getLocusInIGVJSFormat(), "'%s'"%self.strand] + self.profile.getProfileAsArrayForHot()
        return ["'%s'"%self.id] + self.profile.getProfileAsArrayForHot()

    def computeRelativeExpression(self, tool):
        """
        Computes the relative expression of this feature in relation to its parent in a given tool
        :return: the relative expression (float)
        """
        parentExpression = self.parent.profile.tool2NbOfMappedReads[tool]
        featureExpression = self.profile.tool2NbOfMappedReads[tool]
        if parentExpression == 0:
            return 0.0
        else:
            return float(featureExpression)/float(parentExpression)



class Gene(Feature):
    """
    Represents a gene - the same as a Feature, but it has several transcripts (features inside it)
    """
    def __init__(self, id, chromosome, begin, end, strand, tools):
        Feature.__init__(self, id, chromosome, begin, end, strand, tools)
        self.transcriptId2Transcript={}

    def getNbOfIsoformsExpressedInTool(self, tool):
        nbOfIsoformsPresent = 0
        for transcript in self.transcriptId2Transcript.values():
            if transcript.profile.tool2NbOfMappedReads[tool]>0:
                nbOfIsoformsPresent+=1
        return nbOfIsoformsPresent