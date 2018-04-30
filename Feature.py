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
        return [str(self.tool2NbOfMappedReads[tool]) for tool in self.toolsOrder] + [str(self.computeLargestDiscrepancy())]

    def isExpressedInAnyTool(self):
        return sum(self.tool2NbOfMappedReads.values())>0



class Feature:
    """
    Represents a Feature (gene, transcript, exon, intron, etc)
    """
    def __init__ (self, id, chromosome, begin, end, strand, tools):
        self.id=id
        self.chromosome=chromosome
        self.begin=begin
        self.end=end
        self.strand=strand
        self.profile = Profile(tools)

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def getLocusInIGVJSFormat(self):
        return "'%s:%d-%d'"%(self.chromosome, self.begin, self.end)

    def getASArrayForHOT(self):
        return ["'%s'"%self.id, self.getLocusInIGVJSFormat(), "'%s'"%self.strand] + self.profile.getProfileAsArrayForHot()



class Gene(Feature):
    """
    Represents a gene - the same as a Feature, but it has several transcripts (features inside it)
    """
    def __init__(self, id, chromosome, begin, end, strand, tools):
        Feature.__init__(self, id, chromosome, begin, end, strand, tools)
        self.transcriptId2Transcript={}