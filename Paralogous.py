#!/usr/bin/env python
# -*- coding: utf-8 -*-


class ParalogousGroup:
    """
    Represent a paralogous group or family
    """
    def __init__(self, id, group=[]):
        self.id = id #family ID, an int
        self.group=group #the genes in this family, list of Feature

    def addGenesFromOtherGroup(self, otherGroup):
        self.group+=otherGroup.group

    def clearGroup(self):
        del self.group[:]

    def getDescription(self):
        return "Family %d: %s" % (self.id, ", ".join([gene.id for gene in self.group]))

    def getGeneFamilySizeInTool(self, tool):
        paralogousGeneFamilySize = 0
        for gene in self.group:
            if gene.profile.isExpressedInTool(tool):
                paralogousGeneFamilySize += 1
        return paralogousGeneFamilySize

    def isExpressedInTool(self, tool):
        """
        Check if any of the genes is expressed in the given tool
        """
        return any([gene.profile.isExpressedInTool(tool) for gene in self.group])

    def isExpressedInAnyTool(self):
        """
        Check if any of the genes is expressed in any tool
        """
        return any([gene.profile.isExpressedInAnyTool() for gene in self.group])

    def getNbOfIsoformsExpressedInTool(self, tool):
        """
        Get the nb of isoforms expressed in the given tool for all genes of the family
        """
        return sum([gene.getNbOfIsoformsExpressedInTool(tool) for gene in self.group])

    def getDataToShowInPlot(self):
        return self.getDescription()

class Paralogous:
    """
    Represents the groups of paralogous genes
    Each group of paralogous genes has an ID and its genes
    """
    def __init__(self, geneID2gene):
        self.geneID2gene = geneID2gene
        self.paralogousGeneId2IdGroup={geneId:i for i, geneId in enumerate(geneID2gene.keys())} #create a union-find structure
        self.idGroup2paralogousGroup = {i:ParalogousGroup(i, [geneID2gene[geneId]]) for i, geneId in enumerate(geneID2gene.keys())}  # create a union-find structure


    @staticmethod
    def getErrorMessage():
        return "<p style='color: red; font-size: large;'>Paralogous file (--paralogous parameter) was not given, so we did not produce this plot. </p>"


    def __addParalogousRelation(self, geneId1, geneId2):
        """
        Add the information that geneId1 and geneId2 are paralogous
        """
        if geneId1 not in self.paralogousGeneId2IdGroup:
            print "WARNING: %s is in the paralogous file but not in the gtf file! Are you using the same references?" % geneId1
            return

        if geneId2 not in self.paralogousGeneId2IdGroup:
            print "WARNING: %s is in the paralogous file but not in the gtf file! Are you using the same references?" % geneId2
            return


        idGroup1 = self.paralogousGeneId2IdGroup[geneId1]
        idGroup2 = self.paralogousGeneId2IdGroup[geneId2]
        paralogousGroup1 = self.idGroup2paralogousGroup[idGroup1]
        paralogousGroup2 = self.idGroup2paralogousGroup[idGroup2]
        if idGroup1<idGroup2:
            self.paralogousGeneId2IdGroup[geneId2] = idGroup1
            paralogousGroup1.addGenesFromOtherGroup(paralogousGroup2)
            paralogousGroup2.clearGroup()
        elif idGroup1>idGroup2:
            self.paralogousGeneId2IdGroup[geneId1] = idGroup2
            paralogousGroup2.addGenesFromOtherGroup(paralogousGroup1)
            paralogousGroup1.clearGroup()

    def readParalogousFile(self, filename):
        """
        Read a paralogous file and populate this object
        """
        with open(filename) as file:
            header=True
            for line in file:
                if header: #skip the header
                    header=False
                    continue
                lineSplit = line.split()
                self.__addParalogousRelation(lineSplit[0], lineSplit[1])


    def getParalogousGroups(self):
        """
        :return: a list where each element is a ParalogousGroup
        """
        return self.idGroup2paralogousGroup.values()
