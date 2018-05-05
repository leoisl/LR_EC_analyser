#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Paralogous:
    """
    Represents the groups of paralogous genes
    Each group of paralogous genes has an ID and its genes
    """
    def __init__(self, geneID2gene):
        self.paralogousGeneId2IdGroup={geneId:i for i, geneId in enumerate(geneID2gene.keys())} #create a union-find structure
        self.groupSize = len(geneID2gene.keys())


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
        minId = min(idGroup1, idGroup2)
        self.paralogousGeneId2IdGroup[geneId1] = minId
        self.paralogousGeneId2IdGroup[geneId2] = minId


    def __fixGroupsId(self):
        """
        Fix the groups id, i.e., if the groups id are 5, 8 and 9, after this function call, they will be 0, 1, 2
        """
        # TODO: recode this, it is too bad

        #get the sorted unique groups id
        sortedGroupsId = sorted(list(self.paralogousGeneId2IdGroup.values()))
        uniqueGroupsId = []
        for id in sortedGroupsId:
            if id not in uniqueGroupsId:
                uniqueGroupsId.append(id)

        #update the groups id
        oldId2NewId = {id : i for i, id in enumerate(uniqueGroupsId)}
        for geneId in self.paralogousGeneId2IdGroup:
            self.paralogousGeneId2IdGroup[geneId]=oldId2NewId[self.paralogousGeneId2IdGroup[geneId]]

        #update the size
        self.groupSize = max(oldId2NewId.values())+1


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
        self.__fixGroupsId()  # fix the groups' ids


    def getParalogousGroups(self):
        """
        :return: a list where each element is a paralogous group (a list of paralogous genes)
        """
        return [ [geneId for index, geneId in enumerate(self.paralogousGeneId2IdGroup.keys()) if index == groupId] for groupId in xrange(self.groupSize) ]
