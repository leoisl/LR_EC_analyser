# make the plot nb_of_isoforms_lost_or_won_nb_of_genes.80QC_filter.pdf

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import mpld3
from decimal import Decimal

class Category:
    def __init__(self, start, end, step):
        if type(start) is float or type(end) is float or type(step) is float:
            raise Exception("Do not user float as type to class Category, it is error-prone. Use int for integer intervals or decimal.Decimal for real intervals.")

        self.intervals=[]
        while start<end:
            self.intervals.append({"min": start, "max": start+step, "count":0})
            start+=step

    def getCategoriesAsString(self, displayInterval=False):
        """
        transform the intervals list into a list of string that will be the xlabels of the plot
        Basically, transform a vector like [-2, -1, 0, 1, 2] into ["-2+", "-1", "0", "+1", "+2+"]
        :return: list of string
        """
        intervalsAsString = []
        for i, interval in enumerate(self.intervals):
            if not displayInterval:
                lowerBound = interval["min"]
                # set the prefix
                prefix = ""
                if lowerBound > 0:
                    prefix = "+"

                # set the suffix
                suffix = ""
                if i == 0 or i == len(self.intervals) - 1:
                    suffix = "+"

                intervalsAsString.append(prefix + str(lowerBound) + suffix)
            else:
                intervalsAsString.append("[%s,%s)"%(interval["min"], interval["max"]))

        return intervalsAsString

    def addDataPoint(self, point):
        """
        :param point: add a number to a category
        :return:
        """
        nbOfCategoriesItFit = 0

        if point < self.intervals[0]["min"]:
            self.intervals[0]["count"]+=1
            nbOfCategoriesItFit += 1
        if point >= self.intervals[-1]["max"]:
            self.intervals[-1]["count"]+=1
            nbOfCategoriesItFit += 1

        for i, interval in enumerate(self.intervals):
            if point >= interval['min'] and point < interval['max']:
                interval["count"] += 1
                nbOfCategoriesItFit += 1

        if nbOfCategoriesItFit != 1:
            raise Exception("ERROR: %d fit %d categories..." % (point, nbOfCategoriesItFit))

    def __len__(self):
        return len(self.intervals)

    def getIntervalCount(self):
        return [ interval["count"] for interval in self.intervals ]

class Plotter:
    """
    Makes several plots
    """
    def __init__(self, tools):
        self.tools=tools
        self.toolsNoRaw=list(tools)
        self.toolsNoRaw.remove("raw.bam")

    def __producePlotAsHTML(self, tool2DifferenceCategories, blankSpace, xlabel, ylabel, displayInterval=False):
        # produce the plot
        # fig = plt.figure(figsize=(10, 5))
        fig = plt.figure()

        #put the labels
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        #compute the indexes and bar widths
        xAxisIndexes = np.arange(len(tool2DifferenceCategories[self.toolsNoRaw[0]]))
        barWidth = 1.0 / len(self.toolsNoRaw) - blankSpace

        #add each bar
        for index, tool in enumerate(self.toolsNoRaw):
            plt.bar(xAxisIndexes + (index * barWidth), tool2DifferenceCategories[tool].getIntervalCount(),
                    width=barWidth, label=tool)

        #add the x labels
        plt.xticks(xAxisIndexes + (len(self.toolsNoRaw) / 2.0 * barWidth - barWidth / 2.0),
                   tool2DifferenceCategories[self.toolsNoRaw[0]].getCategoriesAsString(displayInterval))

        #add the tool labels
        plt.legend()

        # save plot to html
        return mpld3.fig_to_html(fig)

    def makeDifferenceOnTheNumberOfIsoformsPlot(self, geneID2gene, lowestCategory=-3, highestCategory=3, step=1, blankSpace=0.1):
        """

        :param geneID2gene:
        :param lowestCategory:
        :param highestCategory:
        :param blankSpace:
        :return: a string with html code to be put in the html report
        """


        '''
        Helper functions
        '''
        def addDifferenceOfIsoformNumber(gene, tool, tool2DifferencesCategories):
            nbOfIsoformsExpressedInRaw = gene.getNbOfIsoformsExpressedInTool("raw.bam")
            nbOfIsoformsExpressedInTool = gene.getNbOfIsoformsExpressedInTool(tool)
            if nbOfIsoformsExpressedInRaw > 0 or nbOfIsoformsExpressedInTool > 0:
                # if there is any expression in raw or tool, we add it!
                tool2DifferencesCategories[tool].addDataPoint(nbOfIsoformsExpressedInRaw - nbOfIsoformsExpressedInTool)

        #builds tool2DifferencesCategories - for each tool, we have an array with (raw.expression - tool.expression) for each gene and transcript
        def get_tool2DifferenceCategories():
            tool2DifferenceCategories={tool:Category(lowestCategory, highestCategory+1, step) for tool in self.toolsNoRaw}

            #populate tool2DifferenceCategories
            for gene in geneID2gene.values():
                if gene.profile.isExpressedInAnyTool():
                    for tool in self.toolsNoRaw:
                        addDifferenceOfIsoformNumber(gene, tool, tool2DifferenceCategories)

            return tool2DifferenceCategories
        '''
        Helper functions - END
        '''


        tool2DifferenceCategories = get_tool2DifferenceCategories()
        return self.__producePlotAsHTML(tool2DifferenceCategories, blankSpace, "Difference on the number of isoforms", "Number of genes")


    def makeLostTranscriptInGenesWSP2Plot(self, geneID2gene, blankSpace=0.1):
        '''
        Helper functions
        '''
        # builds tool2RelativeTranscriptOfLostTranscriptCategories
        def get_tool2RelativeTranscriptOfLostTranscriptCategories():
            tool2RelativeTranscriptOfLostTranscriptCategories = {tool: Category(Decimal("0.0"), Decimal("1.0"), Decimal("0.1")) for tool in self.toolsNoRaw}

            # populate tool2RelativeTranscriptOfLostTranscriptCategories
            for gene in geneID2gene.values():
                if gene.profile.isExpressedInAnyTool():
                    for transcript in gene.transcriptId2Transcript.values():
                        if transcript.profile.isExpressedInAnyTool():
                            for tool in self.toolsNoRaw:
                                if transcript.profile.isExpressedInTool("raw.bam") and not transcript.profile.isExpressedInTool(tool):
                                    #if the transcript is in raw, but it is not in the tool, then this transcript "disappeared"
                                    #add the relative transcript coverage in raw dataset
                                    tool2RelativeTranscriptOfLostTranscriptCategories[tool].addDataPoint( float(transcript.profile.tool2NbOfMappedReads["raw.bam"]) / float(gene.profile.tool2NbOfMappedReads["raw.bam"]))


            return tool2RelativeTranscriptOfLostTranscriptCategories

        '''
        Helper functions
        '''
        tool2RelativeTranscriptOfLostTranscriptCategories =get_tool2RelativeTranscriptOfLostTranscriptCategories()
        return self.__producePlotAsHTML(tool2RelativeTranscriptOfLostTranscriptCategories, blankSpace, "Relative transcript coverage in relation to gene coverage", "Number of transcripts", True)