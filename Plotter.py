# make the plot nb_of_isoforms_lost_or_won_nb_of_genes.80QC_filter.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import mpld3
import math
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
    def __init__(self, tools, plotsOutput):
        self.tools=tools
        self.toolsNoRaw=list(tools)
        self.toolsNoRaw.remove("raw.bam")
        self.plotsOutput = plotsOutput

    def __produceBarPlotAsHTML(self, tool2Categories, blankSpace, xlabel, ylabel, displayInterval=False):
        # produce the plot
        fig = plt.figure(figsize=(10, 5))

        #put the labels
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        #compute the indexes and bar widths
        xAxisIndexes = np.arange(len(tool2Categories[self.toolsNoRaw[0]]))
        barWidth = 1.0 / len(self.toolsNoRaw) - blankSpace/len(self.toolsNoRaw)

        #add each bar
        for index, tool in enumerate(self.toolsNoRaw):
            plt.bar(xAxisIndexes + (index * barWidth), tool2Categories[tool].getIntervalCount(),
                    width=barWidth, label=tool)

        #add the x labels
        plt.xticks(xAxisIndexes + (len(self.toolsNoRaw) / 2.0 * barWidth - barWidth / 2.0),
                   tool2Categories[self.toolsNoRaw[0]].getCategoriesAsString(displayInterval))

        #add the tool labels
        plt.legend()

        # save plot to html
        return mpld3.fig_to_html(fig, d3_url="lib/js/d3.v3.min.js", mpld3_url="lib/js/mpld3.v0.3.min.js")

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
                tool2DifferencesCategories[tool].addDataPoint(nbOfIsoformsExpressedInTool - nbOfIsoformsExpressedInRaw)

        #builds tool2DifferencesCategories - for each tool, we have an array with (tool.expression - raw.expression) for each gene and transcript
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
        return self.__produceBarPlotAsHTML(tool2DifferenceCategories, blankSpace, "Difference on the number of isoforms", "Number of genes")


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
                    for tool in self.toolsNoRaw:
                        if gene.profile.isExpressedInTool("raw.bam") and gene.profile.isExpressedInTool(tool): #the gene must be expressed in the raw and in the tool - otherwise we are not looking at the lost transcript, but rather at lost genes
                            for transcript in gene.transcriptId2Transcript.values():
                                if transcript.profile.isExpressedInTool("raw.bam") and not transcript.profile.isExpressedInTool(tool):
                                    #if the transcript is in raw, but it is not in the tool, then this transcript "disappeared"
                                    #add the relative transcript coverage in raw dataset
                                    tool2RelativeTranscriptOfLostTranscriptCategories[tool].addDataPoint(transcript.computeRelativeExpression("raw.bam"))


            return tool2RelativeTranscriptOfLostTranscriptCategories

        '''
        Helper functions
        '''
        tool2RelativeTranscriptOfLostTranscriptCategories =get_tool2RelativeTranscriptOfLostTranscriptCategories()
        return self.__produceBarPlotAsHTML(tool2RelativeTranscriptOfLostTranscriptCategories, blankSpace, "Relative transcript coverage in relation to gene coverage", "Number of transcripts", True)

    def makeDifferencesInRelativeExpressionsBoxPlot(self, geneID2gene, blankSpace=0.1):
        def get_tool2DifferencesInRelativeExpressions():
            tool2DifferencesInRelativeExpressions={tool:[] for tool in self.toolsNoRaw}
            for gene in geneID2gene.values():
                if gene.profile.isExpressedInAnyTool():
                    for tool in self.toolsNoRaw:
                        if gene.profile.isExpressedInTool("raw.bam") or gene.profile.isExpressedInTool(tool):  # the gene must be expressed in the raw or in the tool
                            for transcript in gene.transcriptId2Transcript.values():
                                relativeExpressionBefore = transcript.computeRelativeExpression("raw.bam")
                                relativeExpressionAfter = transcript.computeRelativeExpression(tool)
                                differenceInRelativeExpressions = abs(relativeExpressionBefore-relativeExpressionAfter)
                                tool2DifferencesInRelativeExpressions[tool].append(differenceInRelativeExpressions)
            return tool2DifferencesInRelativeExpressions

        tool2DifferencesInRelativeExpressions = get_tool2DifferencesInRelativeExpressions()

        #make the boxplot
        fig = plt.figure(figsize=(10, 5))

        #put the labels
        plt.ylabel("Tools")
        plt.xlabel("Relative expression")
        data=[tool2DifferencesInRelativeExpressions[tool] for tool in self.toolsNoRaw]
        plt.boxplot(data, labels=self.toolsNoRaw, sym='', vert=False)

        # save plot to html
        return mpld3.fig_to_html(fig, d3_url="lib/js/d3.v3.min.js", mpld3_url="lib/js/mpld3.v0.3.min.js")

    def makeScatterPlotSizeParalogFamilies(self, geneID2gene, paralogous):
        def get_paralogousGenesFamilySizeInTool(paralogousGroups, tool):
            paralogousGenesFamilySize=[]
            for paralogousGroup in paralogousGroups:
                paralogousGeneFamilySize = 0
                for geneId in paralogousGroup:
                    if geneID2gene[geneId].profile.isExpressedInTool(tool):
                        paralogousGeneFamilySize += 1
                paralogousGenesFamilySize.append(paralogousGeneFamilySize)

            return paralogousGenesFamilySize

        paralogousGroups = paralogous.getParalogousGroups()

        paralogousGeneFamilySizeBeforeCorrection = get_paralogousGenesFamilySizeInTool(paralogousGroups, "raw.bam")

        nbOfColumnsInSubplot = 3
        nbRowsInSubplot = int(math.ceil(float(len(self.toolsNoRaw))/nbOfColumnsInSubplot))
        fig = plt.figure(figsize=(5*nbOfColumnsInSubplot, 5*nbRowsInSubplot))
        for toolIndex, tool in enumerate(self.toolsNoRaw):
            paralogousGeneFamilySizeAfterCorrection = get_paralogousGenesFamilySizeInTool(paralogousGroups, tool)

            #removing the gene families where we have 0s in both tools
            paralogousGeneFamilySizeBeforeCorrectionNoZeros=[]
            paralogousGeneFamilySizeAfterCorrectionNoZeros=[]
            for i in xrange(len(paralogousGeneFamilySizeBeforeCorrection)):
                if paralogousGeneFamilySizeBeforeCorrection[i]>0 or paralogousGeneFamilySizeAfterCorrection[i]>0:
                    paralogousGeneFamilySizeBeforeCorrectionNoZeros.append(paralogousGeneFamilySizeBeforeCorrection[i])
                    paralogousGeneFamilySizeAfterCorrectionNoZeros.append(paralogousGeneFamilySizeAfterCorrection[i])

            plt.subplot(nbRowsInSubplot, nbOfColumnsInSubplot, toolIndex+1)
            plt.scatter(paralogousGeneFamilySizeAfterCorrectionNoZeros, paralogousGeneFamilySizeBeforeCorrectionNoZeros, alpha=0.5, label="Gene family")

            #drawing the diagonal line
            maxValue = max(max(paralogousGeneFamilySizeBeforeCorrectionNoZeros), max(paralogousGeneFamilySizeAfterCorrectionNoZeros))
            plt.plot(range(maxValue+1), range(maxValue+1), alpha=0.5, color="black")

            plt.xlabel("Raw")
            plt.ylabel(tool)
            plt.legend()

        plt.savefig(self.plotsOutput+"/scatterPlotSizeParalogFamilies.png")
        return "<img src=plots/scatterPlotSizeParalogFamilies.png />"
        #return mpld3.fig_to_html(fig, d3_url="lib/js/d3.v3.min.js", mpld3_url="lib/js/mpld3.v0.3.min.js")
