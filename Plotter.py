import math
from decimal import Decimal
import traceback
import plotly

class Category:
    def __init__(self, start, end, step):
        if type(start) is float or type(end) is float or type(step) is float:
            raise Exception("Do not user float as type to class Category, it is error-prone. Use int for integer intervals or decimal.Decimal for real intervals.")

        self.intervals=[]
        while start<end:
            self.intervals.append({"min": start, "max": start+step, "count":0})
            start+=step

    def getCategoriesAsString(self, displayInterval=False, displayPlusOnFirstItem=False, displayPlusOnLastItem=False):
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
                prefix = "("
                if lowerBound > 0:
                    prefix += "+"

                # set the suffix
                suffix = ""
                if (displayPlusOnFirstItem and i == 0) or (displayPlusOnLastItem and i == len(self.intervals) - 1):
                    suffix += "+"
                suffix += ")"

                intervalsAsString.append(prefix + str(lowerBound) + suffix)
            else:
                if i < len(self.intervals) - 1 or not displayPlusOnLastItem:
                    intervalsAsString.append("[%s,%s)" % (interval["min"], interval["max"]))
                else:
                    intervalsAsString.append("%s+"%(interval["min"]))

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

    def __buildPlots(self, fig, name):
        """
        build the plots and return what needs to be returned
        """
        plotly.offline.plot(fig, image = 'png', image_filename="%s/%s.png"%(self.plotsOutput, name),
                            filename="%s/%s.html"%(self.plotsOutput, name), auto_open=False)

        return {
            #"imagePlot": "<img src=plots/%s.png />" % name,
            "imagePlot": "TODO",
            "jsPlot": plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
        }

    def __produceBarPlot(self, name, tool2Categories, xlabel, ylabel, displayInterval=False, displayPlusOnFirstItem=False, displayPlusOnLastItem=False):
        #produce the plot
        data = [plotly.graph_objs.Bar(
                x=tool2Categories[tool].getCategoriesAsString(displayInterval, displayPlusOnFirstItem, displayPlusOnLastItem),
                y=tool2Categories[tool].getIntervalCount(),
                name=tool)
                    for tool in self.toolsNoRaw]

        layout = plotly.graph_objs.Layout(
            xaxis={"title": xlabel},
            yaxis={"title": ylabel},
            barmode='group'
        )

        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        return self.__buildPlots(fig, name)

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
        return self.__produceBarPlot("DifferenceOnTheNumberOfIsoformsPlot", tool2DifferenceCategories, "Difference on the number of isoforms", "Number of genes", False, True, True)


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
        return self.__produceBarPlot("LostTranscriptInGenesWSP2Plot", tool2RelativeTranscriptOfLostTranscriptCategories, "Relative transcript coverage in relation to gene coverage", "Number of transcripts", True)

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
        name = "DifferencesInRelativeExpressionsBoxPlot"


        #put the labels
        data=[plotly.graph_objs.Box(y=tool2DifferencesInRelativeExpressions[tool], name=tool) for tool in self.toolsNoRaw]
        layout = plotly.graph_objs.Layout(
            xaxis={"title": "Relative expression"},
            yaxis={"title": "Tools"}
        )
        fig = plotly.graph_objs.Figure(data=data, layout=layout)

        return self.__buildPlots(fig, name)

    def makeScatterPlotSizeParalogFamilies(self, geneID2gene, paralogous, disregardUnchangedGeneFamilies=False, includeOnlyCommonGenes=False):
        try:
            name = "ScatterPlotSizeParalogFamilies"
            if disregardUnchangedGeneFamilies:
                name+="_NoDiagPoints"
            if includeOnlyCommonGenes:
                name += "_OnlyCommonGenes"

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

            #get all the data to plot it
            tool2PlotData={tool:{} for tool in self.toolsNoRaw}
            for tool in self.toolsNoRaw:
                paralogousGeneFamilySizeAfterCorrection = get_paralogousGenesFamilySizeInTool(paralogousGroups, tool)

                #get the data points:
                #x = family size before correction
                #y = family size after correction
                #1/ (0,0) are excluded
                #2/ if disregardUnchangedGeneFamilies is True, then gene families with the same size are disregarded
                #3/ if includeOnlyCommonGenes is True, then only gene families present before and after are considered
                dataPoints=[]
                for i in xrange(len(paralogousGeneFamilySizeBeforeCorrection)):
                    if (includeOnlyCommonGenes and paralogousGeneFamilySizeBeforeCorrection[i]>0 and paralogousGeneFamilySizeAfterCorrection[i]>0) or \
                       (not includeOnlyCommonGenes and (paralogousGeneFamilySizeBeforeCorrection[i]>0 or paralogousGeneFamilySizeAfterCorrection[i]>0)):
                        if not disregardUnchangedGeneFamilies or \
                           (disregardUnchangedGeneFamilies and paralogousGeneFamilySizeBeforeCorrection[i]!=paralogousGeneFamilySizeAfterCorrection[i]):
                            dataPoints.append((paralogousGeneFamilySizeBeforeCorrection[i], paralogousGeneFamilySizeAfterCorrection[i]))

                dataPoint2Count={}
                for dataPoint in dataPoints:
                    if dataPoint not in dataPoint2Count:
                        dataPoint2Count[dataPoint]=dataPoints.count(dataPoint)

                #get the new datapoints
                dataPoints=dataPoint2Count.keys()
                tool2PlotData[tool]["xDataPoints"] = [x for x,y in dataPoints]
                tool2PlotData[tool]["yDataPoints"] = [y for x, y in dataPoints]
                #the counts will be the colors of the scatterplot
                tool2PlotData[tool]["count"] = dataPoint2Count.values()
            largestFamilySize=max([max(max(plotData["xDataPoints"]), max(plotData["yDataPoints"])) for plotData in tool2PlotData.values()])
            largestDatapointCount=max([max(plotData["count"]) for plotData in tool2PlotData.values()])


            #plot the data
            nbOfColumnsInSubplot = 3
            nbRowsInSubplot = int(math.ceil(float(len(self.toolsNoRaw)) / nbOfColumnsInSubplot))
            fig = plotly.tools.make_subplots(rows=nbRowsInSubplot, cols=nbOfColumnsInSubplot)
            for toolIndex, tool in enumerate(self.toolsNoRaw):
                row, col = int(toolIndex / nbOfColumnsInSubplot) + 1, toolIndex % nbOfColumnsInSubplot + 1
                #plot it
                trace = plotly.graph_objs.Scatter(x=tool2PlotData[tool]["xDataPoints"], y=tool2PlotData[tool]["yDataPoints"],
                                                  mode='markers',
                                                  marker={'color': tool2PlotData[tool]["count"],
                                                          'cmax': largestDatapointCount,
                                                          'cmin': 0,
                                                          'colorbar': {'title': 'Count'},
                                                          'colorscale': 'RdBu'}
                                                  )
                fig.append_trace(trace, row, col)

            fig['layout'].update(height=nbRowsInSubplot*600, width=nbOfColumnsInSubplot*600, showlegend=False)

            for toolIndex, tool  in enumerate(self.toolsNoRaw):
                fig['layout']['xaxis%d'%(toolIndex+1)].update(range=[0, largestFamilySize + 1], title="Raw")
                fig['layout']['yaxis%d'%(toolIndex+1)].update(range=[0, largestFamilySize + 1], title=tool)
                fig['layout']['shapes'].append(dict({
                        'xref': "x%d"%(toolIndex+1),
                        'yref': "y%d"%(toolIndex+1),
                        'type': 'line',
                        'x0': 0,
                        'y0': 0,
                        'x1': largestFamilySize,
                        'y1': largestFamilySize,
                        'opacity': 0.5
                    }))


            return self.__buildPlots(fig, name)
        except:
            traceback.print_exc()
            #TODO: treat this better - we report error on computing all plots even if a single plot fails
            return {
                "imagePlot": "<p style='color: red; font-size: large;'>Error on computing this plot...</p>",
                "jsPlot": "<p style='color: red; font-size: large;'>Error on computing this plot...</p>"
            }


    def makeScatterPlotCoverageOfMainIsoform(self, geneID2gene):
        name = "ScatterPlotCoverageOfMainIsoform"

        # get all the data to plot it
        tool2PlotData = {tool: {} for tool in self.toolsNoRaw}
        for tool in self.toolsNoRaw:
            # get the data points:
            # x = coverage before correction
            # y = coverage after correction
            tool2PlotData[tool]["xDataPoints"] = []
            tool2PlotData[tool]["yDataPoints"] = []
            for gene in geneID2gene.values():
                mainIsoform = gene.getMainIsoform()
                if mainIsoform.profile.isExpressedInTool("raw.bam"): #filtering out the cases where the main isoform has 0 reads mapping to it
                    tool2PlotData[tool]["xDataPoints"].append(mainIsoform.profile.tool2NbOfMappedReads["raw.bam"])
                    tool2PlotData[tool]["yDataPoints"].append(mainIsoform.profile.tool2NbOfMappedReads[tool])

        highestExpression = max([max(max(plotData["xDataPoints"]), max(plotData["yDataPoints"])) for plotData in tool2PlotData.values()])

        # plot the data
        nbOfColumnsInSubplot = 3
        nbRowsInSubplot = int(math.ceil(float(len(self.toolsNoRaw)) / nbOfColumnsInSubplot))
        fig = plotly.tools.make_subplots(rows=nbRowsInSubplot, cols=nbOfColumnsInSubplot,
                                         subplot_titles=self.toolsNoRaw)
        for toolIndex, tool in enumerate(self.toolsNoRaw):
            row, col = int(toolIndex / nbOfColumnsInSubplot) + 1, toolIndex % nbOfColumnsInSubplot + 1
            # plot it
            trace = plotly.graph_objs.Scatter(x=tool2PlotData[tool]["xDataPoints"], y=tool2PlotData[tool]["yDataPoints"],
                                              mode='markers',
                                              marker={
                                                  'color': 'black',
                                                  'opacity': 0.2
                                              }
                                              )
            fig.append_trace(trace, row, col)

        fig['layout'].update(height=nbRowsInSubplot * 600, width=nbOfColumnsInSubplot * 600, showlegend=False)

        for toolIndex in xrange(len(self.toolsNoRaw)):
            fig['layout']['xaxis%d' % (toolIndex + 1)].update(range=[0, int(math.ceil(highestExpression*1.1))+1], title="Main coverage before")
            fig['layout']['yaxis%d' % (toolIndex + 1)].update(range=[0, int(math.ceil(highestExpression*1.1))+1], title="Main coverage after")
            fig['layout']['shapes'].append(dict({
                'xref': "x%d"%(toolIndex+1),
                'yref': "y%d"%(toolIndex+1),
                'type': 'line',
                'x0': 0,
                'y0': 0,
                'x1': highestExpression,
                'y1': highestExpression,
                'opacity': 0.5
            }))

        return self.__buildPlots(fig, name)

    def makeBarPlotFromStats(self, statProfiler, metric):
        name = metric

        #produce the plot
        data = [plotly.graph_objs.Bar(x=statProfiler.tools, y=[statProfiler.tool2Stats[tool][metric] for tool in statProfiler.tools])]
        layout = plotly.graph_objs.Layout(
            title=metric,
            xaxis={"title": "Tools"},
            yaxis={"title": metric}
        )

        return {
            "imagePlot": "<img src=plots/%s.png />" % name,
            "jsPlot": plotly.offline.plot({"data": data, "layout": layout}, include_plotlyjs=False, output_type='div')
        }

    def makeBarPlotFromStats(self, statProfiler, metric):
        name = metric

        #produce the plot
        data = [plotly.graph_objs.Bar(x=["Tools"], y=[statProfiler.tool2Stats[tool][metric]], name=tool) for tool in statProfiler.tools]
        layout = plotly.graph_objs.Layout(
            title=metric,
            yaxis={"title": metric},
            barmode="group"
        )
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        return self.__buildPlots(fig, name)

    def makeReadCountPlotDividedBySize(self, statProfiler, feature, title):
        #first we get the labels
        labels = statProfiler.tool2Stats["raw.bam"][feature].getCategoriesAsString(displayInterval=True, displayPlusOnLastItem=True)


        # produce the plot data
        data=[]
        for tool in self.tools:
            data.append(plotly.graph_objs.Scatter(
                x=range(len(labels)),
                y=statProfiler.tool2Stats[tool][feature].getIntervalCount(),
                mode='lines+markers',
                name="%s"%(tool)
                ))


        layout = plotly.graph_objs.Layout(
            title=title,
            xaxis=plotly.graph_objs.XAxis(
                   title="Read lengths",
                   showticklabels=True,
                   tickvals=range(len(labels)),
                   ticktext=labels
                ),
            yaxis={"title": "Read count"}
        )

        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        return self.__buildPlots(fig, "ReadCountPlotDividedBySize_%s"%feature)
