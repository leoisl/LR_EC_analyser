import math
from decimal import Decimal
import plotly
from scipy import stats
import numpy
import re
import json
divIdCapturePattern = re.compile("id=\"(.*?)\"")

class Category:
    def __init__(self, start, end, step):
        if type(start) is float or type(end) is float or type(step) is float:
            raise Exception("Do not user float as type to class Category, it is error-prone. Use int for integer intervals or decimal.Decimal for real intervals.")

        self.intervals=[]
        while start<end:
            self.intervals.append({"min": start, "max": start+step, "data": []})
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

    def addDataPoint(self, point, dataToSave):
        """
        :param point: it is the value of the dataToSave - we use this to infer the category
        :return:
        """
        nbOfCategoriesItFit = 0

        if point < self.intervals[0]["min"]:
            self.intervals[0]["data"].append(dataToSave)
            nbOfCategoriesItFit += 1
        if point >= self.intervals[-1]["max"]:
            self.intervals[-1]["data"].append(dataToSave)
            nbOfCategoriesItFit += 1

        for i, interval in enumerate(self.intervals):
            if point >= interval['min'] and point < interval['max']:
                interval["data"].append(dataToSave)
                nbOfCategoriesItFit += 1

        if nbOfCategoriesItFit != 1:
            raise Exception("ERROR: %d fit %d categories..." % (point, nbOfCategoriesItFit))

    def __len__(self):
        return len(self.intervals)

    def getIntervalCount(self, inPercentage=False):
        if not inPercentage:
            return [ len(interval["data"]) for interval in self.intervals ]
        else:
            total = sum( [ len(interval["data"]) for interval in self.intervals ] )
            return [float(len(interval["data"])) / float(total) * 100 for interval in self.intervals]

    def __repr__(self):
        return self.intervals.__repr__()

    def __str__(self):
        return str(self.intervals)

#TODO: refactor this by making a base class for Category and TextCategory
class TextCategory:
    """
    Class that represent categories, but as text
    """
    def __init__(self, categories):
        self.intervals=[{"category": category, "data": []} for category in categories]

    def getCategoriesAsString(self, displayInterval=False, displayPlusOnFirstItem=False, displayPlusOnLastItem=False):
        return [interval["category"] for interval in self.intervals]

    def addDataPointAndIAlreadyKnowTheCategory(self, category, dataToSave):
        categoriesAsString = self.getCategoriesAsString()
        if category not in categoriesAsString:
            raise Exception("ERROR: non-existing category: %s..." % category)
        else:
            self.intervals[categoriesAsString.index(category)]["data"].append(dataToSave)

    def __len__(self):
        return len(self.intervals)

    def getIntervalCount(self, inPercentage=False):
        if not inPercentage:
            return [ len(interval["data"]) for interval in self.intervals ]
        else:
            total = sum( [ len(interval["data"]) for interval in self.intervals ] )
            return [float(len(interval["data"])) / float(total) * 100 for interval in self.intervals]

    def __repr__(self):
        return self.intervals.__repr__()

    def __str__(self):
        return str(self.intervals)

class Plotter:
    """
    Makes several plots
    """
    def __init__(self, tools, hybridTools, selfTools):
        self.tools=tools
        self.toolsNoRaw=list(tools)
        self.toolsNoRaw.remove("raw.bam")
        self.hybridTools = hybridTools
        self.selfTools = selfTools

    def __getToolCategory(self, tool):
        if tool=="raw.bam":
            return "raw"
        elif tool in self.hybridTools:
            return "hybrid"
        elif tool in self.selfTools:
            return "self"
        else:
            raise Exception("Unknown category for tool %s"%tool)

    def __buildPlots(self, fig, name, label2ToolIndex2Data=None):
        """
        build the plots and return what needs to be returned
        """
        """
        TODO: removed png plots
        plotly.offline.plot(fig, image = 'png', image_filename="%s/%s.png"%(self.plotsOutput, name),
                            filename="%s/%s.html"%(self.plotsOutput, name), auto_open=False)
        """
        fig["layout"]["hovermode"] = "closest"


        #change the html plot code so that we can register mouse events
        #TODO: this is SO BAD, but the only way to make it work...
        #TODO: if it is bad, but works, is it really bad?
        htmlPlotCode = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
        match = re.search(divIdCapturePattern, htmlPlotCode)
        divId = match.group(1)
        htmlPlotCode = htmlPlotCode.replace("</script>", "; myPlot = document.getElementById('%s'); myPlot.on('plotly_click', function(data){ \
            showDataOnModal('%s', data);}); plotInfo['%s']=jQuery.parseJSON(%s);</script>" % (divId, name, name, json.dumps(json.dumps(label2ToolIndex2Data)))) #json.dumps 2 times to encode


        return {
            #"imagePlot": "<img src=plots/%s.png />" % name,
            "imagePlot": "TODO",
            "jsPlot": htmlPlotCode
        }

    def __produceBarPlot(self, name, tool2Categories, xlabel, ylabel, displayInterval=False, displayPlusOnFirstItem=False, displayPlusOnLastItem=False, generateDataToBeShown=False, inPercentage=False):
        #produce the plot
        xLabels = tool2Categories.values()[0].getCategoriesAsString(displayInterval, displayPlusOnFirstItem, displayPlusOnLastItem)

        data = [plotly.graph_objs.Bar(
                x=xLabels,
                y=tool2Categories[tool].getIntervalCount(inPercentage),
                name=tool)
                for tool in self.toolsNoRaw]

        layout = plotly.graph_objs.Layout(
            xaxis={"title": xlabel},
            yaxis={"title": ylabel},
            barmode='group'
        )
        fig = plotly.graph_objs.Figure(data=data, layout=layout)

        # generate label2ToolIndex2Data to put in the plot, if generateDataToBeShown is True
        label2ToolIndex2Data = None
        if generateDataToBeShown:
            label2ToolIndex2Data = {}
            for intervalIndex, xLabel in enumerate(xLabels):
                label2ToolIndex2Data[xLabel]={}
                for toolIndex, tool in enumerate(self.toolsNoRaw):
                    label2ToolIndex2Data[xLabel][toolIndex]="\n".join(tool2Categories[tool].intervals[intervalIndex]["data"])

        return self.__buildPlots(fig, name, label2ToolIndex2Data)

    def __makeDifferenceOnTheNumberOfIsoformsPlotCore(self, geneID2gene, lowestCategory, highestCategory, step, unionOrIntersection):
        """
        :return: a string with html code to be put in the html report
        """
        '''
        Helper functions
        '''
        def addDifferenceOfIsoformNumber(gene, tool, tool2DifferencesCategories, unionOrIntersection):
            nbOfIsoformsExpressedInRaw = gene.getNbOfIsoformsExpressedInTool("raw.bam")
            nbOfIsoformsExpressedInTool = gene.getNbOfIsoformsExpressedInTool(tool)

            if unionOrIntersection == "Union":
                if nbOfIsoformsExpressedInRaw > 0 or nbOfIsoformsExpressedInTool > 0:
                    # if there is any expression in raw or tool, we add it!
                    tool2DifferencesCategories[tool].addDataPoint(nbOfIsoformsExpressedInTool - nbOfIsoformsExpressedInRaw, gene.id)
            elif unionOrIntersection == "Intersection":
                if nbOfIsoformsExpressedInRaw > 0 and nbOfIsoformsExpressedInTool > 0:
                    # if there is expression in both raw and tool, we add it!
                    tool2DifferencesCategories[tool].addDataPoint(nbOfIsoformsExpressedInTool - nbOfIsoformsExpressedInRaw, gene.id)
            else:
                raise Exception("unionOrIntersection in makeDifferenceOnTheNumberOfIsoformsPlot() should be union or intersection")

        #builds tool2DifferencesCategories - for each tool, we have an array with (tool.expression - raw.expression) for each gene and transcript
        def get_tool2DifferenceCategories():
            tool2DifferenceCategories={tool:Category(lowestCategory, highestCategory+1, step) for tool in self.toolsNoRaw}

            #populate tool2DifferenceCategories
            for gene in geneID2gene.values():
                if gene.profile.isExpressedInAnyTool():
                    for tool in self.toolsNoRaw:
                        addDifferenceOfIsoformNumber(gene, tool, tool2DifferenceCategories, unionOrIntersection)

            return tool2DifferenceCategories
        '''
        Helper functions - END
        '''


        tool2DifferenceCategories = get_tool2DifferenceCategories()
        return self.__produceBarPlot("DifferenceOnTheNumberOfIsoforms%sPlot"%unionOrIntersection, tool2DifferenceCategories, "Difference on the number of isoforms", "Number of genes", \
                                     displayPlusOnFirstItem=True, displayPlusOnLastItem=True, generateDataToBeShown=True)


    def makeDifferenceOnTheNumberOfIsoformsPlot(self, geneID2gene, lowestCategory=-3, highestCategory=3, step=1):
        """
        :return: Two plots, one with the union (genes in the raw or in the tool) and the other with the intersection
        """
        return self.__makeDifferenceOnTheNumberOfIsoformsPlotCore(geneID2gene, lowestCategory, highestCategory, step, "Union"), \
               self.__makeDifferenceOnTheNumberOfIsoformsPlotCore(geneID2gene, lowestCategory, highestCategory, step, "Intersection")


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
                                    tool2RelativeTranscriptOfLostTranscriptCategories[tool].addDataPoint(transcript.computeRelativeExpression("raw.bam"), transcript.id)


            return tool2RelativeTranscriptOfLostTranscriptCategories

        '''
        Helper functions
        '''
        tool2RelativeTranscriptOfLostTranscriptCategories = get_tool2RelativeTranscriptOfLostTranscriptCategories()
        return self.__produceBarPlot("LostTranscriptInGenesWSP2Plot", tool2RelativeTranscriptOfLostTranscriptCategories, "Relative transcript coverage in relation to gene coverage", "Number of transcripts", displayInterval=True, generateDataToBeShown=True)

    def makeDifferencesInRelativeExpressionsBoxPlot(self, geneID2gene):
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
        data=[plotly.graph_objs.Box(y=tool2DifferencesInRelativeExpressions[tool], name=tool, boxpoints=False) for tool in self.toolsNoRaw]
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
            labels = ["Shrunk", "Unchanged", "Expanded"]
            tool2GeneralStatsCategories = {tool : TextCategory(labels) for tool in self.toolsNoRaw}
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
                            familyDescription = "Family %d: %s"%(i, ", ".join(paralogousGroups[i]))
                            if paralogousGeneFamilySizeBeforeCorrection[i] < paralogousGeneFamilySizeAfterCorrection[i]:
                                tool2GeneralStatsCategories[tool].addDataPointAndIAlreadyKnowTheCategory("Expanded", familyDescription)
                            elif paralogousGeneFamilySizeBeforeCorrection[i] > paralogousGeneFamilySizeAfterCorrection[i]:
                                tool2GeneralStatsCategories[tool].addDataPointAndIAlreadyKnowTheCategory("Shrunk", familyDescription)
                            else:
                                tool2GeneralStatsCategories[tool].addDataPointAndIAlreadyKnowTheCategory("Unchanged", familyDescription)


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
            nbOfColumnsInSubplot = 4
            nbRowsInSubplot = int(math.ceil(float(len(self.toolsNoRaw)) / nbOfColumnsInSubplot))
            specificPlotFig = plotly.tools.make_subplots(rows=nbRowsInSubplot, cols=nbOfColumnsInSubplot)
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
                specificPlotFig.append_trace(trace, row, col)

                specificPlotFig['layout'].update(height=nbRowsInSubplot*400, width=nbOfColumnsInSubplot*400, showlegend=False)

            for toolIndex, tool  in enumerate(self.toolsNoRaw):
                specificPlotFig['layout']['xaxis%d'%(toolIndex+1)].update(range=[0, largestFamilySize + 1], title="Raw")
                specificPlotFig['layout']['yaxis%d'%(toolIndex+1)].update(range=[0, largestFamilySize + 1], title=tool)
                specificPlotFig['layout']['shapes'].append(dict({
                        'xref': "x%d"%(toolIndex+1),
                        'yref': "y%d"%(toolIndex+1),
                        'type': 'line',
                        'x0': 0,
                        'y0': 0,
                        'x1': largestFamilySize,
                        'y1': largestFamilySize,
                        'opacity': 0.5
                    }))

            #build the plots for the general stats
            return self.__produceBarPlot(name + "General", tool2GeneralStatsCategories, "Tool's behaviour towards the gene family", "Gene family count in %", generateDataToBeShown=True, inPercentage=True), \
                   self.__buildPlots(specificPlotFig, name+"Specific")

        except ValueError as exception:
            if exception.message == "max() arg is an empty sequence":
                #expected error in some situations
                return {
                    "imagePlot": "<p style='color: red; font-size: large;'>Error on computing this plot...</p>",
                    "jsPlot": "<p style='color: red; font-size: large;'>Error on computing this plot...</p>"
                }, {
                    "imagePlot": "<p style='color: red; font-size: large;'>Error on computing this plot...</p>",
                    "jsPlot": "<p style='color: red; font-size: large;'>Error on computing this plot...</p>"
                }
            else:
                raise exception #unexpected error


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
        rSquaredAnnotations=[]
        nbOfColumnsInSubplot = 4
        nbRowsInSubplot = int(math.ceil(float(len(self.toolsNoRaw)) / nbOfColumnsInSubplot))
        fig = plotly.tools.make_subplots(rows=nbRowsInSubplot, cols=nbOfColumnsInSubplot,
                                         subplot_titles=self.toolsNoRaw)
        for toolIndex, tool in enumerate(self.toolsNoRaw):
            row, col = int(toolIndex / nbOfColumnsInSubplot) + 1, toolIndex % nbOfColumnsInSubplot + 1
            # plot it
            trace = plotly.graph_objs.Scatter(x=tool2PlotData[tool]["xDataPoints"], y=tool2PlotData[tool]["yDataPoints"],
											  type='scattergl',
                                              mode='markers',
                                              marker={
                                                  'color': 'black',
                                                  'opacity': 0.2
                                              }
                                              )
            fig.append_trace(trace, row, col)

            # get R squared, from https://plot.ly/python/linear-fits/
            slope, intercept, r_value, p_value, std_err = stats.linregress(tool2PlotData[tool]["xDataPoints"], tool2PlotData[tool]["yDataPoints"])
            rSquared = r_value ** 2
            linearRegressionLine = slope * numpy.array(tool2PlotData[tool]["xDataPoints"]) + intercept
            linearRegressionTrace = plotly.graph_objs.Scatter(
                x=tool2PlotData[tool]["xDataPoints"],
                y=linearRegressionLine,
                mode='lines',
                marker=plotly.graph_objs.Marker(color='rgb(31, 119, 180)'),
                name='Fit'
            )
            fig.append_trace(linearRegressionTrace, row, col)

            rSquaredAnnotations.append(plotly.graph_objs.Annotation(
                xref = "x%d"%(toolIndex + 1),
                yref = "y%d"%(toolIndex + 1),
                x = highestExpression/2,
                y = highestExpression,
                text="R^2 = %f"%rSquared,
                showarrow=False,
                font=plotly.graph_objs.Font(size=14)
            ))



        fig['layout'].update(height=nbRowsInSubplot * 400, width=nbOfColumnsInSubplot * 400, showlegend=False)

        for toolIndex, tool in enumerate(self.toolsNoRaw):
            fig['layout']['xaxis%d' % (toolIndex + 1)].update(range=[0, int(math.ceil(highestExpression*1.1))+1], title="Main isoform coverage before" if toolIndex==0 else "")
            fig['layout']['yaxis%d' % (toolIndex + 1)].update(range=[0, int(math.ceil(highestExpression*1.1))+1], title="Main isoform coverage after" if toolIndex==0 else "")
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

        #
        fig['layout']['annotations'].extend(rSquaredAnnotations)

        return self.__buildPlots(fig, name)

    def makeBarPlotFromStats(self, statProfiler, metric):
        name = statProfiler.getFeatureName(metric)

        #produce the plot
        data = [plotly.graph_objs.Bar(x=["Tools"], y=[statProfiler.getStatsForToolAndMetric(tool, metric)], name=tool) for tool in statProfiler.tools]
        layout = plotly.graph_objs.Layout(
            title=name,
            yaxis={"title": statProfiler.getNiceDescriptionForFeature(metric)},
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
