#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import shutil
from Profiler import *
from ExternalTools import *
from Parsers import *
from Plotter import *
from Paralogous import *
import os
import plotly
import sys

def testOrca():
    try:
        plotly.io.write_image(plotly.graph_objs.Figure(), "/dev/null", format="png")
    except ValueError:
        print "[FATAL ERROR] It seems orca is not working. If you are in a headless server configuration, try:"
        print "1. sudo apt-get install xvfb"
        print "2. cd <LR_EC_analyser_home>"
        print "3. bash fix_orca_headless.sh"
        print "And rerunning LR_EC_analyser. If this does not work, please fill an issue in the gitlab page"
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Long reads error corrector analyser.')
    parser.add_argument("--raw", dest="rawBam", help="The BAM file of the raw reads (i.e. the uncorrected long reads) mapped to the genome (preferably using gmap -n 10 -f samse).", required=True)
    parser.add_argument('--self', metavar='<self.bam>', type=str, nargs='+',
                        help='BAM files of the reads output by the SELF correctors mapped to the genome (preferably using gmap -n 10 -f samse).')
    parser.add_argument('--hybrid', metavar='<hybrid.bam>', type=str, nargs='+',
                        help='BAM files of the reads output by the HYBRID correctors mapped to the genome (preferably using gmap -n 10 -f samse).')

    parser.add_argument("--genome", dest="genome", help="The genome in Fasta file format.", required=True)
    parser.add_argument("--gtf", dest="gtf", help="The transcriptome in GTF file format.", required=True)
    parser.add_argument("--paralogous", help="A file where the first two columns denote paralogous genes (see file GettingParalogs.txt to know how you can get this file).")

    parser.add_argument("-o", dest="output", help="Output folder", default="output/")
    parser.add_argument("-t", dest="threads", type=int, help="Number of threads to use", default=1)
    parser.add_argument('--colours', metavar='<self.colours>', type=str, nargs='+',
                        help='A list of colours in hex encoding to use in the plots. Colour shading is nice to show related corrections (e.g. full-length, trimmed and split outputs from a same tool),'
                             'but unless the analysis is on few tools, it is hard to find a nice automated choice of colour shading. Hand-picking is more laborious but produces better results.'
                             'This parameter allows you to control the colors of each tool. The order of the tools are: raw -> hybrid -> self.'
                             'The hybrid and self ordering are given by parameter --hybrid and --self.'
                             'See an example of this parameter in https://gitlab.com/leoisl/LR_EC_analyser/blob/master/scripts/command_line_paper.sh .')
    parser.add_argument("--skip_bam_process", dest="skip_bam", action="store_true", help="Skips BAM processing (i.e. sorting and indexing BAM files) - assume we had already done this.")
    parser.add_argument("--skip_alignqc", dest="skip_alignqc", action="store_true",
                        help="Skips AlignQC calls - assume we had already done this.")
    parser.add_argument("--skip_copying", dest="skip_copying", action="store_true",
                        help="Skips copying genome and transcriptome to the output folder - assume we had already done this.")
    args=parser.parse_args()

    #test if Orca works
    testOrca()

    #create output dir
    if not os.path.exists(args.output):
        os.makedirs(args.output)


    #get some useful vars
    hybridBams = args.hybrid if args.hybrid != None else []
    selfBams = args.self if args.self != None else []

    if len(hybridBams) + len(selfBams) == 0:
        parser.print_help()
        print "\n[PARAMETER ERROR] - Please specify at least one BAM file for --self or --hybrid."
        sys.exit(1)

    bams = [args.rawBam] + hybridBams + selfBams
    hybridTools = [os.path.basename(bam) for bam in hybridBams]
    selfTools = [os.path.basename(bam) for bam in selfBams]
    tools = ["raw.bam"] + hybridTools + selfTools

    #fix args.colours, if it is not yet fixed
    if args.colours != None:
        for i, colour in enumerate(args.colours):
            if colour[0]!="#":
                args.colours[i] = "#"+colour


    #copy genome to output and index it
    genome = args.output + "/" + os.path.basename(args.genome)
    if not args.skip_copying:
        shutil.copy(args.genome, args.output)
        indexGenome(genome)
    else:
        print "Skipping genome copying and indexing..."

    #get genes from gtf
    geneID2gene = parseGTFToGetGenes(args.gtf, tools)

    # index and copy gtf to output
    gtf = args.output + "/" + os.path.basename(args.gtf)
    if not args.skip_copying:
        shutil.copy(args.gtf, args.output)
        sortAndIndexGTF(gtf)
    else:
        print "Skipping transcriptome copying and indexing..."

    #get the paralogous info, if given
    paralogous = None
    if args.paralogous != None:
        paralogous = Paralogous(geneID2gene)
        paralogous.readParalogousFile(args.paralogous)


    #sort and index bam
    if args.skip_bam:
        print "Skipping bam processing..."
    sortedBams = []
    for bam in bams:
        sortedBam = args.output+"/"+os.path.basename((bam)+".sorted.bam")
        sortedBams.append(sortedBam)
        if not args.skip_bam:
            processBam(bam, sortedBam, args.threads)


    #run AlignQC on the bams
    if not args.skip_alignqc:
        # AlignQC require the genome and gtf gzipped
        gzipFile(genome)
        gzipFile(gtf)

        #run alignQC in each tool
        for tool, bam in zip(tools, sortedBams):
            runAlignQC(tool, bam, genome+".gz", gtf+".gz", args.output, args.threads)
    else:
        print "Skipping AlignQC runs..."



    #create the output by parsing AlignQC results
    print "Running profilers - gathering and processing the data to produce the plots..."

    #stats profiler
    statProfiler = StatProfiler(tools, args.output)

    #gene profiler
    tool2Bam = {tool: os.path.basename(bam) for tool, bam in zip(tools, sortedBams)}
    geneProfiler = FeatureProfiler(geneID2gene, tools, tool2Bam, os.path.basename(genome), os.path.basename(gtf), args.output)
    #populate the gene profiler
    for tool in tools:
        geneProfiler.populateFromAnnotbest(tool, args.output)

    #SS profiler
    splicingSitesProfiler = SplicingSitesProfiler(tools, args.output)

    #ReadSetProfiler
    readSetProfiler = ReadSetProfiler(tools, args.output)

    print "Running profilers - gathering and processing the data to produce the plots - Done!"


    #create the Plotter and the plots
    print "Computing the plots..."
    plotsOutput = args.output+"/plots"
    if not os.path.exists(plotsOutput):
        os.makedirs(plotsOutput)
    plotter = Plotter(tools, hybridTools, selfTools, plotsOutput, args.colours)

    #make all stats plots
    allStatsPlots = {feature: plotter.makeBarPlotFromStats(statProfiler, feature) for feature in statProfiler.allFeatures}

    #make the ReadCuttingPlot
    totalReadsCuttingPlot = plotter.makeReadCountPlotDividedBySize(statProfiler, "TOTAL_SIZE_BINS", "All reads length line plot")
    alignedReadsCuttingPlot = plotter.makeReadCountPlotDividedBySize(statProfiler, "ALIGNED_SIZE_BINS", "Aligned reads length line plot")
    unalignedReadsCuttingPlot = plotter.makeReadCountPlotDividedBySize(statProfiler, "UNALIGNED_SIZE_BINS", "Unaligned reads length line plot")

    htmlDifferenceOnTheNumberOfIsoformsPlotIntersection, htmlDifferenceOnTheNumberOfIsoformsGeneFamilyPlotIntersection \
        = plotter.makeDifferenceOnTheNumberOfIsoformsPlot(geneID2gene, -3, 3, 1, paralogous)
    htmlLostTranscriptInGenesWSP2Plot, htmlLostTranscriptInGenesWSP2NormalizedPlot = plotter.makeLostTranscriptInGenesWSP2Plot(geneID2gene)
    htmlDifferencesInRelativeExpressionsBoxPlot = plotter.makeDifferencesInRelativeExpressionsBoxPlot(geneID2gene)

    #build the coverage scatterplots
    htmlScatterPlotCoverageOfMainIsoforms = plotter.makeScatterPlotCoverage(geneID2gene, "MainIsoforms")
    htmlScatterPlotCoverageOfGenes = plotter.makeScatterPlotCoverage(geneID2gene, "Genes")
    htmlScatterPlotCoverageOfIsoforms = plotter.makeScatterPlotCoverage(geneID2gene, "Isoforms")

    #Build the "Corretion collapsing paralogous gene families plots" plots
    htmlGeneralViewParalogFamilies, htmlScatterPlotSizeParalogFamilies = plotter.makeScatterPlotSizeParalogFamilies(paralogous)
    htmlGeneralViewParalogFamiliesExcluingUnchanged, htmlScatterPlotSizeParalogFamiliesExcluingUnchanged = plotter.makeScatterPlotSizeParalogFamilies(paralogous, True)
    htmlGeneralViewParalogFamiliesExcluingUnchangedCommonGenes, htmlScatterPlotSizeParalogFamiliesExcluingUnchangedCommonGenes = plotter.makeScatterPlotSizeParalogFamilies(paralogous, True, True)
    htmlParalogousGeneFamiliesSizeBarPlot = plotter.getParalogousGeneFamiliesSizeBarPlot(paralogous)

    #Build the splicing sites plots
    htmlCorrectIncorrectSSPlotScalar, htmlCorrectIncorrectSSPlotPercentage,\
    htmlDetailedIncorrectSSPlotScalar, htmlDetailedIncorrectSSPlotPercentage, \
    htmlSpliceSitesDistributionSSPlotScalar, htmlSpliceSitesDistributionSSPlotPercentage = \
        plotter.makeSpliceSitesPlots(splicingSitesProfiler)

    #Build the read connectivity plots
    htmlNbOfIdentifiedExonsLinePlotScalar, htmlNbOfIdentifiedExonsLinePlotPercentage = \
        plotter.buildNbOfIdentifiedExonsLinePlot(readSetProfiler, False), \
        plotter.buildNbOfIdentifiedExonsLinePlot(readSetProfiler, True)
    htmlHighestNbOfConsecutiveExonsLinePlotScalar, htmlHighestNbOfConsecutiveExonsLinePlotPercentage =\
        plotter.buildHighestNbOfConsecutiveExonsLinePlot(readSetProfiler, False), \
        plotter.buildHighestNbOfConsecutiveExonsLinePlot(readSetProfiler, True)
    htmlFullPartialReadsPlotScalar, htmlFullPartialReadsPlotPercentage = \
        plotter.buildFullPartialReadsPlot(readSetProfiler, False),\
        plotter.buildFullPartialReadsPlot(readSetProfiler, True)
    print "Computing the plots - Done!"



    # create the html report
    print "Creating HTML report..."

    def callFunctionAndPopulateTheReports(index, htmlTag, linesHTMLReport, linesHighResHTMLReport, object, methodName=None, **kwargs):
        if (linesHTMLReport!=None and htmlTag in linesHTMLReport[index]) or (linesHighResHTMLReport!=None and htmlTag in linesHighResHTMLReport[index]):
            if methodName!=None: #there is a method passed
                method = getattr(object, methodName)
                plots = method(**kwargs)
            else: #the method is the object itself
                plots = object

            if linesHTMLReport!=None:
                if type(plots) is dict and "imagePlot" in object:
                    linesHTMLReport[index] = linesHTMLReport[index].replace(htmlTag, str(plots["imagePlot"]))
                else:
                    linesHTMLReport[index] = linesHTMLReport[index].replace(htmlTag, str(plots))
            if linesHighResHTMLReport!=None:
                if type(plots) is dict and "jsPlot" in object:
                    linesHighResHTMLReport[index] = linesHighResHTMLReport[index].replace(htmlTag, str(plots["jsPlot"]))
                else:
                    linesHighResHTMLReport[index] = linesHighResHTMLReport[index].replace(htmlTag, str(plots))

    scriptDir = os.path.dirname(os.path.realpath(__file__))
    with open(scriptDir+"/lib/html/index_template.html") as indexTemplateFile:
        linesHTMLReport = indexTemplateFile.readlines()
        linesHighResHTMLReport = list(linesHTMLReport)

    for i in xrange(len(linesHTMLReport)):
        for feature in statProfiler.allFeatures:
            callFunctionAndPopulateTheReports(i, "<html_%s_Plot>"%feature, linesHTMLReport, linesHighResHTMLReport, \
                                              allStatsPlots[feature])

        callFunctionAndPopulateTheReports(i, "<statProfiler.getReadStatsAsJSArrayForHOT()>", linesHTMLReport, linesHighResHTMLReport, \
                                          statProfiler, "getReadStatsAsJSArrayForHOT")
        callFunctionAndPopulateTheReports(i, "<statProfiler.getBaseStatsAsJSArrayForHOT()>", linesHTMLReport, linesHighResHTMLReport, \
                                          statProfiler, "getBaseStatsAsJSArrayForHOT")
        callFunctionAndPopulateTheReports(i, "<statProfiler.getErrorStatsAsJSArrayForHOT()>", linesHTMLReport, linesHighResHTMLReport, \
                                          statProfiler, "getErrorStatsAsJSArrayForHOT")
        callFunctionAndPopulateTheReports(i, "<statProfiler.getAnnotationStatsAsJSArrayForHOT()>", linesHTMLReport, linesHighResHTMLReport, \
                                          statProfiler, "getAnnotationStatsAsJSArrayForHOT")
        callFunctionAndPopulateTheReports(i, "<geneProfiler.geneProfileToJSArrayForHOT()>", None, linesHighResHTMLReport, \
                                          geneProfiler, "geneProfileToJSArrayForHOT")
        callFunctionAndPopulateTheReports(i, "<geneProfiler.transcriptProfileToJSArrayForHOT()>", None, linesHighResHTMLReport, \
                                          geneProfiler, "transcriptProfileToJSArrayForHOT")
        callFunctionAndPopulateTheReports(i, "<htmlDifferenceOnTheNumberOfIsoformsPlotIntersection>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlDifferenceOnTheNumberOfIsoformsPlotIntersection)
        callFunctionAndPopulateTheReports(i, "<htmlDifferenceOnTheNumberOfIsoformsGeneFamilyPlotIntersection>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlDifferenceOnTheNumberOfIsoformsGeneFamilyPlotIntersection)
        callFunctionAndPopulateTheReports(i, "<htmlLostTranscriptInGenesWSP2Plot>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlLostTranscriptInGenesWSP2Plot)
        callFunctionAndPopulateTheReports(i, "<htmlLostTranscriptInGenesWSP2NormalizedPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlLostTranscriptInGenesWSP2NormalizedPlot)
        callFunctionAndPopulateTheReports(i, "<htmlDifferencesInRelativeExpressionsBoxPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlDifferencesInRelativeExpressionsBoxPlot)
        callFunctionAndPopulateTheReports(i, "<tools>", linesHTMLReport, linesHighResHTMLReport, \
                                          tools)
        callFunctionAndPopulateTheReports(i, "<htmlScatterPlotCoverageOfMainIsoforms>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlScatterPlotCoverageOfMainIsoforms)
        callFunctionAndPopulateTheReports(i, "<htmlScatterPlotCoverageOfGenes>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlScatterPlotCoverageOfGenes)
        callFunctionAndPopulateTheReports(i, "<htmlScatterPlotCoverageOfIsoforms>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlScatterPlotCoverageOfIsoforms)

        callFunctionAndPopulateTheReports(i, "<htmlTotalReadsCuttingPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          totalReadsCuttingPlot)
        callFunctionAndPopulateTheReports(i, "<htmlAlignedReadsCuttingPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          alignedReadsCuttingPlot)
        callFunctionAndPopulateTheReports(i, "<htmlUnalignedReadsCuttingPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          unalignedReadsCuttingPlot)

        callFunctionAndPopulateTheReports(i, "<htmlGeneralViewParalogFamilies>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlGeneralViewParalogFamilies)
        callFunctionAndPopulateTheReports(i, "<htmlScatterPlotSizeParalogFamilies>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlScatterPlotSizeParalogFamilies)
        callFunctionAndPopulateTheReports(i, "<htmlScatterPlotSizeParalogFamiliesExcluingUnchanged>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlScatterPlotSizeParalogFamiliesExcluingUnchanged)
        callFunctionAndPopulateTheReports(i, "<htmlScatterPlotSizeParalogFamiliesExcluingUnchangedCommonGenes>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlScatterPlotSizeParalogFamiliesExcluingUnchangedCommonGenes)
        callFunctionAndPopulateTheReports(i, "<htmlParalogousGeneFamiliesSizeBarPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlParalogousGeneFamiliesSizeBarPlot)
        callFunctionAndPopulateTheReports(i, "<htmlParalogousGeneFamiliesSizeBarPlot>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlParalogousGeneFamiliesSizeBarPlot)
        callFunctionAndPopulateTheReports(i, "<htmlCorrectIncorrectSSPlotScalar>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlCorrectIncorrectSSPlotScalar)
        callFunctionAndPopulateTheReports(i, "<htmlCorrectIncorrectSSPlotPercentage>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlCorrectIncorrectSSPlotPercentage)
        callFunctionAndPopulateTheReports(i, "<htmlDetailedIncorrectSSPlotScalar>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlDetailedIncorrectSSPlotScalar)
        callFunctionAndPopulateTheReports(i, "<htmlDetailedIncorrectSSPlotPercentage>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlDetailedIncorrectSSPlotPercentage)
        callFunctionAndPopulateTheReports(i, "<htmlSpliceSitesDistributionSSPlotScalar>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlSpliceSitesDistributionSSPlotScalar)
        callFunctionAndPopulateTheReports(i, "<htmlSpliceSitesDistributionSSPlotPercentage>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlSpliceSitesDistributionSSPlotPercentage)
        callFunctionAndPopulateTheReports(i, "<htmlNbOfIdentifiedExonsLinePlotScalar>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlNbOfIdentifiedExonsLinePlotScalar)
        callFunctionAndPopulateTheReports(i, "<htmlNbOfIdentifiedExonsLinePlotPercentage>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlNbOfIdentifiedExonsLinePlotPercentage)
        callFunctionAndPopulateTheReports(i, "<htmlHighestNbOfConsecutiveExonsLinePlotScalar>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlHighestNbOfConsecutiveExonsLinePlotScalar)
        callFunctionAndPopulateTheReports(i, "<htmlHighestNbOfConsecutiveExonsLinePlotPercentage>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlHighestNbOfConsecutiveExonsLinePlotPercentage)
        callFunctionAndPopulateTheReports(i, "<htmlFullPartialReadsPlotScalar>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlFullPartialReadsPlotScalar)
        callFunctionAndPopulateTheReports(i, "<htmlFullPartialReadsPlotPercentage>", linesHTMLReport, linesHighResHTMLReport, \
                                          htmlFullPartialReadsPlotPercentage)

        #fix the comment tags which will remove some stuff ot the simple report
        callFunctionAndPopulateTheReports(i, "<JSCommentOutSimpleReport>", linesHTMLReport, None, "/*")
        callFunctionAndPopulateTheReports(i, "<JSUncommentOutSimpleReport>", linesHTMLReport, None, "*/")
        callFunctionAndPopulateTheReports(i, "<JSCommentOutSimpleReport>", None, linesHighResHTMLReport, "")
        callFunctionAndPopulateTheReports(i, "<JSUncommentOutSimpleReport>", None, linesHighResHTMLReport, "")
        callFunctionAndPopulateTheReports(i, "<HTMLCommentOutSimpleReport>", linesHTMLReport, None, "<!--")
        callFunctionAndPopulateTheReports(i, "<HTMLUncommentOutSimpleReport>", linesHTMLReport, None, "-->")
        callFunctionAndPopulateTheReports(i, "<HTMLCommentOutSimpleReport>", None, linesHighResHTMLReport, "")
        callFunctionAndPopulateTheReports(i, "<HTMLUncommentOutSimpleReport>", None, linesHighResHTMLReport, "")




    #save the html reports
    with open(args.output+"/report_simple.html", "w") as indexOutFile:
        for line in linesHTMLReport:
            indexOutFile.write(line)
    with open(args.output+"/report_interactive_full.html", "w") as indexOutFile:
        for line in linesHighResHTMLReport:
            indexOutFile.write(line)

    #copy lib to the output
    if os.path.exists(args.output+"/lib"):
        shutil.rmtree(args.output+"/lib")
    shutil.copytree(scriptDir+"/lib", args.output+"/lib")

    # save paralogous gene to /lib/data/, if it exists
    if paralogous != None:
        with open(args.output + "/lib/data/paralogous_families.txt", "w") as paralogousFamiliesFile:
            paralogousFamiliesFile.write(str(paralogous))
    print "Creating HTML report... - Done!"

    print "We are finished!"


if __name__=="__main__":
    main()