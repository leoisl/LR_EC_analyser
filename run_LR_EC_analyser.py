#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import shutil
from Profiler import *
from ExternalTools import *
from Parsers import *
from Plotter import *
from Paralogous import *
import sys


"""
Testing:
python run_LR_EC_analyser.py --genome sample_data/Mus_musculus.GRCm38.dna.chromosome.19.fa --gtf sample_data/Mus_musculus.GRCm38.91.chr19.gtf -t 4 -o sample_data/output --raw sample_data/gmap_CB_1Donly_to_GRCm38_chr19.bam sample_data/good.gmap.chr19.bam sample_data/indels.gmap.chr19.bam sample_data/subs.gmap.chr19.bam



#How we built the index
gmap_build -D /data2/ASTER/nanopore_mouse/Our_mapping_to_mRNA_only/gmap_genome -d GRCm38 /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa

#mapping
GMAP_GENOME_DIR=/data2/ASTER/nanopore_mouse/Our_mapping_to_mRNA_only/gmap_genome
GMAP_GENOME_NAME=GRCm38
gmap -D $GMAP_GENOME_DIR -d $GMAP_GENOME_NAME -n 10 -t $THREADS -f samse $RAW_READS_FILE >  ${OUTPUT_PREFIX}.gmap.sam

#from sam to bam
samtools view -S -o ${OUTPUT_PREFIX}.gmap.bam ${OUTPUT_PREFIX}.gmap.sam
rm ${OUTPUT_PREFIX}.gmap.sam







Input:
python run_LR_EC_analyser.py -g <genome> -gtf <transcriptome> raw_reads.bam LoRDEC.bam LoRMA.bam PBCR.bam ...

MappingInfo is just a class with relevant mapping info to compare 2 EC tools
For now, I think it should have:
    1. # of reads mapped to transcript/gene
    2. overlap size between read and transcript (all the values / average)?
    

0/ Process the gtf and get the gene coordinates
    
1/ For each .bam file:
    1/ Filter the bam to keep only the best hit
    2/ Run AlignQC on this bam
    3/ Process output/data/annotbest.txt and populate gene2ECTool2MappingInfo and transcript2ECTool2MappingInfo
    4/ Create a discrepancy measure (I think now it should be like the largest difference between the # of reads in the raw_reads and one of the tools)
    5/ Create an HTML page sorted by this discrepancy measure and show all the infos we have gathered
    6/ On clicking in one of the genes, we give to igv viewer the genome, transcriptome, all bams, and the gene coordinate



To view results, execute
python run_LR_EC_analyser --view <path_to_results>
    This will start a ftp server on <path_to_results>
    Then will open a browser to view <path_to_results>



"""

def main():
    parser = argparse.ArgumentParser(description='Long read error corrector analyser.')
    parser.add_argument('bams', metavar='<file.bam>', type=str, nargs='+',
                        help='BAM files of the Fastas output by the correctors')
    parser.add_argument("--genome", dest="genome", help="The genome in fasta file")
    parser.add_argument("--gtf", dest="gtf", help="The transcriptome as GTF file")
    parser.add_argument("--paralogous", help="Path to a file where the first two collumns denote paralogous genes (see file GettingParalogs.txt to know how you can get this file)")
    parser.add_argument("--raw", dest="rawBam", help="The BAM file of the raw reads (i.e. the uncorrected long reads file)")
    parser.add_argument("-o", dest="output", help="output folder", default="output/")
    parser.add_argument("-t", dest="threads", type=int, help="Number of threads to use")
    parser.add_argument("--skip_bam_process", dest="skip_bam", action="store_true", help="Skips bam processing - assume we had already done this.")
    parser.add_argument("--skip_alignqc", dest="skip_alignqc", action="store_true",
                        help="Skips AlignQC calls - assume we had already done this.")
    parser.add_argument("--skip_copying", dest="skip_copying", action="store_true",
                        help="Skips copying genome and transcriptome to the output folder.")
    args=parser.parse_args()

    #create output dir
    if not os.path.exists(args.output):
        os.makedirs(args.output)


    #get some useful vars
    bams = [args.rawBam] + args.bams
    tools = ["raw.bam"] + [os.path.basename(bam) for bam in args.bams]

    #copy genome to output and index it
    genome = args.output + "/" + os.path.basename(args.genome)
    if not args.skip_copying:
        shutil.copy(args.genome, args.output)
        indexGenome(genome)
    else:
        print "Skipping genome copying and indexing..."

    #get genes from gtf
    geneID2gene = parseGTFToGetGenes(args.gtf, tools)

    # copy gtf to output
    gtf = args.output + "/" + os.path.basename(args.gtf)
    if not args.skip_copying:
        shutil.copy(args.gtf, args.output)
    else:
        print "Skipping transcriptome copying..."

    #get the paralogous info, if given
    paralogous = None
    if args.paralogous != None:
        paralogous = Paralogous(geneID2gene)
        paralogous.readParalogousFile(args.paralogous)


    #TODO: walk the bam with pysam and get the read lengths to compute mean length of the aligned reads

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
        print "Skipping AlignQC runs../"

    #create the output by parsing AlignQC results
    print "Running Stat profiler..."
    statProfiler = StatProfiler(tools, args.output)
    print "Running Stat profiler - Done!"

    #create the gene profiler
    print "Running Gene profiler..."
    tool2Bam = {tool: os.path.basename(bam) for tool, bam in zip(tools, sortedBams)}
    geneProfiler = FeatureProfiler(geneID2gene, tools, tool2Bam, os.path.basename(genome), os.path.basename(gtf), args.output)
    #populate the gene profiler
    for tool in tools:
        geneProfiler.populateFromAnnotbest(tool, args.output)
    print "Running Gene profiler - Done!"


    #create the Plotter and the plots
    print "Computing the plots..."
    plotsOutput = args.output+"/plots"
    if not os.path.exists(plotsOutput):
        os.makedirs(plotsOutput)
    plotter = Plotter(tools, plotsOutput)
    htmlDifferenceOnTheNumberOfIsoformsPlot = plotter.makeDifferenceOnTheNumberOfIsoformsPlot(geneID2gene, -3, 3)
    htmlLostTranscriptInGenesWSP2Plot = plotter.makeLostTranscriptInGenesWSP2Plot(geneID2gene)
    htmlDifferencesInRelativeExpressionsBoxPlot = plotter.makeDifferencesInRelativeExpressionsBoxPlot(geneID2gene)
    if paralogous != None:
        htmlScatterPlotSizeParalogFamilies = plotter.makeScatterPlotSizeParalogFamilies(geneID2gene, paralogous)
        htmlScatterPlotSizeParalogFamiliesExcluingUnchanged = plotter.makeScatterPlotSizeParalogFamilies(geneID2gene, paralogous, True)
    print "Computing the plots - Done!"



    # create the html report
    print "Creating HTML report..."
    with open("lib/html/index_template.html") as indexTemplateFile:
        indexTemplateLines = indexTemplateFile.readlines()

        for i, line in enumerate(indexTemplateLines):
            if "<statProfiler.getReadStatsAsJSArrayForHOT()>" in line:
                line = line.replace("<statProfiler.getReadStatsAsJSArrayForHOT()>", statProfiler.getReadStatsAsJSArrayForHOT())
            if "<statProfiler.getBaseStatsAsJSArrayForHOT()>" in line:
                line = line.replace("<statProfiler.getBaseStatsAsJSArrayForHOT()>", statProfiler.getBaseStatsAsJSArrayForHOT())
            if "<statProfiler.getErrorStatsAsJSArrayForHOT()>" in line:
                line = line.replace("<statProfiler.getErrorStatsAsJSArrayForHOT()>", statProfiler.getErrorStatsAsJSArrayForHOT())
            if "<geneProfiler.geneProfileToJSArrayForHOT()>" in line:
                line = line.replace("<geneProfiler.geneProfileToJSArrayForHOT()>", geneProfiler.geneProfileToJSArrayForHOT())
            if "<geneProfiler.transcriptProfileToJSArrayForHOT()>" in line:
                line = line.replace("<geneProfiler.transcriptProfileToJSArrayForHOT()>", geneProfiler.transcriptProfileToJSArrayForHOT())
            if "<tools>" in line:
                line = line.replace("<tools>", str(tools))
            if "<htmlDifferenceOnTheNumberOfIsoformsPlot>" in line:
                line = line.replace("<htmlDifferenceOnTheNumberOfIsoformsPlot>", htmlDifferenceOnTheNumberOfIsoformsPlot)
            if "<htmlLostTranscriptInGenesWSP2Plot>" in line:
                line = line.replace("<htmlLostTranscriptInGenesWSP2Plot>", htmlLostTranscriptInGenesWSP2Plot)
            if "<htmlDifferencesInRelativeExpressionsBoxPlot>" in line:
                line = line.replace("<htmlDifferencesInRelativeExpressionsBoxPlot>", htmlDifferencesInRelativeExpressionsBoxPlot)
            if "<htmlScatterPlotSizeParalogFamilies>" in line:
                if paralogous != None:
                    line = line.replace("<htmlScatterPlotSizeParalogFamilies>", htmlScatterPlotSizeParalogFamilies)
                else:
                    line = line.replace("<htmlScatterPlotSizeParalogFamilies>", "<p style='color: red; font-size: large;'>Paralogous file (--paralogous parameter) was not given, so we did not produce this plot. </p>")
            if "<htmlScatterPlotSizeParalogFamiliesExcluingUnchanged>" in line:
                if paralogous != None:
                    line = line.replace("<htmlScatterPlotSizeParalogFamiliesExcluingUnchanged>", htmlScatterPlotSizeParalogFamiliesExcluingUnchanged)
                else:
                    line = line.replace("<htmlScatterPlotSizeParalogFamiliesExcluingUnchanged>", "<p style='color: red; font-size: large;'>Paralogous file (--paralogous parameter) was not given, so we did not produce this plot. </p>")
            indexTemplateLines[i] = line


    #save the html report
    with open(args.output+"/report.html", "w") as indexOutFile:
        for line in indexTemplateLines:
            indexOutFile.write(line)

    #copy lib to the output
    if os.path.exists(args.output+"/lib"):
        shutil.rmtree(args.output+"/lib")
    shutil.copytree("lib", args.output+"/lib")
    print "Creating HTML report... - Done!"



    print "We are finished!"


if __name__=="__main__":
    main()