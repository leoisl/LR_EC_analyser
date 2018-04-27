#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess
import os
import gzip


"""
Testing:
python run_LR_EC_analyser.py --genome AlignQC/example_data/chr21chr22chrM.fa.gz --gtf AlignQC/example_data/chr21chr22chrM.gencode25.gtf.gz --raw AlignQC/example_data/testraw.bam AlignQC/example_data/test1.bam AlignQC/example_data/test2.bam



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

scriptDir = os.path.dirname(os.path.realpath(__file__))

class Gene:
    """
    Represents a gene
    """
    def __init__ (self, geneId, chromosome, begin, end, strand):
        self.geneId=geneId
        self.chromosome=chromosome
        self.begin=begin
        self.end=end
        self.strand=strand

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

#TODO: we should add more data here...
class Profile:
    """
    Represents the profile of all tools to a particular gene - nb of reads mapped to it, quality of the reads, etc...
    """
    def __init__(self, tools):
        # TODO: we should add more data here...
        self.tool2NbOfMappedReads = {tool: 0 for tool in tools}

    def computeLargestDiscrepancy(self):
        """
        :return: the largest value between raw_reads and a tool
        """
        nbOfMappedRawReads = self.tool2NbOfMappedReads["raw"]
        return max([abs(nbOfMappedRawReads-self.tool2NbOfMappedReads[tool]) for tool in self.tool2NbOfMappedReads] )


class GeneProfiler:
    """
    Represents all genes and their profiles
    """
    def __init__(self, genes, tools):
        self.geneId2Gene={gene.geneId:gene for gene in genes}
        self.geneId2Profile={gene.geneId:Profile(tools) for gene in genes}

    def populateFromAnnotbest(self, tool):
        """
        Reads the annotbest.txt.gz file from AlignQC and populates this profiler
        annotbest.txt is the AlignQC file that contains the best mapping of the ANNOTATED READS (reads that we were able to align and that AlignQC was able to assign to a transcript)
    annotbest.txt format:
    columns:
    0: id of the line (?)
    1: read name
    2: gene name
    3: transcript name
    4: partial or full (if it mapped partially or fully in the genome according to AlignQC)
    5: # ?
    6: # most consecutive exons in read
    7: # exons in read
    8: # exons in the transcript
    9: overlap size between read and transcript
    10: read length
    11: transcript length
    12: alignment coordinates
    13: transcript coordinates

    Example:
    5	m150117_043342_42142_c100769800150000001823165407071580_s1_p0/144819/ccs	ENSG00000274276.4	ENST00000624934.3	partial	15	11	16	18	1747	3884	1992	chr21:6446736-6467516	chr21:6445433-6467532	2454
        """
        dataFolder = "alignqc_out_on_%s/data"%tool
        with gzip.open(dataFolder+"/annotbest.txt.gz") as file:
            for line in file:
                lineSplit = line.rstrip().split()
                self.geneId2Profile[lineSplit[2]].tool2NbOfMappedReads[tool] += 1


def parseGTFToGetGenes(gtf):
    print "Parsing %s..."%gtf
    genes=[]
    with gzip.open(gtf) as gtfFile:
        for line in gtfFile:
            if line[0]!="#":
                lineSplit = line.rstrip().split()
                if lineSplit[2]=="gene":
                    genes.append(Gene(lineSplit[lineSplit.index("gene_id")+1][1:-2],
                        lineSplit[0], int(lineSplit[3]), int(lineSplit[4]), lineSplit[6]))
    print "Parsing %s - Done!" % gtf
    return genes

def processBam(bam, threads):
    print "Sorting and indexing %s..." % bam
    commandLine = "bash %s/processBam.sh %s %d"%(scriptDir, bam, threads)
    print "Running %s"%commandLine
    subprocess.check_call(commandLine.split())
    print "Sorting and indexing %s - Done!" % bam


def runAlignQC(tool, bam, genome, gtf, threads):
    print "Running AlignQC for %s..." % bam
    commandLine = "alignqc analyze %s -g %s -t %s --output_folder alignqc_out_on_%s --threads %d" % \
                            (bam, genome, gtf, tool, threads)
    print "Running %s" % commandLine
    subprocess.check_call(commandLine.split())
    print "Running AlignQC for %s - Done!" % bam

def readFileComposedOfPairStringIntToDict(filename):
    stats = {}
    with open(filename) as file:
        for line in file:
            lineSplit = line.split()
            stats[lineSplit[0]] = int(lineSplit[1])
    return stats

def processRarefractionFile(filename, totalReads):
    """
    Goes through the rarefraction file and get the values for rarefraction in the last line
    :param filename:
    :return:
    """
    with open(filename) as file:
        lastLineSplitted = file.readlines()[-1].split()

    nbOfReads = int(lastLineSplitted[0])
    rarefractionLow = int(float(lastLineSplitted[1]))
    rarefractionMed = int(float(lastLineSplitted[2]))
    rarefractionHigh = int(float(lastLineSplitted[3]))

    if rarefractionLow==rarefractionMed and rarefractionLow==rarefractionHigh and nbOfReads==totalReads:
        return rarefractionLow
    else:
        raise Exception("Rarefraction processing error!")


def parseAlignQCOutput(tool):
    dataFolder = "alignqc_out_on_%s/data"%tool
    statsDic = {}

    '''
    Added the following fields:
    TOTAL_READS	1282
    UNALIGNED_READS	200
    ALIGNED_READS	1082
    SINGLE_ALIGN_READS	996
    GAPPED_ALIGN_READS	73
    CHIMERA_ALIGN_READS	13
    TRANSCHIMERA_ALIGN_READS	0
    SELFCHIMERA_ALIGN_READS	13
    TOTAL_BASES	4561236
    UNALIGNED_BASES	2257927
    ALIGNED_BASES	2303309
    SINGLE_ALIGN_BASES	2222158
    GAPPED_ALIGN_BASES	58561
    CHIMERA_ALIGN_BASES	22590
    TRANSCHIMERA_ALIGN_BASES	0
    SELFCHIMERA_ALIGN_BASES	22590
    '''
    statsDic.update(readFileComposedOfPairStringIntToDict(dataFolder+"/alignment_stats.txt"))

    '''
    Added the following fields:
    ALIGNMENT_COUNT	501
    ALIGNMENT_BASES	1097266
    ANY_ERROR	118039
    MISMATCHES	40866
    ANY_DELETION	28083
    COMPLETE_DELETION	13374
    HOMOPOLYMER_DELETION	14709
    ANY_INSERTION	49090
    COMPLETE_INSERTION	28668
    HOMOPOLYMER_INSERTION	20422
    '''
    statsDic.update(readFileComposedOfPairStringIntToDict(dataFolder + "/error_stats.txt"))


    statsDic["GENES_DETECTED_ANY_MATCH"] = processRarefractionFile(dataFolder+"/gene_rarefraction.txt", statsDic["TOTAL_READS"])
    statsDic["GENES_DETECTED_FULL_MATCH"] = processRarefractionFile(dataFolder + "/gene_full_rarefraction.txt", statsDic["TOTAL_READS"])

    return statsDic


def main():
    parser = argparse.ArgumentParser(description='Long read error corrector analyser.')
    parser.add_argument('bams', metavar='<file.bam>', type=str, nargs='+',
                        help='BAM files of the Fastas output by the correctors')
    parser.add_argument("--genome", dest="genome", help="The genome in gzipped fasta file (this file needs to be a .gz)")
    parser.add_argument("--gtf", dest="gtf", help="The transcriptome as gzipped GTF file (this file needs to be a .gz)")
    parser.add_argument("--raw", dest="rawBam", help="The BAM file of the raw reads (i.e. the uncorrected long reads file)")
    parser.add_argument("-o", dest="output", help="output file name prefix", default="output")
    parser.add_argument("-t", dest="threads", type=int, help="Number of threads to use")
    args=parser.parse_args()

    bams = [args.rawBam] + args.bams

    #get the tools names
    tools = ["raw"] + [os.path.basename(bam) for bam in args.bams]

    #get genes from gtf
    genes = parseGTFToGetGenes(args.gtf)

    #create the gene profiler
    geneProfiler = GeneProfiler(genes, tools)

    #TODO: walk the bam with pysam and get the read lengths to compute mean length of the aligned reads

    #sort and index bam
    for bam in bams:
        processBam(bam, args.threads)

    #update the bam files
    sortedBams = [bamfile+".sorted.bam" for bamfile in bams]

    #run AlignQC on the bams
    for tool, bam in zip(tools, sortedBams):
        runAlignQC(tool, bam, args.genome, args.gtf, args.threads)

    #create the output by parsing AlignQC results
    tool2Stats={tool:parseAlignQCOutput(tool) for tool in tools}

    #populate the gene profiler
    for tool in tools:
        geneProfiler.populateFromAnnotbest(tool)

    #with open(args.output+".html", "w") as htmlOutput:









if __name__=="__main__":
    main()