# About LR_EC_analyser
LR_EC_analyser stands for Long Read Error Correction analyser. It is a python script that analyses the output of
long reads error correctors, like daccord, LoRDEC, LoRMA, MECAT, NaS, PBcR, pbdagcon, proovreads, etc. It does so by
running AlignQC (https://github.com/jason-weirather/AlignQC) on the BAMs built by the mapping the output of error correction
tools to a reference genome and parsing its output, and creating other custom plots, and then putting all the relevant information
in a HTML report. It also makes use of IGV.js (https://github.com/igvteam/igv.js) for an in-depth gene and transcript analysis.

# Motivation
RNA long reads from third generation sequencers like PacBio or ONT are becoming widely used due to the fact that such reads are usually
long enough to cover entire transcripts, removing the need for transcriptome assembly and allowing for the discovery of the full
transcript structure and long-range exon coupling. One of the main issues on using long reads is their high error rates. As such,
many error correctors were developed in the past few years, that are broadly classified into two categories: non-hybrid (correct long reads
using only long reads) and hybrid (correct long reads using short reads). However, the big majority of these error correctors are tailored for
the genomic context, and transcriptomic characteristics like the presence of alternative splicing, which can create long gaps between common sequences,
the highly dynamic isoform expression levels, which create regions with high coverage and low coverage, can pose difficulties to these tools, since
they are not present in the DNA level, and thus probably not modelled. Ultimately, it is not simple to choose which error corrector to use, mainly if
all your options are tailored for correcting genomic long reads. This tool aims at filling this gap, by evaluating the quality of correction of different
error correctors when applied to a RNA long reads dataset. This tool is not suited to evaluate how well the tools perform on correction DNA long reads.
It can help users choose which error corrector to use on a specific transcritome dataset.

# Installing

For installing and running the script, please have pre-installed
```
python 2.7
virtualenv
samtools (v1.8+ http://www.htslib.org/download/)
R
Java 8
```

Then you need to download the source files and setup the virtual environment. This can be done with:
```
git clone https://gitlab.com/leoisl/LR_EC_analyser
cd LR_EC_analyser
bash setup.sh
```


# Running on a sample example

TODO

# For viewing the results
## Viewing the html page:
    Any browser (Google Chrome is recommended, since it seems to be the fastest to view IGV plots)
## Viewing IGV plots (Gene stats and Transcript stats):
In order to serve the bams for the browser, you need to start a web server to serve the required files. To do so, you first need to install the tool (see [Installing section](#installing)) and then:
```
cd LR_EC_analyser
source venv/bin/activate
python serve_files.py <output_folder>
```

# Parameters
```
usage: run_LR_EC_analyser.py [-h] --genome GENOME --gtf GTF
                             [--paralogous PARALOGOUS] --raw RAWBAM
                             [-o OUTPUT] [-t THREADS] [--skip_bam_process]
                             [--skip_alignqc] [--skip_copying]
                             <file.bam> [<file.bam> ...]

Long reads error corrector analyser.

positional arguments:
  <file.bam>            BAM files of the reads output by the correctors mapped
                        to the genome (preferably using gmap -n 10 -f samse).

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       The genome in Fasta file format.
  --gtf GTF             The transcriptome in GTF file format.
  --paralogous PARALOGOUS
                        A file where the first two columns denote paralogous
                        genes (see file GettingParalogs.txt to know how you
                        can get this file).
  --raw RAWBAM          The BAM file of the raw reads (i.e. the uncorrected
                        long reads) mapped to the genome.
  -o OUTPUT             Output folder
  -t THREADS            Number of threads to use
  --skip_bam_process    Skips BAM processing (i.e. sorting and indexing BAM
                        files) - assume we had already done this.
  --skip_alignqc        Skips AlignQC calls - assume we had already done this.
  --skip_copying        Skips copying genome and transcriptome to the output
                        folder - assume we had already done this.
```

# Thirdparties
1. AlignQC (https://github.com/jason-weirather/AlignQC)
2. IGV.js (https://github.com/igvteam/igv.js)
3. IGVTools (https://software.broadinstitute.org/software/igv/igvtools)
4. samtools (http://www.htslib.org)
5. handsontable (https://handsontable.com/)
6. matplotlib (https://matplotlib.org/) and mpld3 (http://mpld3.github.io/)
7. node.js (https://nodejs.org/)
8. npm (https://www.npmjs.com/)
9. http-server (https://www.npmjs.com/package/http-server)
10. nodeenv (https://github.com/ekalinin/nodeenv)
11. querystring (https://github.com/jgallen23/querystring)
12. jQuery (http://jquery.com/)
