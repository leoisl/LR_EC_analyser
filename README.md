# NEWS

v0.0.12 released. Changes:
* Added `--pdf` parameter;
* Fixed several bugs;

v0.0.11 released. Changes:
* Added `--colours` parameter;
* Fixed a bug in the interactive plot.

v0.0.10 released. Changes:
* Added long reads connectivity plots;
* Created a simple report which has plots as static images and has no Gene/Transcript read alignment viewer, making it easily shareable;

Software currently under development and testing. The latest version, found in this main branch, should be stable and working, but some bugs might be expected.

If you find one, please create an issue in the [https://gitlab.com/leoisl/LR_EC_analyser/issues](Issues page).

A pre-print with application of this method can be found at https://www.biorxiv.org/content/early/2018/11/23/476622 .

# About LR_EC_analyser
LR_EC_analyser stands for Long Read Error Correction analyser. It is a python script that analyses the output of
long reads error correctors, like LoRDEC, NaS, PBcR, proovread, canu, daccord, LoRMA, MECAT, pbdagcon, etc. It does so by
running AlignQC (https://github.com/jason-weirather/AlignQC) on the BAMs built by the mapping the output of error correction
tools to a reference genome (using for example gmap or minimap2) and parsing its output, and creating other custom plots, and then putting all the relevant information
in a HTML report. It also makes use of IGV.js (https://github.com/igvteam/igv.js) for an in-depth gene and transcript analysis.

# Motivation
Long-read sequencing technologies offer promising alternatives to high-throughput short read sequencing,
especially in the context of RNA-sequencing.
However these technologies are currently hindered by high error rates that affect analyses such as the identification of isoforms,
exon boundaries, open reading frames, and the creation of gene catalogues. Due to the novelty of such data,
computational methods are still actively being developed and options for the error-correction of RNA-sequencing long reads remain limited.
LR_EC_analyser can be applied to evaluate the extent to which existing long-read DNA error correction methods are capable of correcting
long reads.
It not only reports classical error-correction metrics but also the effect of correction on long read connectivity (impacts the inference of transcript structure and exon coupling),
gene families, isoform diversity, bias toward the major isoform, and splice site detection.
An application of this method to cDNA Nanopore reads can be found at https://www.biorxiv.org/content/early/2018/11/23/476622 .

# Example on real datasets
To know if LR_EC_analyser is useful for you, you can check the latest tool reports here: http://leoisl.gitlab.io/LR_EC_analyser_support .

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
usage: run_LR_EC_analyser.py [-h] --raw RAWBAM
                             [--self <self.bam> [<self.bam> ...]]
                             [--hybrid <hybrid.bam> [<hybrid.bam> ...]]
                             --genome GENOME --gtf GTF
                             [--paralogous PARALOGOUS] [-o OUTPUT]
                             [-t THREADS]
                             [--colours <self.colours> [<self.colours> ...]]
                             [--pdf] [--skip_bam_process] [--skip_alignqc]
                             [--skip_copying]

Long reads error corrector analyser.

optional arguments:
  -h, --help            show this help message and exit
  --raw RAWBAM          The BAM file of the raw reads (i.e. the uncorrected
                        long reads) mapped to the genome (preferably using
                        gmap -n 10 -f samse).
  --self <self.bam> [<self.bam> ...]
                        BAM files of the reads output by the SELF correctors
                        mapped to the genome (preferably using gmap -n 10 -f
                        samse).
  --hybrid <hybrid.bam> [<hybrid.bam> ...]
                        BAM files of the reads output by the HYBRID correctors
                        mapped to the genome (preferably using gmap -n 10 -f
                        samse).
  --genome GENOME       The genome in Fasta file format.
  --gtf GTF             The transcriptome in GTF file format.
  --paralogous PARALOGOUS
                        A file where the first two columns denote paralogous
                        genes (see file GettingParalogs.txt to know how you
                        can get this file).
  -o OUTPUT             Output folder
  -t THREADS            Number of threads to use
  --colours <self.colours> [<self.colours> ...]
                        A list of colours in hex encoding to use in the plots.
                        Colour shading is nice to show related corrections
                        (e.g. full-length, trimmed and split outputs from a
                        same tool),but unless the analysis is on few tools, it
                        is hard to find a nice automated choice of colour
                        shading. Hand-picking is more laborious but produces
                        better results.This parameter allows you to control
                        the colors of each tool. The order of the tools are:
                        raw -> hybrid -> self.The hybrid and self ordering are
                        given by parameter --hybrid and --self.See an example
                        of this parameter in https://gitlab.com/leoisl/LR_EC_a
                        nalyser/blob/master/scripts/command_line_paper.sh .
  --pdf                 Produce .pdf files of the plots in the
                        <output_folder>/plots directory.
  --skip_bam_process    Skips BAM processing (i.e. sorting and indexing BAM
                        files) - assume we had already done this.
  --skip_alignqc        Skips AlignQC calls - assume we had already done this.
  --skip_copying        Skips copying genome and transcriptome to the output
                        folder - assume we had already done this.
```

# Thirdparties
1. AlignQC (https://github.com/jason-weirather/AlignQC)
2. IGV.js (https://github.com/igvteam/igv.js)
3. Plotly (https://plot.ly/python/)
4. IGVTools (https://software.broadinstitute.org/software/igv/igvtools)
5. samtools (http://www.htslib.org)
6. handsontable (https://handsontable.com/)
7. matplotlib (https://matplotlib.org/) and mpld3 (http://mpld3.github.io/)
8. node.js (https://nodejs.org/)
9. npm (https://www.npmjs.com/)
10. http-server (https://www.npmjs.com/package/http-server)
11. nodeenv (https://github.com/ekalinin/nodeenv)
12. querystring (https://github.com/jgallen23/querystring)
13. jQuery (http://jquery.com/)
14. iziModal (http://izimodal.marcelodolce.com/)
15. jQuery BlockUI (http://malsup.com/jquery/block/)
16. anchorific.js (https://github.com/renettarenula/anchorific.js/)
17. psutil (https://github.com/giampaolo/psutil)
18. orca (https://github.com/plotly/orca)