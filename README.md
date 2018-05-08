# About LR_EC_analyser
LR_EC_analyser stands for Long Read Error Correction analyser. It is a python script that analyses the output of
long reads error correctors, like daccord, LoRDEC, LoRMA, MECAT, NaS, PBcR, pbdagcon, proovreads, etc. It does so by
running AlignQC (https://github.com/jason-weirather/AlignQC) on the BAMs built by the mapping the output of error correction
tools to a reference genome and parsing its output, and creating other custom plots, and then putting all the relevant information
in a HTML report. It also makes use of IGV.js (https://github.com/igvteam/igv.js) for an in-depth gene and transcript analysis.

# WARNING
Use Firefox or Safari to view the reports for the moment.

Google Chrome browser has a bug with IGV.js for now.

# Installing

For installing and running the script, please have pre-installed
```
python 2.7
virtualenv
samtools (v1.8+ http://www.htslib.org/download/)
R
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
    Firefox or Safari (bugged on Chrome for the moment)
## Viewing IGV plots (Gene stats and Transcript stats):
In order to serve the bams for the browser, you need to start a web server to serve the required files. To do so, you first need to install the tool and then:
```
cd LR_EC_analyser
source venv/bin/activate
python serve_files.py <output_folder>
```

# Parameters
```
usage: run_LR_EC_analyser.py [-h] [--genome GENOME] [--gtf GTF]
                             [--paralogous PARALOGOUS] [--raw RAWBAM]
                             [-o OUTPUT] [-t THREADS] [--skip_bam_process]
                             [--skip_alignqc] [--skip_copying]
                             <file.bam> [<file.bam> ...]

Long read error corrector analyser.

positional arguments:
  <file.bam>            BAM files of the Fastas output by the correctors

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       The genome in fasta file
  --gtf GTF             The transcriptome as GTF file
  --paralogous PARALOGOUS
                        Path to a file where the first two collumns denote
                        paralogous genes (see file GettingParalogs.txt to know
                        how you can get this file)
  --raw RAWBAM          The BAM file of the raw reads (i.e. the uncorrected
                        long reads file)
  -o OUTPUT             output folder
  -t THREADS            Number of threads to use
  --skip_bam_process    Skips bam processing - assume we had already done
                        this.
  --skip_alignqc        Skips AlignQC calls - assume we had already done this.
  --skip_copying        Skips copying genome and transcriptome to the output
                        folder.

```

# Thirdparties
1. AlignQC (https://github.com/jason-weirather/AlignQC)
2. IGV.js (https://github.com/igvteam/igv.js)
3. samtools (http://www.htslib.org)
4. handsontable (https://handsontable.com/)
5. matplotlib (https://matplotlib.org/)
6. mpld3 (http://mpld3.github.io/)
7. Twisted (https://twistedmatrix.com/trac/)
8. querystring (https://github.com/jgallen23/querystring)
9. jQuery (http://jquery.com/)
