# About LR_EC_analyser
LR_EC_analyser stands for Long Read Error Correction analyser. It is a simple python script that analyses the output of
long reads error correctors, like daccord, LoRDEC, LoRMA, MECAT, NaS, PBcR, pbdagcon, proovreads, etc. It does so by
running AlignQC (https://github.com/jason-weirather/AlignQC) on the BAMs built by mapping the output of error correction
tools to a reference genome and parsing its output, then putting all the relevant information in a HTML report. It also
makes use of IGV.js (https://github.com/igvteam/igv.js) for an in-depth gene and transcript analysis.

# Requirements

## For running the script
```
python 2.7
virtualenv
samtools
```

## For viewing the results
### Viewing the html page:
    Any browser
### Viewing IGV plots (Gene stats):
    1/ Install https://www.npmjs.com/package/http-server
    2/ Go to the output folder
    3/ http-server --cors -p 19974


# Installing
```
git clone --recursive https://gitlab.com/leoisl/LR_EC_analyser
```

# Running on a sample example
TODO (virtual env and sample example)

# Parameters
```
usage: run_LR_EC_analyser.py [-h] [--genome GENOME] [--gtf GTF] [--raw RAWBAM]
                             [-o OUTPUT] [-t THREADS] [--skip_bam_process]
                             [--skip_alignqc]
                             <file.bam> [<file.bam> ...]

Long read error corrector analyser.

positional arguments:
  <file.bam>          BAM files of the Fastas output by the correctors

optional arguments:
  -h, --help          show this help message and exit
  --genome GENOME     The genome in fasta file
  --gtf GTF           The transcriptome as GTF file
  --raw RAWBAM        The BAM file of the raw reads (i.e. the uncorrected long
                      reads file)
  -o OUTPUT           output folder
  -t THREADS          Number of threads to use
  --skip_bam_process  Skips bam processing - assume we had already done this.
  --skip_alignqc      Skips AlignQC calls - assume we had already done this.
```

# Thirdparties
1. AlignQC (https://github.com/jason-weirather/AlignQC)
2. IGV.js (https://github.com/igvteam/igv.js)
3. handsontable (https://handsontable.com/)
4. jQuery (http://jquery.com/)
5. querystring (https://github.com/jgallen23/querystring)