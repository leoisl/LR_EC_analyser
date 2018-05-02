# About LR_EC_analyser
LR_EC_analyser stands for Long Read Error Correction analyser. It is a simple python script that analyses the output of
long reads error correctors, like daccord, LoRDEC, LoRMA, MECAT, NaS, PBcR, pbdagcon, proovreads, etc. It does so by
running AlignQC (https://github.com/jason-weirather/AlignQC) on the BAMs built by mapping the output of error correction
tools to a reference genome and parsing its output, then putting all the relevant information in a HTML report. It also
makes use of IGV.js (https://github.com/igvteam/igv.js) for an in-depth gene and transcript analysis.

# WARNING
Use Firefox or Safari to view the reports for the moment.

# Requirements

## For running the script, please have installed
```
python 2.7
virtualenv
samtools (v1.8+ http://www.htslib.org/download/)
R
```

## For viewing the results
### Viewing the html page:
    Firefox or Safari for the moment (bugged on Chrome for the moment)
### Viewing IGV plots (Gene stats and Transcript stats):
```
cd LR_EC_analyser
python serve_files.py --output <output_folder>
```


# Installing
```
git clone --recursive https://gitlab.com/leoisl/LR_EC_analyser
```

# Running on a sample example
This sample example contains some real long reads as raw reads, and subsets of the raw reads with some substitutions and
indels artificially inserted as error-corrected reads. It is just for a proof of concept, no real error correction tools
were ran in these reads. To run the tool on this sample example, do:

```
python run_LR_EC_analyser.py --genome sample_data/Mus_musculus.GRCm38.dna.chromosome.19.fa --gtf sample_data/Mus_musculus.GRCm38.91.chr19.gtf -t 4 -o sample_data/output --raw sample_data/gmap_CB_1Donly_to_GRCm38_chr19.bam sample_data/good.gmap.chr19.bam sample_data/indels.gmap.chr19.bam sample_data/subs.gmap.chr19.bam
```

The output can be found here: https://www.dropbox.com/s/phdsty29yszn47j/output_sample_data.tar.gz?dl=1

# Parameters
```
usage: run_LR_EC_analyser.py [-h] [--genome GENOME] [--gtf GTF] [--raw RAWBAM]
                             [-o OUTPUT] [-t THREADS] [--skip_bam_process]
                             [--skip_alignqc] [--skip_copying]
                             [--skip_splitting_bam]
                             <file.bam> [<file.bam> ...]

Long read error corrector analyser.

positional arguments:
  <file.bam>            BAM files of the Fastas output by the correctors

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       The genome in fasta file
  --gtf GTF             The transcriptome as GTF file
  --raw RAWBAM          The BAM file of the raw reads (i.e. the uncorrected
                        long reads file)
  -o OUTPUT             output folder
  -t THREADS            Number of threads to use
  --skip_bam_process    Skips bam processing - assume we had already done
                        this.
  --skip_alignqc        Skips AlignQC calls - assume we had already done this.
  --skip_copying        Skips copying genome and transcriptome to the output
                        folder.
  --skip_splitting_bam  Skips splitting the bams by genes and transcripts.
```

# Thirdparties
1. AlignQC (https://github.com/jason-weirather/AlignQC)
2. IGV.js (https://github.com/igvteam/igv.js)
3. samtools (http://www.htslib.org)
4. handsontable (https://handsontable.com/)
5.matplotlib (https://matplotlib.org/)
6. mpld3 (http://mpld3.github.io/)
7. Twisted (https://twistedmatrix.com/trac/)
8. querystring (https://github.com/jgallen23/querystring)
9. jQuery (http://jquery.com/)
