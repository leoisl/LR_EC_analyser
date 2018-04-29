#usage: bash getOnlyBestHitsFromBam.sh <bam_file> <output_file> <threads>
set -eux

#sort and index the bam
samtools sort -@ $3 $1 -o $2
samtools index $2