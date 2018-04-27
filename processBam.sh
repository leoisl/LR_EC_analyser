#usage: bash getOnlyBestHitsFromBam.sh <bam_file> <threads>

#sort and index the bam
samtools sort -@ ${2} ${1}.bam -o ${1}.sorted.bam
samtools index ${1}.sorted.bam