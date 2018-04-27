#usage: bash getOnlyBestHitsFromBam.sh <bam_file> <threads>

#put the best hits in ${1}.best_hits.bam
samtools view -b -F 0x900 -o ${1}.best_hits.bam ${1}

#sort and index the bam
samtools sort -@ ${2} ${1}.best_hits.bam -o ${1}.best_hits.sorted.bam
rm ${1}.best_hits.bam
samtools index ${1}.best_hits.sorted.bam