samtools view -b -F 0x900 ${1} > ${1}.best_hits.bam
samtools view -H ${1} > ${1}.header
samtools view ${1}.best_hits.bam > ${1}.best_hits.sam
cat ${1}.header ${1}.best_hits.sam > ${1}.best_hits.sam2
rm ${1}.header ${1}.best_hits.sam
mv ${1}.best_hits.sam2 ${1}.best_hits.sam
samtools view -b -S ${1}.best_hits.sam > ${1}.best_hits.bam
rm ${1}.best_hits.sam