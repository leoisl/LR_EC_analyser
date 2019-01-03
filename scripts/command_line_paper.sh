#1D cDNA command-line
source venv/bin/activate
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.gtf --paralogous /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.paralog_genes.80_percent_plus_identity.txt -o output -t 4 --raw /data2/ASTER/error_correction/raw_1D_only_data_no_mito_no_rDNA/gmap_bam/raw.bam --self /data2/ASTER/error_correction/canu/gmap_bam/canu.bam /data2/ASTER/error_correction/daccord/gmap_bam/daccord.bam /data2/ASTER/error_correction/daccord_trimmed/gmap_bam/daccord_trimmed.bam /data2/ASTER/error_correction/LoRMA/gmap_bam/LoRMA.bam /data2/ASTER/error_correction/MECAT/gmap_bam/MECAT.bam /data2/ASTER/error_correction/pbdagcon/gmap_bam/pbdagcon.bam /data2/ASTER/error_correction/pbdagcon_trimmed/gmap_bam/pbdagcon_trimmed.bam --hybrid /data2/ASTER/error_correction/HALC/gmap_bam/HALC.bam /data2/ASTER/error_correction/HALC_trimmed/gmap_bam/HALC_trimmed.bam /data2/ASTER/error_correction/HALC_split/gmap_bam/HALC_split.bam /data2/ASTER/error_correction/LoRDEC/gmap_bam/LoRDEC.bam /data2/ASTER/error_correction/LoRDEC_trimmed/gmap_bam/LoRDEC_trimmed.bam /data2/ASTER/error_correction/LoRDEC_split/gmap_bam/LoRDEC_split.bam /data2/ASTER/error_correction/LSC/gmap_bam/LSC.bam /data2/ASTER/error_correction/LSC_trimmed/gmap_bam/LSC_trimmed.bam /data2/ASTER/error_correction/NaS/gmap_bam/NaS.bam /data2/ASTER/error_correction/PBcR/gmap_bam/PBcR.bam /data2/ASTER/error_correction/proovread/gmap_bam/proovread.bam /data2/ASTER/error_correction/proovread_trimmed/gmap_bam/proovread_trimmed.bam --colours 000000 FF0000 F16666 FFCCCC 00009C 6666FF CCCCFF 006400 7FE37F FFFF00 00FFFF 654321 E4C2A0 FFA500 FF00FF FF7FFF BF00FF D68A59 696969 CFCFCF --skip_bam_process --skip_alignqc --skip_copying
deactivate


#direct RNA command-line
source venv/bin/activate
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/error_correction/Direct_RNA_human/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /data2/ASTER/error_correction/Direct_RNA_human/reference/Homo_sapiens.GRCh38.94.gtf --paralogous /data2/ASTER/error_correction/Direct_RNA_human/reference/Homo_sapiens.GRCh38.p12.94.paralog_genes.80_percent_plus_identity.txt -o output_RNA_direct -t 4 --raw /data2/ASTER/error_correction/Direct_RNA_human/raw.bam --self /data2/ASTER/error_correction/Direct_RNA_human/daccord/daccord.bam /data2/ASTER/error_correction/Direct_RNA_human/daccord_trimmed/daccord_trimmed.bam /data2/ASTER/error_correction/Direct_RNA_human/LoRMA/LoRMA.bam /data2/ASTER/error_correction/Direct_RNA_human/MECAT/MECAT.bam /data2/ASTER/error_correction/Direct_RNA_human/pbdagcon/pbdagcon.bam /data2/ASTER/error_correction/Direct_RNA_human/pbdagcon_trimmed/pbdagcon_trimmed.bam
#if you need just to recompute the plots, use the --skips:
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/error_correction/Direct_RNA_human/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /data2/ASTER/error_correction/Direct_RNA_human/reference/Homo_sapiens.GRCh38.94.gtf --paralogous /data2/ASTER/error_correction/Direct_RNA_human/reference/Homo_sapiens.GRCh38.p12.94.paralog_genes.80_percent_plus_identity.txt -o output_RNA_direct -t 4 --raw /data2/ASTER/error_correction/Direct_RNA_human/raw.bam --self /data2/ASTER/error_correction/Direct_RNA_human/daccord/daccord.bam /data2/ASTER/error_correction/Direct_RNA_human/daccord_trimmed/daccord_trimmed.bam /data2/ASTER/error_correction/Direct_RNA_human/LoRMA/LoRMA.bam /data2/ASTER/error_correction/Direct_RNA_human/MECAT/MECAT.bam /data2/ASTER/error_correction/Direct_RNA_human/pbdagcon/pbdagcon.bam /data2/ASTER/error_correction/Direct_RNA_human/pbdagcon_trimmed/pbdagcon_trimmed.bam --skip_bam_process --skip_alignqc --skip_copying
deactivate



#about colors:
raw:
black ("#000000") - raw

#trimmed and split:
red ("#FF0000", "#F16666", "#FFCCCC") - HALC
blue ("#00009C", "#6666FF", "#CCCCFF") lordec

#trimmed
green ("#006400", "#7FE37F") lsc
brown ("#654321", "#E4C2A0") proovread
magenta ("#FF00FF", "#FF7FFF") daccord
gray ("#696969", "#CFCFCF") pbdagcon

#untrimmed (5):
yellow ("#FFFF00") nas
cyan ("#00FFFF") pbcr
orange("#FFA500") canu
purple ("#BF00FF") lorma
? ("#D68A59") mecat



trimming/splitting details:
details:
LORMA IS IN TRIMMED + SPLIT BECAUSE:
	"The script lorma.sh runs LoRDEC iteratively, trims and splits the resulting initial corrected reads, and finally runs LoRMA on these trimmed and split reads.", from https://www.cs.helsinki.fi/u/lmsalmel/LoRMA/README.txt

PBCR = WHERE TO PUT???? ("The consensus module may trim/split/discard sequences if there is insufficient coverage to achieve a minimum QV for the corrected sequence.", from http://wgs-assembler.sourceforge.net/wiki/index.php/PBcR)
	FOR NOW, IN THE UNTRIMMED_UNSPLIT