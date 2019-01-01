source venv/bin/activate
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.gtf --paralogous /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.paralog_genes.80_percent_plus_identity.txt -o output -t 4 --raw /data2/ASTER/error_correction/raw_1D_only_data_no_mito_no_rDNA/gmap_bam/raw.bam --self /data2/ASTER/error_correction/canu/gmap_bam/canu.bam /data2/ASTER/error_correction/daccord/gmap_bam/daccord.bam /data2/ASTER/error_correction/daccord_trimmed/gmap_bam/daccord_trimmed.bam /data2/ASTER/error_correction/LoRMA/gmap_bam/LoRMA.bam /data2/ASTER/error_correction/MECAT/gmap_bam/MECAT.bam /data2/ASTER/error_correction/pbdagcon/gmap_bam/pbdagcon.bam /data2/ASTER/error_correction/pbdagcon_trimmed/gmap_bam/pbdagcon_trimmed.bam --hybrid /data2/ASTER/error_correction/HALC/gmap_bam/HALC.bam /data2/ASTER/error_correction/HALC_trimmed/gmap_bam/HALC_trimmed.bam /data2/ASTER/error_correction/HALC_split/gmap_bam/HALC_split.bam /data2/ASTER/error_correction/LoRDEC/gmap_bam/LoRDEC.bam /data2/ASTER/error_correction/LoRDEC_trimmed/gmap_bam/LoRDEC_trimmed.bam /data2/ASTER/error_correction/LoRDEC_split/gmap_bam/LoRDEC_split.bam /data2/ASTER/error_correction/LSC/gmap_bam/LSC.bam /data2/ASTER/error_correction/LSC_trimmed/gmap_bam/LSC_trimmed.bam /data2/ASTER/error_correction/NaS/gmap_bam/NaS.bam /data2/ASTER/error_correction/PBcR/gmap_bam/PBcR.bam /data2/ASTER/error_correction/proovread/gmap_bam/proovread.bam /data2/ASTER/error_correction/proovread_trimmed/gmap_bam/proovread_trimmed.bam --skip_bam_process --skip_alignqc --skip_copying

#separate runs (WAY better, I think...)

UNTRIMMED_UNSPLIT:
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.gtf --paralogous /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.paralog_genes.80_percent_plus_identity.txt -o output_untrimmed_unsplit -t 4 --raw /data2/ASTER/error_correction/raw_1D_only_data_no_mito_no_rDNA/gmap_bam/raw.bam --self /data2/ASTER/error_correction/canu/gmap_bam/canu.bam /data2/ASTER/error_correction/daccord/gmap_bam/daccord.bam /data2/ASTER/error_correction/MECAT/gmap_bam/MECAT.bam /data2/ASTER/error_correction/pbdagcon/gmap_bam/pbdagcon.bam --hybrid /data2/ASTER/error_correction/HALC/gmap_bam/HALC.bam /data2/ASTER/error_correction/LoRDEC/gmap_bam/LoRDEC.bam /data2/ASTER/error_correction/LSC/gmap_bam/LSC.bam /data2/ASTER/error_correction/NaS/gmap_bam/NaS.bam /data2/ASTER/error_correction/PBcR/gmap_bam/PBcR.bam /data2/ASTER/error_correction/proovread/gmap_bam/proovread.bam --skip_bam_process --skip_alignqc --skip_copying


TRIMMED:
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.gtf --paralogous /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.paralog_genes.80_percent_plus_identity.txt -o output_trimmed -t 4 --raw /data2/ASTER/error_correction/raw_1D_only_data_no_mito_no_rDNA/gmap_bam/raw.bam --self /data2/ASTER/error_correction/daccord_trimmed/gmap_bam/daccord_trimmed.bam /data2/ASTER/error_correction/pbdagcon_trimmed/gmap_bam/pbdagcon_trimmed.bam --hybrid /data2/ASTER/error_correction/HALC_trimmed/gmap_bam/HALC_trimmed.bam /data2/ASTER/error_correction/LoRDEC_trimmed/gmap_bam/LoRDEC_trimmed.bam /data2/ASTER/error_correction/LSC_trimmed/gmap_bam/LSC_trimmed.bam /data2/ASTER/error_correction/proovread_trimmed/gmap_bam/proovread_trimmed.bam --skip_bam_process --skip_alignqc --skip_copying


TRIMMED + SPLIT:
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.gtf --paralogous /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.paralog_genes.80_percent_plus_identity.txt -o output_trimmed_split -t 4 --raw /data2/ASTER/error_correction/raw_1D_only_data_no_mito_no_rDNA/gmap_bam/raw.bam --self /data2/ASTER/error_correction/LoRMA/gmap_bam/LoRMA.bam --hybrid /data2/ASTER/error_correction/HALC_split/gmap_bam/HALC_split.bam /data2/ASTER/error_correction/LoRDEC_split/gmap_bam/LoRDEC_split.bam --skip_bam_process --skip_alignqc --skip_copying


BEST_OF_ALL_CATEGORIES:
	SHALL WE DO THIS, SO THAT WE CAN COMPARE NORMAL W/ TRIMMED W/ TRIMMED+SPLIT???


SOME EXPLANATIONS:

LORMA IS IN TRIMMED + SPLIT BECAUSE:
	"The script lorma.sh runs LoRDEC iteratively, trims and splits the resulting initial corrected reads, and finally runs LoRMA on these trimmed and split reads.", from https://www.cs.helsinki.fi/u/lmsalmel/LoRMA/README.txt

PBCR = WHERE TO PUT???? ("The consensus module may trim/split/discard sequences if there is insufficient coverage to achieve a minimum QV for the corrected sequence.", from http://wgs-assembler.sourceforge.net/wiki/index.php/PBcR)
	FOR NOW, IN THE UNTRIMMED_UNSPLIT


deactivate
