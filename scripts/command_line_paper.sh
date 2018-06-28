source venv/bin/activate
python /data2/ASTER/error_correction/LR_EC_analyser/run_LR_EC_analyser.py --genome /data2/ASTER/nanopore_mouse/genomic_data/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.gtf --raw /data2/ASTER/error_correction/raw_1D_only_data_no_mito_no_rDNA/gmap_bam/raw.bam -o output -t 4 --self /data2/ASTER/error_correction/canu/gmap_bam/canu.bam /data2/ASTER/error_correction/daccord/gmap_bam/daccord.bam /data2/ASTER/error_correction/daccord_trimmed/gmap_bam/daccord_trimmed.bam /data2/ASTER/error_correction/LoRMA/gmap_bam/LoRMA.bam /data2/ASTER/error_correction/MECAT/gmap_bam/MECAT.bam /data2/ASTER/error_correction/pbdagcon/gmap_bam/pbdagcon.bam /data2/ASTER/error_correction/pbdagcon_trimmed/gmap_bam/pbdagcon_trimmed.bam --hybrid /data2/ASTER/error_correction/LoRDEC/gmap_bam/LoRDEC.bam /data2/ASTER/error_correction/NaS/gmap_bam/NaS.bam /data2/ASTER/error_correction/PBcR/gmap_bam/PBcR.bam /data2/ASTER/error_correction/proovread_untrimmed/gmap_bam/proovread_untrimmed.bam /data2/ASTER/error_correction/proovread_trimmed/gmap_bam/proovread_trimmed.bam --skip_bam_process --skip_alignqc --skip_copying --paralogous /data2/ASTER/nanopore_mouse/transcriptomic_data/Mus_musculus.GRCm38.87.paralog_genes.80_percent_plus_identity.txt
deactivate
