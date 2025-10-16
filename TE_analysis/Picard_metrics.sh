

#Index files
cat sample_list.txt | while read sample; do
	samtools index BAM_files_multi_trimmed/${sample}Aligned.sortedByCoord.out.bam
done

#Output metrics
mkdir Picard_metrics
cat sample_list.txt | while read sample; do
picard CollectRnaSeqMetrics \
    I=BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam \
    O=Picard_metrics/${sample}_output_rna_metrics.txt \
    REF_FLAT='/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Database_GTF_files_TE/hg38.refFlat.txt' \
    STRAND_SPECIFICITY=NONE \
    R='/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Database_GTF_files_TE/Reference_LINE1/human/GRCh38.p13.genome.fa'
done

#Picard metric strand specific
mkdir Picard_metrics_first_strand_spec
cat sample_list.txt | while read sample; do
picard CollectRnaSeqMetrics \
    I=BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam \
    O=Picard_metrics_first_strand_spec/${sample}_output_rna_metrics.txt \
    REF_FLAT='/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Database_GTF_files_TE/mm10.refFlat.txt' \
    STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
    R='/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Database_GTF_files_TE/Reference_LINE1/mouse/mm10.fa'
done


