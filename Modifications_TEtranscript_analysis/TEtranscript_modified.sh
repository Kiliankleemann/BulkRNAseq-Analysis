#Running TEtranscript
#Download required files (for mouse)
mkdir GTF_files_TEtranscript
cd GTF_files_TEtranscript
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm38_GENCODE_rmsk_TE.gtf.gz.download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz

#filtering out only instances of full transcripts
/media/kilian/OS/BulkRNAseq-Analysis/Modifications_TEtranscript_analysis/filter_gtf.sh mm10.refGene.gtf > mm10.transcripts_filtered.gtf

#Removing all TE counts from repeatmasker which overlap with refGene.gtf
bedtools subtract -s -A -a GRCm38_GENCODE_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > GRCm38_GENCODE_rmsk_TE_filtered.gtf

#Finding overlaping 
bedtools intersect -s -a GRCm38_GENCODE_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > GRCm38_rmsk_TE_transcript_overlap.gtf


#Multimapping
#TEtranscript
mkdir TETranscripts_multi
TEtranscripts -t  \
                     -c 
                     --GTF /media/kilian/OS/GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE /media/kilian/OS/GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE_filtered.gtf \
                     --mode multi \
                     --outdir TETranscripts_multi \
                     --sortByPos



#Quantification of TEs with Stringtie
mkdir stringtie_output_TE_RMSK
cat sample_list.txt | while read sample; 
	do stringtie -e -B -p 8 -G reference/GRCm38_GENCODE_rmsk_TE.gtf -o stringtie_output_TE_RMSK/${sample}/${sample}.gtf BAM_files/${sample}Aligned.sortedByCoord.out.bam 
done

mkdir stringtie_output_TE_RMSK
cat sample_list.txt | while read sample; 
	do stringtie -e -B -p 8 -G reference/GRCm38_GENCODE_rmsk_TE.gtf -o stringtie_output_TE_RMSK/${sample}/${sample}.gtf BAM_files/${sample}Aligned.sortedByCoord.out.bam 
done

#Running Python script to extract raw counts
cd stringtie_output_TE_RMSK
python /media/kilian/OS/Melani_THP1_August_2023/scripts_analysis/prepDE.py3



