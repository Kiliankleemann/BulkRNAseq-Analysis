#Running TEtranscript
#Download required files (for mouse)
mkdir GTF_files_TEtranscript
cd GTF_files_TEtranscript
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm38_GENCODE_rmsk_TE.gtf.gz.download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz

#filtering out only instances of full transcripts MOUSE
/media/kilian/OS/BulkRNAseq-Analysis/Modifications_TEtranscript_analysis/filter_gtf.sh gencode.vM10.annotation.gtf > gencode.vM10.transcripts_filtered.gtf
#filtering out only instances of full transcripts HUMAN
/media/kilian/OS/BulkRNAseq-Analysis/Modifications_TEtranscript_analysis/filter_gtf.sh GRCh38.refGene.gtf > GRCh38.transcripts_filtered.gtf

#Removing all TE counts from repeatmasker which overlap with refGene.gtf
bedtools subtract -A -a GRCm38_GENCODE_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > GRCm38_GENCODE_rmsk_TE_filtered.gtf

bedtools subtract -A -a GRCh38_GENCODE_rmsk_TE.gtf -b GRCh38.transcripts_filtered.gtf > GRCm38_GENCODE_rmsk_TE_filtered.gtf


#Only L1Md family members
#Filter out only L1Md_T elements
awk '$0 ~ /gene_id "L1Md/' GRCm38_GENCODE_rmsk_TE.gtf > GRCm38_GENCODE_L1Md.gtf


#Filter out Exon overlapping


#Filter out L1Md overlapping annotated genes
bedtools subtract -A -a GRCm38_GENCODE_L1Md.gtf -b mm10.transcripts_filtered.gtf > GRCm38_GENCODE_rmsk_TE_filtered.gtf



convert2bed -i gtf < GRCm38_GENCODE_rmsk_TE.gtf > GRCm38_GENCODE_rmsk_TE.bed
convert2bed -i gtf < gencode.vM10.transcripts_filtered.gtf > gencode.vM10.annotation.bed

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}' gencode.vM10.annotation.bed > simplified_annotation.bed
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}' GRCm38_GENCODE_rmsk_TE.bed > simplified_TE.bed

sort -k1,1 -k2,2n simplified_TE.bed > sorted_TE.bed
sort -k1,1 -k2,2n simplified_annotation.bed > sorted_annotation.bed


bedtools subtract -A -a mm10_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > mm10_TE_intronfiltered.gtf

bedtools intersect -a mm10_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > mm10_TE_intron_intersect.gtf

#Finding overlaping regions between repeatmakser and gene annotations MOUSE 
bedtools intersect -s -a GRCm38_GENCODE_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > GRCm38_rmsk_TE_transcript_overlap.gtf
bedtools intersect -s -wa -a GRCm38_GENCODE_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > GRCm38_rmsk_TE_transcript_overlap_wa.gtf
bedtools intersect -wa -a GRCm38_GENCODE_rmsk_TE.gtf -b mm10.transcripts_filtered.gtf > GRCm38_rmsk_TE_transcript_overlap_wa_nostrandedness.gtf

#Finding overlaping regions between repeatmakser and gene annotations HUMAN 
bedtools intersect -s -a GRCh38_GENCODE_rmsk_TE.gtf -b GRCh38.transcripts_filtered.gtf > GRCh38_rmsk_TE_transcript_overlap.gtf
bedtools intersect -s -wa -a GRCh38_GENCODE_rmsk_TE.gtf -b GRCh38.transcripts_filtered.gtf > GRCh38_rmsk_TE_transcript_overlap_wa.gtf
bedtools intersect -wa -a GRCh38_GENCODE_rmsk_TE.gtf -b GRCh38.transcripts_filtered.gtf > GRCh38_rmsk_TE_transcript_overlap_wa_nostrandedness.gtf


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



