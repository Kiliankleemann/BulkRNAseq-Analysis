

### STAR ALIGNMENT
#STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm10 --genomeFastaFiles mm10_STAR/mm10.fa  --sjdbGTFfile mm10_STAR/mm10.refGene.gtf
#CHANGE genomeDIR to STAR_index_mm10 for mouse!!!!
mkdir BAM_files_multi
cat sample_list.txt | while read sample; do
	STAR --runThreadN 12 \
	--readFilesIn fastq_files/trimmed_reads/${sample}.trimmed.fq.gz \
	--genomeDir /media/kilian/OS/References/STAR_index_hg38 \
	--outSAMtype BAM SortedByCoordinate  \
	--runMode alignReads \
	--outFilterMultimapNmax 100 \
	--outMultimapperOrder Random \
	--winAnchorMultimapNmax 100 \
	--outFilterMismatchNmax 3 \
	--alignEndsType EndToEnd \
	--alignIntronMax 1 \
	--alignMatesGapMax 350 \
	--outFileNamePrefix BAM_files_multi/${sample} \
	--readFilesCommand zcat 
done


cat sample_list.txt | while read sample; do
	samtools index BAM_multi/${sample}Aligned.sortedByCoord.out.bam
done

#Running TEtranscript
#Download required files
mkdir GTF_files_TEtranscript
cd GTF_files_TEtranscript
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm38_GENCODE_rmsk_TE.gtf.gz.download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz

#Multimapping
mkdir TETranscripts_multi
TEtranscripts -t /media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE06NGSLABR103437_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE17NGSLABR103433_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE22NGSLABR103440_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE28NGSLABR103429_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE33NGSLABR103436_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE38NGSLABR103443_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE44NGSLABR103432_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE49NGSLABR103439_1Aligned.sortedByCoord.out.bam \
                     -c /media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE55NGSLABR103428_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE60NGSLABR103435_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE65NGSLABR103442_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE71NGSLABR103431_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE76NGSLABR103438_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE87NGSLABR103434_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE92NGSLABR103441_1Aligned.sortedByCoord.out.bam \
/media/kilian/DATA/AvM_Aprl_2024_240405_P2024-039-LEX-LB/BAM_files_multi/DE98NGSLABR103430_1Aligned.sortedByCoord.out.bam \
                     --GTF /media/kilian/OS/GTF_files_TEtranscript/GRCh38.refGene.gtf \
                     --TE /media/kilian/OS/GTF_files_TEtranscript/GRCh38_GENCODE_rmsk_TE.gtf \
                     --mode multi \
                     --outdir TETranscripts_multi \
                     --sortByPos
                     
                    


#Quantification with TElocal
mkdir TElocal_multi
cat sample_list.txt | while read sample; do
 TElocal -b BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam  \
		--GTF /media/kilian/OS/GTF_files_TEtranscript/GRCh38.refGene.gtf \
		--TE /media/kilian/OS/GTF_files_TEtranscript/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
		--mode multi \
		--project TElocal_multi/${sample} \
		--sortByPos
done



#Quantification of TEs with Stringtie
mkdir stringtie_output_TE_RMSK
cat sample_list.txt | while read sample; 
	do stringtie -e -B -p 8 -G reference/GRCm38_GENCODE_rmsk_TE.gtf -o stringtie_output_TE_RMSK/${sample}/${sample}.gtf BAM_files/${sample}Aligned.sortedByCoord.out.bam 
done

#Running Python script to extract raw counts
cd stringtie_output_TE_RMSK
python /media/kilian/OS/Melani_THP1_August_2023/scripts_analysis/prepDE.py3


#Trying with full length bed file
scripts/bed2gtf.py -s reference/L1base_FL_mm10.bed  

