#Create custom mm10 reference star index


#merge gene and rmsk mm10 gtf files
cat mm10.refGene.gtf mm10_rmsk.gtf > mm10_refGene_rmsk_combined.gtf


STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm10_Gene_TE_merged --genomeFastaFiles STAR_reference_files/mm10.fa  --sjdbGTFfile STAR_reference_files/mm10_refGene_rmsk_combined.gtf --sjdbOverhang 99 \
     --limitSjdbInsertNsj 5000000


mkdir BAM_files_multi_index_merged_custom
cat sample_list.txt | while read sample; do
   STAR --runThreadN 8 \
	--readFilesIn fastq_files/trimmed_reads/${sample}.1.trimmed.fq.gz fastq_files/trimmed_reads/${sample}.2.trimmed.fq.gz \
	--genomeDir /media/kilian/OS/References/STAR_index_mm10_Gene_TE_merged \
	--outSAMtype BAM SortedByCoordinate \
	--runMode alignReads \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100 \
	--outMultimapperOrder Random \
	--runRNGseed 777 \
	--outFilterType BySJout \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outFileNamePrefix BAM_files_multi_index_merged_custom/${sample} \
	--readFilesCommand zcat
done



