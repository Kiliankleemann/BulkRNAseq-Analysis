mkdir BAM_files_multi
cat sample_list.txt | while read sample; do
	STAR --runThreadN 12 \
	--readFilesIn fastq_files/trimmed_reads/${sample}_1.trimmed.fq.gz fastq_files/trimmed_reads/${sample}_2.trimmed.fq.gz\
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
