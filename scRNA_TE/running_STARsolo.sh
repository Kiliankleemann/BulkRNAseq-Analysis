#Running STARsolo

echo 'Detected fastq samples'
find fastq -name "*_1.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '_' -f1 | sort -u > sample_list.txt
cat sample_list.txt


cat sample_list.txt | while read sample; do
	STAR --runMode alignReads \
	--soloType Droplet \
	--genomeDir STAR_index_hg38 \
	--readFilesIn fastq/${sample}_2.fastq.gz fastq/${sample}_1.fastq.gz \
	--soloCBwhitelist reference/3M-february-2018.txt \
	--outSAMtype BAM Unsorted \
	--outSAMattributes NH HI AS nM CR CY UR UY \
	--readFilesCommand zcat \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100 \
	--outMultimapperOrder Random \
	--runRNGseed 777 \
	--outSAMmultNmax 1 \
	--twopassMode Basic \
	--soloUMIlen 10 \
	--soloMultiMappers EM \
	--runThreadN 16 \
	--outFileNamePrefix ${sample}/ \
	--outReadsUnmapped Fastx \
	--outSAMunmapped Within \
	--soloBarcodeReadLength 0
done