#mkdir reference
# cd reference
 #wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz 
 #wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
 #gzip -d *.gz
 #cd ..

 mkdir STAR_index_hg38
 STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_hg38 --genomeFastaFiles reference/GRCh38.p13.genome.fa  --sjdbGTFfile reference/GRCh38.refGene.gtf


 #echo 'Detected fastq samples'
 find ./lps_fastq -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
 cat sample_list.txt

mkdir BAM_files

cat sample_list.txt | while read sample; do
STAR --runThreadN 8 \
	--readFilesIn lps_fastq/${sample}.fastq.gz \
	--genomeDir /media/kilian/OS/STAR_index_hg38 \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix BAM_files/${sample} \
	--readFilesCommand zcat \
	--runMode alignReads \
	--outFilterMultimapNmax 1000 \
	--outMultimapperOrder Random \
	--winAnchorMultimapNmax 1000 \
	--outFilterMismatchNmax 3 \
	--alignEndsType EndToEnd \
	--alignIntronMax 1 \
	--alignMatesGapMax 350
done
