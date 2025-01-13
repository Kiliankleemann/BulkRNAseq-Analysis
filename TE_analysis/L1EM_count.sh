#Installing and running L1EM in mouse and human

###Mouse
# Have to run SAM unmapped within to retain unmapped reads.

#Generate mouse index
STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm39 --genomeFastaFiles GTF_files_TEtranscript/mm39.fa  #--sjdbGTFfile mm10_STAR/mm39.refGene.gtf

mkdir BAM_files_multi_mm39
cat sample_list.txt | while read sample; do
   STAR --runThreadN 8 \
	--readFilesIn fastq_files/trimmed_reads/${sample}.2.trimmed.fq.gz fastq_files/trimmed_reads/${sample}.1.trimmed.fq.gz \
	--genomeDir /media/kilian/OS/References/STAR_index_mm39 \
	--outSAMtype BAM SortedByCoordinate \
	--runMode alignReads \
	--outFilterMultimapNmax 100 \
	--outMultimapperOrder Random \
	--winAnchorMultimapNmax 100 \
	--outFilterMismatchNmax 3 \
	--alignEndsType EndToEnd \
	--alignIntronMax 1 \
	--alignMatesGapMax 350 \
	--outFileNamePrefix BAM_files_multi_mm39/${sample} \
	--readFilesCommand zcat \
	--outSAMunmapped Within

   samtools index BAM_files_multi_mm39/${sample}Aligned.sortedByCoord.out.bam
done

#Generate mouse genome index
bash generate_mm39_L1EM_fasta_and_index.sh /media/kilian/OS/References/GTF_files_TEtranscript/mm39.fa
bash /home/kilian/L1EM/run_L1EM_mm39.sh /media/kilian/DATA/Chronic_Acute/BMDM_IFNa4KI_PE_241217_P2024-266-RNA/BAM_files_multi_mm39/DE07NGSUKBR147099Aligned.sortedByCoord.out.bam /home/kilian/L1EM /media/kilian/OS/GTF_files_TEtranscript/BWA_mm39/mm39.fa

cat sample_list.txt | while read sample; do
	mkdir L1EM_${sample}
	cd L1EM_${sample}
	bash /home/kilian/L1EM/run_L1EM_mm39.sh \
/media/kilian/DATA/Chronic_Acute/BMDM_IFNa4KI_PE_241217_P2024-266-RNA/BAM_files_multi_mm39/${sample}Aligned.sortedByCoord.out.bam \
/home/kilian/L1EM /media/kilian/OS/GTF_files_TEtranscript/BWA_mm39/mm39.fa
	cd ..
done


####Human
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
zcat hg38.fa.gz > hg38.fa
bwa index hg38.fa

bash generate_L1EM_fasta_and_index.sh /media/kilian/OS/References/BWA_hg38/hg38.fa
