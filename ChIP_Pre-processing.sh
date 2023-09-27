#!/bin/bash
echo 'Detected following samples'
find ./fastq -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d'_' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#echo 'Unzipping files'
#gunzip fastq/*.fastq.gz

echo 'Checking quality of untrimmed reads'
mkdir QC
mkdir QC/untrimmed
fastqc fastq/*.fastq.gz -t 2 -o 'QC/untrimmed'
multiqc QC/untrimmed -o QC/untrimmed

echo 'Trimming reads - stand by'
mkdir fastq/trimmed_reads
cat sample_list.txt | while read sample; do
	echo $sample
	cutadapt --cores 0 --minimum-length 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o fastq/trimmed_reads/${sample}.1.trimmed.fastq.gz -p fastq/trimmed_reads/${sample}.2.trimmed.fastq.gz fastq/${sample}_R1_001.fastq fastq/${sample}_R2_001.fastq
done

echo 'Checking quality of trimmed_reads'
mkdir QC/trimmed
fastqc fastq/trimmed_reads/*.fastq.gz -t 4 -o 'QC/trimmed'
multiqc QC/trimmed -o QC/trimmed 

echo 'Performing Bowtie2 alignment to SAM files'
mkdir SAM_files
cat sample_list.txt | while read sample; do
	bowtie2 -p 4 -q --local -x ~/Desktop/mm10/mm10 -1 fastq/trimmed_reads/${sample}.1.trimmed.fastq.gz -2 fastq/trimmed_reads/${sample}.2.trimmed.fastq.gz -S SAM_files/${sample}.unsorted.sam
done

#echo 'Zipping untrimmed files'
#gzip fastq/*.fastq

echo "Creating BAM files from SAM files"
mkdir BAM_files
mkdir BAM_files/unsorted
cat sample_list.txt | while read sample; do
	samtools view -h -S -b -o BAM_files/unsorted/${sample}.unsorted.bam SAM_files/${sample}.unsorted.sam
done

echo "Sorting BAM files"
mkdir BAM_files/sorted
cat sample_list.txt | while read sample; do
	sambamba sort -t 4 -o BAM_files/sorted/${sample}.sorted.bam BAM_files/unsorted/${sample}.unsorted.bam
done

echo "Filtering sorted BAM files"
mkdir BAM_files/final
cat sample_list.txt | while read sample; do
	sambamba view -h -t 4 -f bam -F "[XS] == null and not unmapped and not duplicate"  BAM_files/sorted/${sample}.sorted.bam > BAM_files/final/${sample}.final.bam
done

echo "Indexing final BAM files"
cat sample_list.txt | while read sample; do
	samtools index BAM_files/final/${sample}.final.bam
done

echo "Alignment statistics"
cat sample_list.txt | while read sample; do
	samtools flagstat BAM_files/sorted/${sample}.sorted.bam
done

echo"remove mitochondrial reads"
cat sample_list.txt | while read sample; do
	samtools  view -h BAM_files/sorted/${sample}.sorted.bam  |  removeChrom - - chrM  |  samtools view -b - > BAM_files/sorted/${sample}.sorted.bam
done


echo "Making BIGWIG_files"
mkdir BIGWIG_files
cat sample_list.txt | while read sample; do
	bamCoverage -b BAM_files/final/${sample}.final.bam -o BIGWIG_files/${sample}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads --numberOfProcessors 4 --skipNonCoveredRegions
done





















