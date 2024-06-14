#!/bin/bash
#SBATCH --account=ag_ukikckp_behrendt
#SBATCH --job-name=Align_DNA_damage_STAR_mm10
#SBATCH --output=Align_DNA_damage_STAR_mm10.out
#SBATCH --ntasks=2
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5400M
#SBATCH --time=2-00:00:00


module load STAR/2.7.10b-GCC-11.3.0

# mkdir reference
# cd reference
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz 
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
# gzip -d *.gz
# cd ..

# mkdir STAR_index_hg38
# STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_hg38 --genomeFastaFiles reference/GRCh38.p13.genome.fa  --sjdbGTFfile reference/GRCh38.refGene.gtf


# #echo 'Detected fastq samples'
# find ./raw_fastq -name "*.fq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
# cat sample_list.txt

mkdir BAM_files

cat sample_list.txt | while read sample; do
STAR --runThreadN 16 \
	--readFilesIn raw_fastq/${sample}.fq.gz raw_fastq/${sample}.fq.gz \
	--genomeDir ~/STAR_index_hg38 \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix BAM_files/${sample}.sorted.bam \
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


 
	


