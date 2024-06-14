#!/bin/bash
#SBATCH --account=ag_ukikckp_behrendt
#SBATCH --job-name=Align_L1_STAR
#SBATCH --output=Align_L1_STAR.out
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5400M
#SBATCH --time=2-00:00:00

cat raw_fastq/sample_list.txt | while read sample; do
	STAR --runThreadN 4 \
	--readFilesIn raw_fastq/${sample}.fastq.gz \
	--genomeDir mm10_index/ \
	--outSAMtype BAM SortedByCoordinate > ${sample}.bam \
	--runMode alignReads \
	--outFilterMultimapNmax 1000 \
	--outMultimapperOrder Random \
	--winAnchorMultimapNmax 1000 \
	--outFilterMismatchNmax 3 \
	--alignEndsType EndToEnd \
	--alignIntronMax 1 \
	--alignMatesGapMax 350 \
	--outFileNamePrefix BAM_files/
done
