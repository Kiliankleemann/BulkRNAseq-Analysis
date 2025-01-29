#Running Bedtools to check overlap of reads with L1base bed annotations

# Step 1: Filter BAM for properly paired, non-chimeric reads
mkdir BAM_files_filtered
cat sample_list.txt | while read sample; do
    samtools view -h -f 2 -F 2048 BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam > BAM_files_filtered/${sample}filtered.bam

    samtools sort -o BAM_files_filtered/${sample}filtered.bam

    samtools index BAM_files_filtered/${sample}filtered.bam

done

#Fractional change
samtools view input.bam | awk '{for (i=12; i<=NF; i++) if ($i ~ /^NH:i:/) {split($i, nh, ":"); printf("%s\t%s\t%s\t%.2f\n", $1, $2, $3, 1/nh[3]); next}}' > fractional.bam

# Step 3: Run strand-specific coverage
bedtools coverage -a features.bed -b filtered_sorted.bam -s -counts > strand_specific_counts.txt

# Step 4: Handle reverse strand (if needed)
bedtools coverage -a features.bed -b filtered_sorted.bam -S -counts > reverse_strand_counts.txt

bedtools coverage -a features.bed -b filtered.bam -s -f 1.0 > counts.txt


mkdir bedtools_out_mmorf2l1_8438
cat sample_list.txt | while read sample; do
    bedtools coverage -a /media/kilian/OS/References/GTF_files_featureCounts_mm10/mmorf2l1_8438.bed\
     -b BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam -s -counts -f 0.8 > bedtools_out_mmorf2l1_8438/${sample}_bedtool_counts_output.txt
done

