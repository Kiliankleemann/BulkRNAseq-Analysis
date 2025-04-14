#!/bin/bash
echo "MACS2 Peak Calling"
mkdir MACS2_outputs
mkdir MACS2_outputs/broad
mkdir MACS2_outputs/narrow
cat sample_list.txt | while read sample; do
  macs2 callpeak -t BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam -f BAMPE --keep-dup all --mfold 5 50 -p 0.001 -g mm -n ${sample} --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/${sample}.log
done


#Making master peak file ATAC
echo 'Merging peak files and editing with bedtools'
#H3K9ac &H3K27ac
mkdir Master_peak_files
mkdir Master_peak_files/control
mkdir Master_peak_files/Sirt1_antagonist
mkdir Master_peak_files/Sirt1_agonists
multiIntersectBed -i  MACS2_outputs/narrow/control/*_peaks.narrowPeak  > Master_peak_files/control/Master_peak_merged.bed
multiIntersectBed -i  MACS2_outputs/narrow/Sirt1_antagonist/*_peaks.narrowPeak  > Master_peak_files/Sirt1_antagonist/Master_peak_merged.bed
multiIntersectBed -i  MACS2_outputs/narrow/Sirt1_agonists/*_peaks.narrowPeak  > Master_peak_files/Sirt1_agonists/Master_peak_merged.bed

bedtools merge -i Master_peak_files/control/Master_peak_merged.bed -d 200 -c 3 -o count > Master_peak_files/control/Master_peak_d200_3_merged.bed
bedtools merge -i Master_peak_files/Sirt1_antagonist/Master_peak_merged.bed -d 200 -c 3 -o count > Master_peak_files/Sirt1_antagonist/Master_peak_d200_3_merged.bed
bedtools merge -i Master_peak_files/Sirt1_agonists/Master_peak_merged.bed -d 200 -c 3 -o count > Master_peak_files/Sirt1_agonists/Master_peak_d200_3_merged.bed

#Single master peak file
multiIntersectBed -i  MACS2_outputs/narrow/control/*_peaks.narrowPeak MACS2_outputs/narrow/Sirt1_antagonist/*_peaks.narrowPeak MACS2_outputs/narrow/Sirt1_agonists/*_peaks.narrowPeak > Master_peak_files/Master_peak_merged.bed
bedtools merge -i Master_peak_files/Master_peak_merged.bed -d 200 -c 5 -o count > Master_peak_files/Master_peak_d200_5_merged.bed



### Call peaks with macs 3 and atac optimized
mkdir BED_files_sorted
cat sample_list.txt | while read sample; do
  macs3 filterdup --keep-dup all -f BAMPE -i BAM_files_final/${sample}.final.bam -o BED_files_sorted/${sample}.bedpe
done

mkdir macs3_output
cat sample_list.txt | while read sample; do
  macs3 hmmratac -i BED_files_sorted/${sample}.bedpe -f BEDPE -n macs3_output/${sample}
done

#Try macs2 with extra setting

mkdir MACS2_outputs/narrow_extra
cat sample_list.txt | while read sample; do
  macs2 callpeak -t BED_files_sorted/${sample}.bedpe -f BED --keep-dup all --nomodel --shift 100 --extsize 200 -g mm -n ${sample} --outdir MACS2_outputs/narrow_extra
done

#Single master peak file
multiIntersectBed -i  MACS2_outputs/narrow_extra/*_peaks.narrowPeak > Master_peak_files/Master_peak_merged_macs2_extra.bed
bedtools merge -i Master_peak_files/Master_peak_merged_macs2_extra.bed -d 300 -c 3 -o count > Master_peak_files/Master_peak_merged_macs2_extra_d300.bed

#Filter merged bedfile
awk '$4 > 3' Master_peak_files/Master_peak_merged_macs2_extra_d300.bed > Master_peak_files/Master_peak_merged_macs2_extra_d300_filtered.bed


#Extracting Data
#counts for control
mkdir Extracted_counts
mkdir Extracted_counts/control
bedtools multicov -bams BAM_files_sorted/control/*Aligned.sortedByCoord.out.bam -bed Master_peak_files/control/Master_peak_d200_3_merged.bed > Extracted_counts/control/Extracted_counts_d200_c3.txt
#counts for Sirt1_antagonist
mkdir Extracted_counts/Sirt1_antagonist
bedtools multicov -bams BAM_files_sorted/Sirt1_antagonist/*Aligned.sortedByCoord.out.bam -bed Master_peak_files/Sirt1_antagonist/Master_peak_d200_3_merged.bed > Extracted_counts/Sirt1_antagonist/Extracted_counts_d200_c3.txt
#counts for Sirt1_agonists
mkdir Extracted_counts/Sirt1_agonists
bedtools multicov -bams BAM_files_sorted/Sirt1_agonists/*Aligned.sortedByCoord.out.bam -bed Master_peak_files/Sirt1_agonists/Master_peak_d200_3_merged.bed > Extracted_counts/Sirt1_agonists/Extracted_counts_d200_c3.txt

#Extract counts from all bam files
mkdir Extracted_counts/all
bedtools multicov -bams BAM_files_final/*.bam -bed Master_peak_files/Master_peak_merged_macs2_extra_d300_filtered.bed > Extracted_counts/all/Extracted_counts_d200_c3.txt



#Peakannotation with Homer
#edit peak file for annotation
awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t+"}' Master_peak_files/Master_peak_merged_macs2_extra_d300_filtered.bed > Master_peak_files/Master_peak_merged_macs2_extra_d300_filtered_annotation.bed



#Individual groups
mkdir Annotated_peak
annotatePeaks.pl Master_peak_files/control/Master_peak_d200_3_merged.bed mm10 > Annotated_peak/Annotated_file_control_master_peak_merged_d200_c3.txt
annotatePeaks.pl Master_peak_files/Sirt1_antagonist/Master_peak_d200_3_merged.bed mm10 > Annotated_peak/Annotated_file_Sirt1_antagonist_master_peak_merged_d200_c3.txt
annotatePeaks.pl Master_peak_files/Sirt1_agonists/Master_peak_d200_3_merged.bed mm10 > Annotated_peak/Annotated_file_Sirt1_agonists_master_peak_merged_d200_c3.txt
#All
mkdir Annotated_peak/all
annotatePeaks.pl Master_peak_files/Master_peak_merged_macs2_extra_d300_filtered_annotation.bed mm10 > Annotated_peak/all/Master_peak_merged_macs2_extra_d300_filtered.txt


#BigwigMerge for vizualization UCSCtools
bigWigMerge BIGWIG_files/control/*.bw BIGWIG_files/control/Merged_bigwig_control.bedGraph
bigWigMerge BIGWIG_files/Sirt1_antagonist/*.bw BIGWIG_files/Sirt1_antagonist/Merged_bigwig_H3K9ac_Sirt1_antagonist.bedGraph
bigWigMerge BIGWIG_files/Sirt1_agonists/*.bw BIGWIG_files/Sirt1_agonists/Merged_bigwig_Sirt1_agonists.bedGraph

#control
for i in `ls BIGWIG_files/control/*bedGraph`
do
(head -n 1 $i && tail -n +2 $i | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $i.sort
bedGraphToBigWig $i.sort reference/mm10.chrom.sizes $i.bw 
done

#Sirt1_antagonist
for i in `ls BIGWIG_files/Sirt1_antagonist/*bedGraph`
do
(head -n 1 $i && tail -n +2 $i | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $i.sort
bedGraphToBigWig $i.sort reference/mm10.chrom.sizes $i.bw 
done

#Sirt1_agonists
for i in `ls BIGWIG_files/Sirt1_agonists/*bedGraph`
do
(head -n 1 $i && tail -n +2 $i | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $i.sort
bedGraphToBigWig $i.sort reference/mm10.chrom.sizes $i.bw 
done

#Plotting heatmap using deeptools
#H3K9
#500
computeMatrix reference-point --referencePoint center -S BIGWIG_files/control/Merged_bigwig_control.bedGraph.bw -R Master_peak_files/control/Master_peak_d200_3_merged.bed -p max/2 -a 500 -b 500 -o Matrix_d200_control_500.gz
computeMatrix reference-point --referencePoint center -S BIGWIG_files/Sirt1_antagonist/Merged_bigwig_Sirt1_antagonist.bedGraph.bw -R Master_peak_files/Sirt1_antagonist/Master_peak_d200_3_merged.bed -p max/2 -a 500 -b 500 -o Matrix_d200_Sirt1_antagonist_500.gz 
computeMatrix reference-point --referencePoint center -S BIGWIG_files/Sirt1_agonists/Merged_bigwig_Sirt1_agonists.bedGraph.bw -R Master_peak_files/Sirt1_agonists/Master_peak_d200_3_merged.bed -p max/2 -a 500 -b 500 -o Matrix_d200_Sirt1_agonists_500.gz

#1000
computeMatrix reference-point --referencePoint center -S BIGWIG_files/control/Merged_bigwig_control.bedGraph.bw -R Master_peak_files/control/Master_peak_d200_3_merged.bed -p max/2 -a 1000 -b 1000 -o Matrix_d200_control_1000.gz
computeMatrix reference-point --referencePoint center -S BIGWIG_files/Sirt1_antagonist/Merged_bigwig_Sirt1_antagonist.bedGraph.bw -R Master_peak_files/Sirt1_antagonist/Master_peak_d200_3_merged.bed -p max/2 -a 1000 -b 1000 -o Matrix_d200_Sirt1_antagonist_1000.gz
computeMatrix reference-point --referencePoint center -S BIGWIG_files/Sirt1_agonists/Merged_bigwig_Sirt1_agonists.bedGraph.bw -R Master_peak_files/Sirt1_agonists/Master_peak_d200_3_merged.bed -p max/2 -a 1000 -b 1000 -o Matrix_d200_Sirt1_agonists_1000.gz

#Plot Heatmap
plotHeatmap -m Matrix_d200_Sirt1_agonists_500.gz \
     -out Heatmap_Matrix_d200_Sirt1_agonists_500_white_red.png \
     --colorList white,red,darkred \
     --whatToShow 'heatmap and colorbar' \
     --averageTypeSummaryPlot mean \
     --heatmapHeight 30\
     --heatmapWidth 10\
     --missingDataColor white \
     --zMax 300 

plotHeatmap -m Matrix_H3K9ac_WT_1000.gz \
     -out Heatmap_Matrix_H3K9ac_WT_1000_white_red.png \
     --colorList white,red,darkred \
     --whatToShow 'heatmap and colorbar' \
     --averageTypeSummaryPlot mean \
     --heatmapHeight 30\
     --heatmapWidth 10\
     --missingDataColor white 

#Plot overall profile
plotProfile -m Matrix_d200_Sirt1_agonists_1000.gz \
              -out Matrix_d200_Sirt1_agonists_1000.png \
              --yMax 90

plotProfile -m Matrix_H3K9ac_1000.gz \
              -out Matrix_H3K9ac_1000.png \
              --yMax 70 


#ATAC
#ATAC
computeMatrix reference-point --referencePoint center -S BIGWIG_files/control/Merged_bigwig_control.bedGraph.bw -R Master_peak_files/control/Master_peak_d200_3_merged.bed -p max/2 -a 1000 -b 1000 -o Matrix_ATAC_d200_control_1000.gz

computeMatrix reference-point --referencePoint center -S Merged_bigwig_ATAC_WT.bedGraph.bw -R MACS2_outputs/master_peak/Master_ATAC_peak_merged_d300_c7_all.bed -p max/2 -a 1000 -b 1000 -o Matrix_ATAC_WT_1000.gz

#Plot Heatmap
plotHeatmap -m Matrix_ATAC_WT_1000.gz \
     -out Heatmap_Matrix_ATAC_WT_1000_white_red.png \
     --colorList white,blue,darkblue \
     --whatToShow 'heatmap and colorbar' \
     --averageTypeSummaryPlot mean \
     --heatmapHeight 30\
     --heatmapWidth 10\
     --missingDataColor white \
     --zMax 1000 

plotHeatmap -m Matrix_ATAC_KO_1000.gz \
     -out Heatmap_Matrix_ATAC_KO_1000_white_red.png \
     --colorList white,blue,darkblue \
     --whatToShow 'heatmap and colorbar' \
     --averageTypeSummaryPlot mean \
     --heatmapHeight 30\
     --heatmapWidth 10\
     --missingDataColor white \
     --zMax 1000 

#Plot overall profile
plotProfile -m Matrix_ATAC_WT_1000.gz \
              -out Matrix_ATAC_WT_1000.png \
              --yMax 330

plotProfile -m Matrix_ATAC_KO_1000.gz \
              -out Matrix_ATAC_KO_1000.png \
              --yMax 330


#Plotting PCA
#H3K9ac
multiBigwigSummary bins -b BIGWIG_files/control/H3K9ac3.bw BIGWIG_files/H3K9ac4.bw BIGWIG_files/H3K9ac5.bw -o H3K9_KO.npz
multiBigwigSummary bins -b BIGWIG_files/H3K9ac1.bw BIGWIG_files/H3K9ac2.bw BIGWIG_files/H3K9ac6.bw  BIGWIG_files/H3K9ac7.bw -o H3K9_WT.npz

#ATAC
multiBigwigSummary bins -b BIGWIG_files/ATA4.bw BIGWIG_files/ATA5.bw BIGWIG_files/ATA7.bw BIGWIG_files/ATA9.bw BIGWIG_files/ATA11.bw BIGWIG_files/ATA13.bw BIGWIG_files/ATA14.bw -o ATAC_KO.npz
multiBigwigSummary bins -b BIGWIG_files/ATA1.bw BIGWIG_files/ATA2.bw BIGWIG_files/ATA3.bw BIGWIG_files/ATA6.bw BIGWIG_files/ATA8.bw BIGWIG_files/ATA10.bw BIGWIG_files/ATA12.bw -o ATAC_WT.npz




