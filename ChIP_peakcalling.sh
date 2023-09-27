#!/bin/bash

echo "MACS2 Peak Calling"
mkdir MACS2_outputs
mkdir MACS2_outputs/broad
mkdir MACS2_outputs/narrow
macs2 callpeak -t BAM_files/final/H3K9ac1.final.bam -c BAM_files/final/Input1.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac1 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac1_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac2.final.bam -c BAM_files/final/Input2.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac2 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac2_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac3.final.bam -c BAM_files/final/Input3.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac3 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac3_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac4.final.bam -c BAM_files/final/Input4.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac4 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac4_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac5.final.bam -c BAM_files/final/Input5.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac5 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac5_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac6.final.bam -c BAM_files/final/Input6.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac6 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac6_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac7.final.bam -c BAM_files/final/Input7.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac7 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac7_macs2.log



#ATAC
macs2 callpeak -t BAM_files/final/ATA1.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA1 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA1_macs2.log
macs2 callpeak -t BAM_files/final/ATA2.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA2 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA2_macs2.log
macs2 callpeak -t BAM_files/final/ATA3.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA3 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA3_macs2.log
macs2 callpeak -t BAM_files/final/ATA4.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA4 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA4_macs2.log
macs2 callpeak -t BAM_files/final/ATA5.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA5 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA5_macs2.log
macs2 callpeak -t BAM_files/final/ATA6.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA6 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA6_macs2.log
macs2 callpeak -t BAM_files/final/ATA7.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA7 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA7_macs2.log
macs2 callpeak -t BAM_files/final/ATA8.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA8 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA8_macs2.log
macs2 callpeak -t BAM_files/final/ATA9.final.bam  -f BAMPE --mfold 5 50 -p 0.001 -g mm -n ATA9 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/ATA9_macs2.log



#H3K27ac
macs2 callpeak -t BAM_files/final/H3K27ac1.final.bam -c BAM_files/final/Input1.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac1 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac1_macs2.log
macs2 callpeak -t BAM_files/final/H3K27ac2.final.bam -c BAM_files/final/Input2.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac2 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac2_macs2.log
macs2 callpeak -t BAM_files/final/H3K27ac3.final.bam -c BAM_files/final/Input3.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac3 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac3_macs2.log
macs2 callpeak -t BAM_files/final/H3K27ac4.final.bam -c BAM_files/final/Input4.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac4 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac4_macs2.log
macs2 callpeak -t BAM_files/final/H3K27ac5.final.bam -c BAM_files/final/Input5.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac5 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac5_macs2.log
macs2 callpeak -t BAM_files/final/H3K27ac6.final.bam -c BAM_files/final/Input6.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac6 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac6_macs2.log
macs2 callpeak -t BAM_files/final/H3K27ac7.final.bam -c BAM_files/final/Input7.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K27ac7 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K27ac7_macs2.log



#Making master peak file ATAC
echo 'Merging peak files and editing with bedtools'
multiIntersectBed -i  MACS2_outputs/narrow/ATA1_peaks.narrowPeak MACS2_outputs/narrow/ATA2_peaks.narrowPeak MACS2_outputs/narrow/ATA3_peaks.narrowPeak MACS2_outputs/narrow/ATA4_peaks.narrowPeak MACS2_outputs/narrow/ATA5_peaks.narrowPeak MACS2_outputs/narrow/ATA6_peaks.narrowPeak MACS2_outputs/narrow/ATA7_peaks.narrowPeak MACS2_outputs/narrow/ATA8_peaks.narrowPeak MACS2_outputs/narrow/ATA9_peaks.narrowPeak > Master_ATAC_peak_merged.bed
bedtools merge -i Master_ATAC_peak_merged.bed -d 300 -c 7 -o count > Master_ATAC_peak_merged_d300_c7_all.bed

#Extracting Data
bedtools multicov -bams BAM_files/final/ATA1.final.bam BAM_files/final/ATA2.final.bam BAM_files/final/ATA3.final.bam BAM_files/final/ATA4.final.bam BAM_files/final/ATA5.final.bam BAM_files/final/ATA6.final.bam BAM_files/final/ATA7.final.bam BAM_files/final/ATA8.final.bam BAM_files/final/ATA9.final.bam BAM_files/final/ATA10.final.bam BAM_files/final/ATA11.final.bam BAM_files/final/ATA12.final.bam BAM_files/final/ATA13.final.bam BAM_files/final/ATA14.final.bam -bed Master_peak/Master_peak_merged_d300_c7_all.bed > Extracted_counts_d300_c7_all


#Making master peak file H3K9ac
echo 'Merging peak files and editing with bedtools'
multiIntersectBed -i  MACS2_outputs/narrow/H3K9ac1_peaks.narrowPeak MACS2_outputs/narrow/H3K9ac2_peaks.narrowPeak MACS2_outputs/narrow/H3K9ac3_peaks.narrowPeak MACS2_outputs/narrow/H3K9ac4_peaks.narrowPeak MACS2_outputs/narrow/H3K9ac5_peaks.narrowPeak MACS2_outputs/narrow/H3K9ac6_peaks.narrowPeak MACS2_outputs/narrow/H3K9ac7_peaks.narrowPeak > Master_peak_merged_H3K9ac.bed  
bedtools merge -i Master_peak_merged_H3K9ac.bed -d 300 -c 3 -o count > Master_peak_merged_H3K9ac_d300_c3_all.bed  

#Extracting Data
bedtools multicov -bams BAM_files/final/H3K9ac1.final.bam BAM_files/final/H3K9ac2.final.bam BAM_files/final/H3K9ac3.final.bam BAM_files/final/H3K9ac4.final.bam BAM_files/final/H3K9ac5.final.bam BAM_files/final/H3K9ac6.final.bam BAM_files/final/H3K9ac7.final.bam -bed Master_peak_merged_H3K9ac_d300_c3_all.bed  > Extractec_counts_H3K9ac_d300_c3_all 

#Peakannotation with Homer
annotatePeaks.pl Annotation_file_ATAC_master_peak_merged_d300_c7.txt mm10 > Annotated_file_ATAC_master_peak_merged_d300_c7
annotatePeaks.pl Annotation_file_H3K9ac_master_peak_merged_d300_c3.txt mm10 > Annotated_file_H3K9ac_master_peak_merged_d300_c3


annotatePeaks.pl Annotation_file_ATAC_padj05.txt mm10 > Annotated_file_ATAC_padj05
annotatePeaks.pl Annotation_file_H3K9ac_padj05.txt mm10 > Annotated_file_H3K9ac_padj05


#Motif discovery ATAC using HOMER
findMotifsGenome.pl Annotation_file_DPs_UPREGULATED_cre-vs-wt_ATAC_padj05.txt mm10 Motif_discovery/UPREGULATED_cre-vs-wt_ATAC_padj05/ -size given -mask
findMotifsGenome.pl Annotation_file_DPs_DOWNREGULATED_cre-vs-wt_ATAC_padj05.txt mm10 Motif_discovery/DOWNREGULATED_cre-vs-wt_ATAC_padj05/ -size given -mask


#Motif discovery H3K9ac using HOMER
findMotifsGenome.pl Annotation_file_DPs_UPREGULATED_cre-vs-wt_H3K9ac_padj05.txt mm10 Motif_discovery/UPREGULATED_cre-vs-wt_H3K9ac_padj05/ -size given -mask
findMotifsGenome.pl Annotation_file_DPs_DOWNREGULATED_cre-vs-wt_H3K9ac_padj05.txt mm10 Motif_discovery/DOWNREGULATED_cre-vs-wt_H3K9ac_padj05/ -size given -mask

#BigwigMerge for vizualization UCSCtools
#H3K9ac
bigWigMerge BIGWIG_files/H3K9ac3.bw BIGWIG_files/H3K9ac4.bw BIGWIG_files/H3K9ac5.bw  Merged_bigwig_H3K9ac.bedGraph
bigWigMerge BIGWIG_files/H3K9ac1.bw BIGWIG_files/H3K9ac2.bw BIGWIG_files/H3K9ac6.bw  BIGWIG_files/H3K9ac7.bw Merged_bigwig_H3K9ac_WT.bedGraph

#ATAC
bigWigMerge BIGWIG_files/ATA1.bw BIGWIG_files/ATA2.bw BIGWIG_files/ATA3.bw BIGWIG_files/ATA6.bw BIGWIG_files/ATA8.bw BIGWIG_files/ATA10.bw BIGWIG_files/ATA12.bw Merged_bigwig_ATAC_WT.bedGraph
bigWigMerge BIGWIG_files/ATA4.bw BIGWIG_files/ATA5.bw BIGWIG_files/ATA7.bw BIGWIG_files/ATA9.bw BIGWIG_files/ATA11.bw BIGWIG_files/ATA13.bw BIGWIG_files/ATA14.bw Merged_bigwig_ATAC_KO.bedGraph


for i in `ls *bedGraph`
do
(head -n 1 $i && tail -n +2 $i | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $i.sort
bedGraphToBigWig $i.sort mm10.chrom.size.txt $i.bw 
done


#Plotting heatmap using deeptools
#H3K9
#500
computeMatrix reference-point --referencePoint center -S Merged_bigwig_H3K9ac.bedGraph.bw -R MACS2_outputs/master_peak/Master_peak_merged_H3K9ac_d300_c3_all.bed -p max/2 -a 500 -b 500 -o Matrix_H3K9ac_500.gz
computeMatrix reference-point --referencePoint center -S Merged_bigwig_H3K9ac_WT.bedGraph.bw -R MACS2_outputs/master_peak/Master_peak_merged_H3K9ac_d300_c3_all.bed -p max/2 -a 500 -b 500 -o Matrix_H3K9ac_WT_500.gz

#1000
computeMatrix reference-point --referencePoint center -S BIGWIG_files/Merged_bigwig_H3K9ac.bedGraph.bw -R MACS2_outputs/master_peak/Master_peak_merged_H3K9ac_d300_c3_all.bed -p max/2 -a 1000 -b 1000 -o Matrix_H3K9ac_1000.gz
computeMatrix reference-point --referencePoint center -S BIGWIG_files/Merged_bigwig_H3K9ac_WT.bedGraph.bw -R MACS2_outputs/master_peak/Master_peak_merged_H3K9ac_d300_c3_all.bed -p max/2 -a 1000 -b 1000 -o Matrix_H3K9ac_WT_1000.gz

#Plot Heatmap
plotHeatmap -m Matrix_H3K9ac_1000.gz \
     -out Heatmap_Matrix_H3K9ac_1000_white_red.png \
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
plotProfile -m Matrix_H3K9ac_WT_1000.gz \
              -out Matrix_H3K9ac_WT_1000.png \
              --yMax 70

plotProfile -m Matrix_H3K9ac_1000.gz \
              -out Matrix_H3K9ac_1000.png \
              --yMax 70 


#ATAC
#ATAC
computeMatrix reference-point --referencePoint center -S Merged_bigwig_ATAC_KO.bedGraph.bw -R MACS2_outputs/master_peak/Master_ATAC_peak_merged_d300_c7_all.bed -p max/2 -a 1000 -b 1000 -o Matrix_ATAC_KO_1000.gz
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
multiBigwigSummary bins -b BIGWIG_files/H3K9ac3.bw BIGWIG_files/H3K9ac4.bw BIGWIG_files/H3K9ac5.bw -o H3K9_KO.npz
multiBigwigSummary bins -b BIGWIG_files/H3K9ac1.bw BIGWIG_files/H3K9ac2.bw BIGWIG_files/H3K9ac6.bw  BIGWIG_files/H3K9ac7.bw -o H3K9_WT.npz

#ATAC
multiBigwigSummary bins -b BIGWIG_files/ATA4.bw BIGWIG_files/ATA5.bw BIGWIG_files/ATA7.bw BIGWIG_files/ATA9.bw BIGWIG_files/ATA11.bw BIGWIG_files/ATA13.bw BIGWIG_files/ATA14.bw -o ATAC_KO.npz
multiBigwigSummary bins -b BIGWIG_files/ATA1.bw BIGWIG_files/ATA2.bw BIGWIG_files/ATA3.bw BIGWIG_files/ATA6.bw BIGWIG_files/ATA8.bw BIGWIG_files/ATA10.bw BIGWIG_files/ATA12.bw -o ATAC_WT.npz




