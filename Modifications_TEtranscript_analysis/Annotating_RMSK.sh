#Annotating repeatmasker file in human and mouse
convert2bed -i gtf < GRCh38_GENCODE_rmsk_TE.gtf > GRCh38_GENCODE_rmsk_TE.bed
annotatePeaks.pl GRCh38_GENCODE_rmsk_TE.bed hg38 > Annotated_GRCh38_GENCODE_rmsk_TE.gtf


convert2bed -i gtf < GRCm38_GENCODE_rmsk_TE.gtf > GRCm38_GENCODE_rmsk_TE.bed
annotatePeaks.pl GRCm38_GENCODE_rmsk_TE.gtf mm10 > Annotated_GRCm38_GENCODE_rmsk_TE.gtf
