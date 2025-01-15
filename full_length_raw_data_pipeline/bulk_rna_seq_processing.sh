mkdir fastq_files
mv *.fastq.gz fastq_files

#Edit fastq  ending of file and symbol to cut the name 
find ./fastq_files -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '_' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#Running fastq for all files
mkdir QC
fastqc fastq_files/*.fastq.gz -t 12 -o QC

#Move all fastqc reports into the QC folder 
mv fastq_files/*fastqc* QC/

#Make quality overview
multiqc QC/



#Running cutadapt for adapter trimming 
mkdir fastq_files/trimmed_reads
cat sample_list.txt | while read sample; do
	cutadapt --cores 12 --minimum-length 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o fastq_files/trimmed_reads/${sample}.1.trimmed.fastq.gz -p fastq_files/trimmed_reads/${sample}.2.trimmed.fastq.gz fastq_files/${sample}_L2_1.fq.gz fastq_files/${sample}_L2_2.fq.gz
done

#Running fastqc on the filtered reads
fastqc -t 12 fastq_files/trimmed_reads/*.fq.gz 
multiqc fastq_files/trimmed_reads/
#Check that adapter have been trimmed

#Running salmon against transcriptome
mkdir transcript_quant
cat sample_list.txt | while read sample; 
	do salmon quant -i /media/kilian/OS/References/mm10_salmon/ -l A -1 fastq_files/${sample}_L1_1.fq.gz -2 fastq_files/${sample}_L1_2.fq.gz --validateMappings -o transcript_quant/${sample}_quant --thread 8
done

### STAR ALIGNMENT for TE analysis
#STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm10 --genomeFastaFiles mm10_STAR/mm10.fa  #--sjdbGTFfile mm10_STAR/mm10.refGene.gtf
mkdir BAM_files_multi_trimmed
cat sample_list.txt | while read sample; do
	STAR --runThreadN 8 \
	--readFilesIn fastq_files/trimmed_reads/${sample}.1.trimmed.fastq.gz fastq_files/trimmed_reads/${sample}.2.trimmed.fastq.gz \
	--genomeDir /media/kilian/OS/STAR_index_mm10 \
	--outSAMtype BAM SortedByCoordinate  \
	--runMode alignReads \
	--outFilterMultimapNmax 100 \
	--outMultimapperOrder Random \
	--winAnchorMultimapNmax 100 \
	--outFilterMismatchNmax 3 \
	--alignEndsType EndToEnd \
	--alignIntronMax 1 \
	--alignMatesGapMax 350 \
	--outFileNamePrefix BAM_files_multi_trimmed/${sample} \
	--readFilesCommand zcat 
done

#Index files
cat sample_list.txt | while read sample; do
	samtools index BAM_files_multi_trimmed/${sample}Aligned.sortedByCoord.out.bam
done

#Output metrics
mkdir Picard_metrics
cat sample_list.txt | while read sample; do
picard CollectRnaSeqMetrics \
    I=BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam \
    O=Picard_metrics/${sample}_output_rna_metrics.txt \
    REF_FLAT=/media/kilian/OS/References/GTF_files_TEtranscript/mm10_refFlat.txt \
    STRAND_SPECIFICITY=NONE \
    R=/media/kilian/OS/References/GTF_files_TEtranscript/mm10.fa
done


#Running TEtranscript
#Download required files (for mouse)
mkdir GTF_files_TEtranscript
cd GTF_files_TEtranscript
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm38_GENCODE_rmsk_TE.gtf.gz.download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz

#Multimapping
#TEtranscript
mkdir TETranscripts_multi
TEtranscripts -t  \
                     -c 
                     --GTF /media/kilian/OS/GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE /media/kilian/OS/GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf \
                     --mode multi \
                     --outdir TETranscripts_multi \
                     --sortByPos

#TElocal                   
cat sample_list.txt | while read sample; do
TElocal -b BAM_files_multi_trimmed/${sample}Aligned.sortedByCoord.out.bam \
        --GTF /media/kilian/OS/GTF_files_TEtranscript/mm10.refGene.gtf \
        --TE /media/kilian/OS/GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf.locInd \
        --project ${sample} \
	--mode multi \
	--sortByPos
done


#Unique mapping 
mkdir TETranscripts_uniq
TEtranscripts -t BAM_files/DE12NGSUKBR137309_1.final.bam BAM_files/DE55NGSUKBR137311_1.final.bam \
                     -c BAM_files/DE39NGSUKBR137308_1.final.bam BAM_files/DE82NGSUKBR137310_1.final.bam \
                     --GTF GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf \
                     --mode uniq \
                     --outdir TETranscripts_uniq

#TElocal
cat sample_list.txt | while read sample; do
TElocal -b BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam \
        --GTF /media/kilian/OS/GTF_files_TEtranscript/mm10.refGene.gtf \
        --TE /media/kilian/OS/GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf.locInd \
        --project ${sample} \
	--mode uniq \
	--sortByPos
done



#Quantification of TEs with Stringtie
mkdir stringtie_output_TE_RMSK
cat sample_list.txt | while read sample; 
	do stringtie -e -B -p 8 -G reference/GRCm38_GENCODE_rmsk_TE.gtf -o stringtie_output_TE_RMSK/${sample}/${sample}.gtf BAM_files/${sample}Aligned.sortedByCoord.out.bam 
done

mkdir stringtie_output_TE_RMSK
cat sample_list.txt | while read sample; 
	do stringtie -e -B -p 8 -G reference/GRCm38_GENCODE_rmsk_TE.gtf -o stringtie_output_TE_RMSK/${sample}/${sample}.gtf BAM_files/${sample}Aligned.sortedByCoord.out.bam 
done

#Running Python script to extract raw counts
cd stringtie_output_TE_RMSK
python /media/kilian/OS/Melani_THP1_August_2023/scripts_analysis/prepDE.py3


#Trying with full length bed file
scripts/bed2gtf.py -s reference/L1base_FL_mm10.bed  

