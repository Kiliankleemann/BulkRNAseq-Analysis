#Alignment with Minimap2
#Issue with x64 M1 chip !! do not copy the notes into commad line
CONDA_SUBDIR=osx-64 conda create -n Bulk_RNA_TE_pipeline 
conda activate Bulk_RNA_TE_pipeline
conda env config vars set CONDA_SUBDIR=osx-64
conda deactivate
conda activate Bulk_RNA_TE_pipeline
conda install multiqc
conda install fastqc
conda install bowtie2
conda install samtools
conda install salmon
conda install fastp

#Downgrading python for the environment
conda deactivate 
conda activate Bulk_RNA_TE_pipeline
conda install python=3.7
conda install tetranscripts

#Detect all samples in fastq_files folder
echo 'Detected fastq samples'
find ./fastq_files -name "*.fq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#Running fastq for all files
mkdir fastq_files/QC
cd fastq_files
fastqc *.fq.gz -t 8

#Move all fastqc reports into the QC folder 
mv *fastqc* QC/

#Check quality of samples
cd QC
multiqc .
cd ..

#run fastp on all samples 
mkdir fastq_files/fastp_sorted
cat sample_list.txt | while read sample; 
	do fastp -i fastq_files/${sample}.fq.gz -o fastq_files/fastp_sorted/${sample}.fq.gz --thread 8
done


mkdir fastq_files/trimmed_reads
cat sample_list.txt | while read sample; do
echo $sample
cutadapt --cores 8 --minimum-length 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o fastq_files/trimmed_reads/${sample}.trimmed.fq.gz fastq_files/${sample}.fq.gz --poly-a
done


#Running fastqc on the filtered reads
#mkdir fastq_files/fastp_sorted/fastqc_reports
#cat sample_list.txt | while read sample; 
#	do fastqc -t 8 fastq_files/fastp_sorted/${sample}.fq.gz -o fastq_files/fastp_sorted/fastqc_reports/
#done

#Check quality of samples
#cd fastq_files/fastp_sorted/fastqc_reports/
#multiqc .
#cd ..


#Running salmon against transcriptome
mkdir transcript_quant
cat sample_list.txt | while read sample; 
	do salmon quant -i /media/kilian/OS/mm10_salmon/ -l A -r fastq_files/fastp_sorted/${sample}.fq.gz --validateMappings -o transcript_quant/${sample}_quant --threads 8
done

### STAR ALIGNMENT
STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm10 --genomeFastaFiles mm10_STAR/mm10.fa  --sjdbGTFfile mm10_STAR/mm10.refGene.gtf

cat sample_list.txt | while read sample; do
	STAR --runThreadN 12 \
	--readFilesIn fastq_files/${sample}.fq.gz \
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
	--outFileNamePrefix ${sample} \
	--readFilesCommand zcat 
done


cat sample_list.txt | while read sample; do
	samtools index BAM_files/${sample}.final.bam
done

#Running TEtranscript
#Download required files
mkdir GTF_files_TEtranscript
cd GTF_files_TEtranscript
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm38_GENCODE_rmsk_TE.gtf.gz.download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz

#Multimapping
mkdir TETranscripts_multi
TEtranscripts -t  BAM_multi_trimmed/DE14NGSUKBR140280_1Aligned.sortedByCoord.out.bam \
                     -c BAM_multi_trimmed/DE26NGSUKBR140258_1Aligned.sortedByCoord.out.bam \
                     --GTF /media/kilian/OS/GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE /media/kilian/OS/GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf \
                     --mode multi \
                     --outdir TETranscripts_multi
 
                     
                     
                     
#Unique mapping 
mkdir TETranscripts_uniq
TEtranscripts -t BAM_files/DE12NGSUKBR137309_1.final.bam BAM_files/DE55NGSUKBR137311_1.final.bam \
                     -c BAM_files/DE39NGSUKBR137308_1.final.bam BAM_files/DE82NGSUKBR137310_1.final.bam \
                     --GTF GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf \
                     --mode uniq \
                     --outdir TETranscripts_uniq



#Running salmon against transcriptome
mkdir transcript_quant
cat sample_list.txt | while read sample; 
	do salmon quant -i /Volumes/OhneTitel/mm10 -l A -r fastq_files/fastp_sorted/${sample}.fq.gz --validateMappings -o transcript_quant/${sample}_quant
done



#Quantification of TEs with Stringtie
mkdir stringtie_output_TE_RMSK
cat sample_list.txt | while read sample; 
	do stringtie -e -B -p 8 -G reference/GRCm38_GENCODE_rmsk_TE.gtf -o stringtie_output_TE_RMSK/${sample}/${sample}.gtf BAM_files/${sample}Aligned.sortedByCoord.out.bam 
done

#Running Python script to extract raw counts
cd stringtie_output_TE_RMSK
python /media/kilian/OS/Melani_THP1_August_2023/scripts_analysis/prepDE.py3


#Trying with full length bed file
scripts/bed2gtf.py -s reference/L1base_FL_mm10.bed  

