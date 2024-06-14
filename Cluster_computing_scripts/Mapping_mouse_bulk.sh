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
conda install tetranscripts

#Detect all samples in fastq_files folder
echo 'Detected fastq samples'
find ./fastq_files -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#Running fastq for all files
mkdir fastq_files/QC
cd fastq_files
fastqc *.fastq.gz

#Move all fastqc reports into the QC folder 
mv *fastqc* QC/


#Check quality of samples
cd QC
multiqc .
cd ..

#run fastp on all samples 
mkdir fastq_files/fastp_sorted
cat sample_list.txt | while read sample; 
	do fastp -i fastq_files/${sample}.fq.gz -o fastq_files/fastp_sorted/${sample}.fq.gz
done

#Running fastqc on the filtered reads
mkdir fastq_files/fastp_sorted/fastqc_reports
cat sample_list.txt | while read sample; 
	do fastqc -t 4 fastq_files/fastp_sorted/${sample}.fq.gz -o fastq_files/fastp_sorted/fastqc_reports/
done

#Check quality of samples
cd fastq_files/fastp_sorted/fastqc_reports/
multiqc .
cd ..

#Running bowtie (alignment)
#First running sortmeRNA (only on linux)
mkdir SAM_files
cat sample_list.txt | while read sample; do
	bowtie2 -N 1 -X 1000 -p 8  -x ~/Desktop/mm10/mm10 -U fastq_files/fastp_sorted/${sample}.fq.gz  -S SAM_files/${sample}.sam 
done

mkdir BAM_files
echo "Sorting SAM files and output BAM"
cat sample_list.txt | while read sample; do
	samtools sort SAM_files/${sample}.sam -o BAM_files/${sample}.final.bam
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
TEtranscripts -t  BAM_files/DE12NGSUKBR137309_1.final.bam \
                     -c  BAM_files/DE39NGSUKBR137308_1.final.bam \
                     --GTF GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf \
                     --mode multi \
                     --outdir TETranscripts_multi

#Unique mapping 
TEtranscripts -t  BAM_files/DE12NGSUKBR137309_1.final.bam \
                     -c  BAM_files/DE39NGSUKBR137308_1.final.bam \
                     --GTF GTF_files_TEtranscript/mm10.refGene.gtf \
                     --TE GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf \
                     --mode uniq \
                     --outdir TETranscripts_multi



#####----- Script for the server ----#####
#Running allignement with STAR


#Random mode
mkdir BAM_files
cat sample_list.txt | while read sample; do
	STAR --runThreadN 4 \
	--readFilesIn fastq_files/${sample}.fastq.gz \
	--genomeDir indexes/ \
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

#TEtranscripts
singularity pull tetranscripts.sif docker://mhammelllab/tetranscripts:latest
singularity exec tetranscripts.sif TEtranscripts -t <treatment sample> -c <control sample> --GTF <genic-GTF-file> --TE <TE-GTF-file>




#Quantification with ExplorATE
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.out.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz


CONDA_SUBDIR=osx-64 conda create -n ExplorATE python
  conda activate ExplorATE
  conda config --env --set subdir osx-64


bash ExplorATE_shell_script/ExplorATE mo \
	-p 12 \
	-f mm10.fa \
	-s /Users/kiliankleemann/opt/anaconda3/envs/SALMON \
	-b /Users/kiliankleemann/opt/anaconda3/envs/SALMON \
	-g mm10.refGene.gtf \
	-r mm10.fa.out \
	-e 'se' \
	-l fastq_files/ \
	-o out_mo


#Running salmon against transcriptome
mkdir transcript_quant
cat sample_list.txt | while read sample; 
	do salmon quant -i /Volumes/OhneTitel/mm10 -l A -r fastq_files/fastp_sorted/${sample}.fq.gz --validateMappings -o transcript_quant/${sample}_quant
done


