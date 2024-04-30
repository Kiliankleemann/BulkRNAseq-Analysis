#Make fastq_folder
mkdir fastq_files
mv *.fq.gz fastq_files

#Detect all samples in fastq_files folder
echo 'Detected fastq samples'
find ./fastq_files -name "*.fq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#Running fastq for all files
mkdir QC
mv QZ.zip QC
unzip QC.zip

#Create multiqc report
multiqc QC/ -o QC/

#run cutadapt on all samples 
mkdir fastq_files/trimmed_reads
cat sample_list.txt | while read sample; do
echo $sample
	cutadapt --cores 12 --minimum-length 15 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o fastq_files/trimmed_reads/${sample}.trimmed.fq.gz fastq_files/${sample}.fq.gz --poly-a
done

#Running fastqc on the filtered reads
mkdir QC/trimmed_reads
cat sample_list.txt | while read sample; 
	do fastqc -t 8 fastq_files/trimmed_reads/${sample}.trimmed.fq.gz -o QC/trimmed_reads
done

#Check quality of samples
multiqc QC/trimmed_reads -o QC/trimmed_reads

#Running salmon against transcriptome
mkdir transcript_quant
cat sample_list.txt | while read sample; 
	do salmon quant -i /media/kilian/OS/mm10_salmon/ -l A -r fastq_files/trimmed_reads/${sample}.trimmed.fq.gz --validateMappings -o transcript_quant/${sample}_quant --threads 12
done


