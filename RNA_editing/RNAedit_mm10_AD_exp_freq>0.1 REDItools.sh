conda create -n nature_protocol python=2.7 anaconda
conda activate nature_protocol
$ conda install -n nature_protocol bcftools==1.9
$ conda install -n nature_protocol bedtools==2.28.0
$ conda install -n nature_protocol bzip2==1.0.6
$ conda install -n nature_protocol bwa==0.7.17
$ conda install -n nature_protocol bx-python==0.8.2
$ conda install -n nature_protocol fastp==0.20.0
$ conda install -n nature_protocol fastqc==0.11.8
$ conda install -n nature_protocol fisher==0.1.4 (optional)
$ conda install -n nature_protocl git==2.21.0
$ conda install -n nature_protocol gmap==2018.07.04
$ conda install -n nature_protocol htslib==1.9
$ conda install -n nature_protocol libdeflate==1.0
$ conda install -n nature_protocol numpy==1.16.2
$ conda install -n nature_protocol pysam==0.15.2
$ conda install -n nature_protocol rseqc==2.6.4
$ conda install -n nature_protocol samtools==1.9
$ conda install -n nature_protocol scipy==1.2.1
$ conda install -n nature_protocol star==2.7.0f
$ conda install -n nature_protocol wget==1.20.1
#########################
##install REDItools from github Now using REDItools2
cd /media/xinwei/6886A4590BAF2BC0/SharedFolder/RNA_seq1/N2412481_80-1540194018_2024-06-05/240603-A00199A
mkdir rna_editing_protocol
cd rna_editing_protocol
git clone https://github.com/BioinfoUNIBA/REDItools.git
cd REDItools
##check if installed
python -c 'import pysam' #don't work
conda list | grep pysam #work
##conda install pysam==0.15.2
#cd main #test
#python REDItoolDnaRna.py -h # can work

####################################################download#########################################3
#install pblat.git for multithread
cd rna_editing_protocol
git clone https://github.com/icebert/pblat.git
cd pblat/
make
cd ..

#Download and unzip the human reference genome hg19 in FASTA format ● Timing ~1 min:
mkdir genome_hg19
cd genome_hg19/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip GRCh37.primary_assembly.genome.fa.gz
cd ..

#download and unzip hg19 ##or mm10 from
#https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCM38.genome.fa.gz for Release M35 (GRCM38)
#https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz for Release 46 (GRCh38.p14)
mkdir Gencode_annotation
cd Gencode_annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gtf.gz
gunzip gencode.v30lift37.annotation.gtf.gz
cd ..

#ownload and unzip hg19 RefSeq annotations in bed format for strand detection ● Timing ~30 s
mkdir Strand_detection
cd Strand_detection
#wget --no-check-certiﬁcate https://sourceforge.net/projects/rseqc/ﬁles/BED/Human_Homo_sapiens/hg19_RefSeq.bed.gz
wget https://sourceforge.net/projects/rseqc/ﬁles/BED/Human_Homo_sapiens/hg19_RefSeq.bed.gz
gunzip hg19_RefSeq.bed.gz
cd ..
#ownload and unzip mm10 RefSeq annotations in bed format for strand detection ● Timing ~30 s
mkdir Strand_detection
cd Strand_detection
#wget --no-check-certiﬁcate https://sourceforge.net/projects/rseqc/ﬁles/BED/Human_Homo_sapiens/hg19_RefSeq.bed.gz
# wget https://sourceforge.net/projects/rseqc/ﬁles/BED/Mouse_Mus_musculus/mm10_RefSeq.bed.gz
gunzip mm10_RefSeq.bed.gz
cd ..

#Download and unzip RepeatMasker annotations ● Timing ~1 min
mkdir rmsk
cd rmsk
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
gunzip rmsk.txt.gz
cd ..

#Download and unzip dbSNP annotations ● Timing ~23 min
mkdir snp142
cd snp142
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp142.txt.gz
gunzip snp142.txt.gz
cd ..

#Download and unzip REDIportal annotations ● Timing ~1 min
mkdir rediportal
cd rediportal
wget http://srv00.recas.ba.infn.it/webshare/rediportalDownload/table1_full.txt.gz
gunzip table1_full.txt.gz
cd ..

#########################################continue Preparation for required data#################
#build a reference genome index for BWA 
mkdir genome_mm10
cd genome_mm10/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/GRCm38.p4.genome.fa.gz
gunzip GRCm38.p4.genome.fa.gz
bwa index GRCm38.p4.genome.fa
cd..
cd rmsk
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
cd snp142
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/snp142.txt.gz

$ mkdir Gencode_annotation
$ cd Gencode_annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.vM35.annotation.gtf.gz
gunzip gencode.vM35.annotation.gtf.gz
cd ..

#build a reference genome index for STAR
mkdir STAR
cd STAR
mkdir STAR_genome_index_ucsc
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_hg19 --genomeFastaFiles genome_hg19/GRCh37.primary_assembly.genome.fa --sjdbGTFfile Gencode_annotation/gencode.v30lift37.annotation.gtf
#STAR --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm10 --genomeFastaFiles mm10_STAR/mm10.fa --sjdbGTFfile mm10_STAR/mm10.refGene.gtf

#prepare RepeatMasker annotations for REDItools
cd rmsk
awk 'OFS="\t"{print $6,"rmsk_mm10",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}' rmsk.txt > rmsk.gtf
sort -k1,1 -k4,4n rmsk.gtf > rmsk.sorted.gtf
bgzip rmsk.sorted.gtf
tabix -p gff rmsk.sorted.gtf.gz
cd ..

#prepare dbSNP annotations for REDItools
cd snp142
awk 'OFS="\t"{if ($11=="genomic" && $12=="single") print $2,"ucsc_snp142_mm10","snp",$4,$4,".",$7,".","gene_id \""$5"\"; transcript_id \""$5"\";"}' snp142.txt > snp142.gtf
sort -k1,1 -k4,4n snp142.gtf > snp142.sorted.gtf
bgzip snp142.sorted.gtf
tabix -p gff snp142.sorted.gtf.gz
cd ..

#prepare splice sites annotations for REDItools
cd Gencode_annotation
gtf_splicesites gencode.vM35.annotation.gtf > splicesites
awk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' splicesites > gencode.vM35.splicesites.txt
cd ..

#prepare REDIportal annotations for REDItools and extract recoding events
cd rediportal
awk 'OFS="\t"{sum+=1; print $1,"rediportal","ed",$2,$2,".",$5,".","gene_id \""sum"\"; transcript_id \""sum"\";"}' TABLE1_mm10.txt > atlas.gtf
bgzip atlas.gtf
tabix -p gff atlas.gtf.gz
python ../REDItools/accessory/rediportal2recoding.py TABLE1_mm10.txt > atlas_recoding.gff
sort -V -k1,1 -k4,4n atlas_recoding.gff > srtd_atlas_recoding.gff
bgzip srtd_atlas_recoding.gff
tabix -p gff srtd_atlas_recoding.gff.gz
cd ..

#create the nochr file for REDItools
cd genome_mm10/
grep ">" GRCm38.p4.genome.fa | awk '{if (substr($1,1,3)==">GL") print $2}' > nochr
cd ..

#index  the reference genome for REDItools
cd genome_mm10/
samtools faidx GRCm38.p4.genome.fa

##################################Procedure 1: RNA editing detection in the NA12878 cell line#####################
#Edit fastq  ending of file and symbol to cut the name 
cd RNA_mm10
find ./fastq_files -name "*.fastq.gz" -type f -exec basename "{}" \; |  cut -d'_' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#Running fastq for all files
mkdir QC
fastqc fastq_files/*.fastq.gz -t 12 -o QC

#Make quality overview
multiqc QC/ -o QC

#run cutadapt on all samples 
mkdir fastq_files/trimmed_reads
cat sample_list.txt | while read sample; do
	echo $sample
    fastp -i fastq_files/${sample}_1.fastq.gz -I fastq_files/${sample}_2.fastq.gz -o fastq_files/trimmed_reads/${sample}.1.trimmed.fq.gz -O fastq_files/trimmed_reads/${sample}.2.trimmed.fq.gz -q 25 -u 10 -l 50 -y -x -w 4
done
#Running fastqc on the filtered reads
mkdir QC/trimmed_reads
fastqc -t 12 fastq_files/trimmed_reads/*.fq.gz -o QC/trimmed_reads
multiqc QC/trimmed_reads -o QC/trimmed_reads

#Running salmon against transcriptome
#mkdir transcript_quant
#cat sample_list.txt | while read sample; 
#	do salmon quant -i /media/kilian/OS/References/mm10_salmon/ -l A -1 fastq_files/${sample}_L1_1.fq.gz -2 fastq_files/${sample}_L1_2.fq.gz --validateMappings -o transcript_quant/${sample}_quant --thread 8
#done

### STAR ALIGNMENT for TE analysis
#STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_mm10 --genomeFastaFiles mm10_STAR/mm10.fa

#cat sample_list.txt | while read sample; do
#	STAR --runThreadN 12 \
#	--readFilesIn fastq_files/trimmed_reads/${sample}.1.trimmed.fq.gz fastq_files/trimmed_reads/${sample}.2.trimmed.fq.gz  \
#	--genomeDir /media/xinwei/6886A4590BAF2BC0/reference/STAR_index_hg19 \
#	--outSAMtype BAM SortedByCoordinate  \
#	--runMode alignReads \
#	--outFileNamePrefix BAM_files_sorted/${sample} \
#	--readFilesCommand zcat \
#    --outTmpDir /media/xinwei/6886A4590BAF2BC0/SharedFolder/STAR_tmp
#done
cat sample_list.txt | while read sample; do
STAR --runThreadN 12 \
     --genomeDir /media/xinwei/6886A4590BAF2BC0/reference/STAR_index_mm10 \
     --genomeLoad NoSharedMemory \
     --outFileNamePrefix BAM_files_sorted/${sample} \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMstrandField intronMotif \
     --outSAMattributes All \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --runMode alignReads \
     --outFilterMultimapNmax 1 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --readFilesIn  fastq_files/trimmed_reads/${sample}.1.trimmed.fq.gz fastq_files/trimmed_reads/${sample}.2.trimmed.fq.gz  \
     --outTmpDir /media/xinwei/6886A4590BAF2BC0/SharedFolder/STAR_tmp
done
#Index files
cat sample_list.txt | while read sample; do
	samtools index BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam
done
cd ..
#Detection of the strand orientation of RNAseq reads
cd Strand_detection/
cat sample_list.txt | while read sample; do
infer_experiment.py -r mm10_RefSeq.bed -s 2000000 -i ../RNA_mm10/BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam > ../RNA_mm10/${sample}strandinfo.txt
done
cd ..



#REDItools2 running
mkdir Reditool_mm10
cd Reditool_mm10
cat sample_list.txt | while read sample; do
SOURCE_BAM_FILE="../RNA_mm10/BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam"
REFERENCE="../genome_mm10/GRCm38.p4.genome.fa"
SIZE_FILE="../genome_mm10/GRCm38.p4.genome.fa.fai"

NUM_CORES=8
mkdir test_results_${sample}
OUTPUT_FILE="test_results_${sample}/output/parallel_table.txt.gz"
TEMP_DIR="test_results_${sample}/temp"
COVERAGE_FILE="test_results_${sample}/coverage/${sample}Aligned.sortedByCoord.out.cov"
COVERAGE_DIR="test_results_${sample}/coverage/"

../REDItools2.0/extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
mpirun -np $NUM_CORES ../REDItools2.0/src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR
mv test_results_${sample}/tempfiles.txt test_results_${sample}/temp/files.txt
mv test_results_${sample}/tempgroups.txt test_results_${sample}/temp/groups.txt
mv test_results_${sample}/temptimes.txt test_results_${sample}/temp/times.txt
../REDItools2.0/merge.sh $TEMP_DIR $OUTPUT_FILE $NUM_CORES
done

##Filtering of DNA-RNA variants
#TEST for REDItools2
cat sample_list.txt | while read sample; do
cd test_results_${sample}/output
gunzip parallel_table.txt.gz
awk 'FS="\t" {if ($8!="-" && $5>=10 && $13=="-") print}' parallel_table.txt > parallel_table.txt_all_chr.out
#Annotate positions using RepeatMasker and dbSNP annotations
python ../../../REDItools/accessory/AnnotateTable.py -a /media/xinwei/2TPassport/RNA_editing/mm10/rmsk/rmsk.sorted.gtf.gz -n rmsk -i parallel_table.txt_all_chr.out -o parallel_table.txt_all_chr.out.rmsk -u
python ../../../REDItools/accessory/AnnotateTable.py -a /media/xinwei/2TPassport/RNA_editing/mm10/snp142/snp142.sorted.gtf.gz -n snp142 -i parallel_table.txt_all_chr.out.rmsk -o parallel_table.txt_all_chr.out.rmsk.snp -u
#Create a first set of positions selecting sites supported by at least five RNAseq reads and a single mismatch
python ../../../REDItools/accessory/selectPositions.py -i parallel_table.txt_all_chr.out.rmsk.snp -c 5 -v 1 -f 0.0 -o parallel_table.txt_all_chr.out.rmsk.snp.sel1
#Create a second set of positions selecting sites supported by ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1
python ../../../REDItools/accessory/selectPositions.py -i parallel_table.txt_all_chr.out.rmsk.snp -c 10 -v 3 -f 0.1 -o parallel_table.txt_all_chr.out.rmsk.snp.sel2

#Select ALU sites from the first set of positions
awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-"&& $8!="-" && $13=="-") print}' parallel_table.txt_all_chr.out.rmsk.snp.sel1 > parallel_table.txt_all_chr.out.rmsk.snp.alu

#Select REP NON ALU sites from the second set of positions, excluding sites in Simple repeats or Low complexity regions
awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $9>=0.1) print}' parallel_table.txt_all_chr.out.rmsk.snp.sel2 > parallel_table.txt_all_chr.out.rmsk.snp.nonalu

#Select NON REP sites from the second set of positions
awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $9>=0.1) print}' parallel_table.txt_all_chr.out.rmsk.snp.sel2 > parallel_table.txt_all_chr.out.rmsk.snp.nonrep

#Annotate ALU, REP NON ALU and NON REP sites using known editing events from REDIportal
python ../../../REDItools/accessory/AnnotateTable.py -a ../../../rediportal/atlas.gtf.gz -n ed -k R -c 1 -i parallel_table.txt_all_chr.out.rmsk.snp.alu -o parallel_table.txt_all_chr.out.rmsk.snp.alu.ed -u
python ../../../REDItools/accessory/AnnotateTable.py -a ../../../rediportal/atlas.gtf.gz -n ed -k R -c 1 -i parallel_table.txt_all_chr.out.rmsk.snp.nonalu -o parallel_table.txt_all_chr.out.rmsk.snp.nonalu.ed -u
python ../../../REDItools/accessory/AnnotateTable.py -a ../../../rediportal/atlas.gtf.gz -n ed -k R -c 1 -i parallel_table.txt_all_chr.out.rmsk.snp.nonrep -o parallel_table.txt_all_chr.out.rmsk.snp.nonrep.ed -u
cd ..
cd ..
done

#Extract known editing events from ALU, REP NON ALU and NON REP sites
#mv parallel_table.txt_all_chr.out.rmsk.snp.alu.ed alu../../../Alignment/${sample}Aligned.sortedByCoord.out.bam -f ../../../genome_mm10/GRCm38.p4.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w /media/xinwei/2TPassport/RNA_editing/mm10/Gencode_annotation/gencode.vM35.splicesites.txt -k ../../../genome_mm10/nochr --reads -R --addP -o first
#done
cat sample_list.txt | while read sample; do
cd test_results_${sample}/output
mv parallel_table.txt_all_chr.out.rmsk.snp.alu.ed alu
mv parallel_table.txt_all_chr.out.rmsk.snp.nonalu.ed nonalu
mv parallel_table.txt_all_chr.out.rmsk.snp.nonrep.ed nonrep
cat alu nonalu nonrep > alu-nonalu-nonrep
awk 'FS="\t" {if ($19=="ed") print}' alu-nonalu-nonrep > knownEditing

#Convert editing candidates in REP NON ALU and NON REP sites in GFF format for further filtering
cat nonalu nonrep > nonalu-nonrep
awk 'FS="\t" {if ($19!="ed") print}' nonalu-nonrep > pos.txt
python ../../../REDItools/accessory/TableToGFF.py -i pos.txt -s -t -o pos.gff

#Convert editing candidates in ALU sites in GFF format for further filtering:
awk 'FS="\t" {if ($19!="ed") print}' alu > posalu.txt
#Launch REDItoolDnaRna.py on ALU sites using stringent criteria to recover potential editing candidates:
#cat /media/xinwei/6886A4590BAF2BC0/RNA_editing/RNA_mm10/BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam -f ../../../genome_mm10/GRCm38.p4.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w /media/xinwei/2TPassport/RNA_editing/mm10/Gencode_annotation/gencode.vM35.splicesites.txt -k ../../../genome_mm10/nochr --reads -R --addP -o first
#doneting/rna_editing_protocol/RNA_mm10sample_list.txt | while read sample; do
cat sample_list.txt | while read sample; do
cd test_results_${sample}/output
python ../../../REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 12 -i /BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam -f ../../../genome_mm10/GRCm38.p4.genome.fa -c 5,5 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 1 -n 0.0 -e -T posalu.sorted.gff.gz -w /media/xinwei/2TPassport/RNA_editing/mm10/Gencode_annotation/gencode.vM35.splicesites.txt -k ../../../genome_mm10/nochr -R -o firstalu

#Launch REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria to recover RNAseq reads harboring reference mismatches:
#cat sample_list.txt | while read sample; do
python ../../../REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 12 -i /BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam -f ../../../genome_mm10/GRCm38.p4.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w /media/xinwei/2TPassport/RNA_editing/mm10/Gencode_annotation/gencode.vM35.splicesites.txt -k ../../../genome_mm10/nochr --reads -R --addP -o first
cd ..
cd ..
done
#Launch pblat on RNAseq reads harboring reference mismatches from Step 22 and select multi-mapping reads:
cat sample_list.txt | while read sample; do
find ./test_results_${sample}/output/first -name "DnaRna_*" -maxdepth 1 -exec basename "{}" \; |  cut -d '_' -f2 > ID.txt
cat ID.txt
done


######extract and merge ID.txt to get ID_all.txt
paste Reditools_mm10/ID_all.txt sample_list.txt | while IFS=$'\t' read -r ID sample; do
  cd test_results_${sample}/output
  ../../../pblat/pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 ../../../genome_mm10/GRCm38.p4.genome.fa first/DnaRna_${ID}/outReads_${ID} reads.psl
  python ../../../REDItools/accessory/readPsl.py reads.psl badreads.txt

  #Extract RNAseq reads harboring reference mismatches from Step 22 and remove duplicates:
  sort -k1,1 -k2,2n -k3,3n first/DnaRna_${ID}/outPosReads_${ID} | mergeBed > bed
  cd ..
  cd ..
done
#sample=AD_Ect_F1
cat sample_list.txt | while read sample; do
cd test_results_${sample}/output
samtools view -@ 12 -L bed -h -b /BAM_files_sorted/${sample}Aligned.sortedByCoord.out.bam > ${sample}_bed.bam
samtools sort -@ 12 -n ${sample}_bed.bam -o ${sample}_bed_ns.bam
samtools fixmate -@ 12 -m ${sample}_bed_ns.bam ${sample}_bed_ns_fx.bam
samtools sort -@ 12 ${sample}_bed_ns_fx.bam -o ${sample}_bed_ns_fx_st.bam
samtools markdup -r -@ 12 ${sample}_bed_ns_fx_st.bam ${sample}_bed_dedup.bam
samtools index ${sample}_bed_dedup.bam

#Re-run REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria,deduplicated reads and mis-mapping info:
python ../../../REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 4 -i ${sample}_bed_dedup.bam -f ../../../genome_mm10/GRCm38.p4.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w /media/xinwei/2TPassport/RNA_editing/mm10/Gencode_annotation/gencode.vM35.splicesites.txt -R -k ../../../genome_mm10/nochr -b badreads.txt --rmIndels -o second

#Collect filtered ALU, REP NON ALU and NON REP sites:
python ../../../REDItools/NPscripts/collect_editing_candidates.py
sort -k1,1 -k2,2n editing.txt > editing_sorted.txt

#Inspect the distribution of editing candidates to look at A-to-I enrichment:
python ../../../REDItools/NPscripts/get_Statistics.py
cd ..
cd ..
done
#sample=AD_Ect_F1
cat sample_list.txt | while read sample; do
cd /media/xinwei/6886A4590BAF2BC0/RNA_editing/Rworking
mkdir tables
cd tables
cd /media/xinwei/6886A4590BAF2BC0/RNA_editing/rna_editing_protocol/Reditools_mm10/test_results_${sample}/output
cp editing_sorted.txt /media/xinwei/6886A4590BAF2BC0/RNA_editing/Rworking/tables/${sample}.txt
done

#extract statistics file
cat sample_list.txt | while read sample; do
cd /media/xinwei/6886A4590BAF2BC0/RNA_editing/Rworking
mkdir tables
cd tables
cd /media/xinwei/6886A4590BAF2BC0/RNA_editing/rna_editing_protocol/Reditools_mm10/test_results_${sample}/output
cp editingStats.txt /media/xinwei/6886A4590BAF2BC0/RNA_editing/Rworking/tables/${sample}_Stats.txt
done
#Create a comma-separated information file
#Sample,Status
#SRR3306823,DIS
#SRR3306824,DIS
#SRR3306825,DIS
#SRR3306826,DIS
#SRR3306827,DIS
#SRR3306828,DIS
#SRR3306829,DIS
#SRR3306830,CTRL
#SRR3306831,CTRL
#SRR3306832,CTRL
#SRR3306833,CTRL
#SRR3306834,CTRL
conda activate DE ###step in QEedit in github
#Detect differential RNA editing:
python /media/xinwei/6886A4590BAF2BC0/RNA_editing/QEdit-master/scripts/sample_status_file_creator.py sample_information.csv CTRL AD
#python /media/xinwei/6886A4590BAF2BC0/RNA_editing/QEdit-master/scripts/sample_status_file_creator.py sample_information_file.csv DIS CTRL
#python /media/xinwei/6886A4590BAF2BC0/RNA_editing/QEdit-master/scripts/sample_path_folder_creator.py DIS_vs_CTRL.sif
python /media/xinwei/6886A4590BAF2BC0/RNA_editing/QEdit-master/scripts/sample_path_folder_creator.py CTRL_vs_AD.sif

#python /media/xinwei/6886A4590BAF2BC0/RNA_editing/rna_editing_protocol/REDItools/accessory/get_DE_events.py -cpval 2 -input_file CTRL_vs_Apobec.sif -sig no -linear
#python get_DE_events_TC.py -cpval 2 -input_file CTRL_vs_AD.sif -sig no -linear

python get_DE_events_filter_same.py -cpval 2 -input_file CTRL_vs_AD.sif -sig yes -linear

python get_DE_events_filter_same.py -cpval 2 -input_file CTRL_vs_Apobec.sif -sig yes -linear
