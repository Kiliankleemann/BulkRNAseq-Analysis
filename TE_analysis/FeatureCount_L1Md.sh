# Download RepeatMasker tracks for LINE1
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
gunzip rmsk.txt.gz

### COUNTING WITH FEATURECOUNTS 

#MOUSE
#Filter out only L1Md_T elements
awk '$0 ~ /gene_id "L1Md/' GRCm38_GENCODE_rmsk_TE.gtf > GRCm38_GENCODE_L1Md.gtf
awk '$0 ~ /gene_id "L1Md/' mm10_rmsk_TE.gtf > mm10_L1Md.gtf

#filter out 
bedtools subtract -A -a mm10_L1Md.gtf -b mm10.transcripts_filtered.gtf > mm10_L1Md_filtered.gtf


#Modify the gtf file to suit featureCounts counting!
#!/bin/bash
# Input and output file paths
INPUT_GTF="GRCm38_GENCODE_L1Md.gtf"
OUTPUT_GTF="GRCm38_GENCODE_L1Md_unique.gtf"

# Initialize a counter
counter=1

# Create or overwrite the output file
> "$OUTPUT_GTF"

# Read the input GTF line by line
while IFS= read -r line; do
    # Skip comment lines
    if [[ $line == \#* ]]; then
        echo "$line" >> "$OUTPUT_GTF"
    else
        # Append a unique counter to gene_id
        modified_line=$(echo "$line" | sed -E "s/(gene_id \"L1Md)(\";)/\1_$counter\2/")
        echo "$modified_line" >> "$OUTPUT_GTF"
        counter=$((counter + 1))
    fi
done < "$INPUT_GTF"

echo "Modified GTF file saved to $OUTPUT_GTF"


#Filter out gene overlapping LINE1
mkdir featureCounts_out
featureCounts -p -B -C \
  -a /media/kilian/OS/References/GTF_files_TEtranscript/GRCm38_GENCODE_L1Md.gtf \
  -o featureCounts_out/DE07NGSUKBR147099Aligned_featurecounts_output.txt \
/media/kilian/DATA/Chronic_Acute/BMDM_IFNa4KI_PE_241217_P2024-266-RNA/BAM_files_multi/DE07NGSUKBR147099Aligned.sortedByCoord.out.bam

mkdir featureCounts_out

cat sample_list.txt | while read sample; do
	featureCounts -p -B -C -M --fraction -s 2 -T 8 -a /media/kilian/OS/References/GTF_files_TEtranscript/GRCm38_GENCODE_L1Md_unique.gtf \
 -o featureCounts_out/${sample}_featurecounts_output.txt BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam
done



#HUMAN
#Filter out only L1HS elements in repeatmasker file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip rmsk.txt.gz


awk '$0 ~ /gene_id "L1HS/' hg38_rmsk.gtf > hg38_rmsk_L1HS.gtf


#filter out transcript overlapping L1HS
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gunzip hg38.refGene.gtf.gz
#filtering out only instances of full transcripts HUMAN
/media/kilian/OS/BulkRNAseq-Analysis/Modifications_TEtranscript_analysis/filter_gtf.sh hg38.refGene.gtf > hg38.transcripts_filtered.gtf

bedtools subtract -A -a GRCh38_GENCODE_rmsk_L1HS.gtf -b hg38.transcripts_filtered.gtf > GRCh38_GENCODE_rmsk_L1HS_filtered.gtf


#Modify the gtf file to suit featureCounts counting!
#!/bin/bash
# Input and output file paths
INPUT_GTF="GRCh38_GENCODE_rmsk_L1HS.gtf"
OUTPUT_GTF="GRCh38_GENCODE_rmsk_L1HS_unique.gtf"

# Initialize a counter
counter=1

# Create or overwrite the output file
> "$OUTPUT_GTF"

# Read the input GTF line by line
while IFS= read -r line; do
    # Skip comment lines
    if [[ $line == \#* ]]; then
        echo "$line" >> "$OUTPUT_GTF"
    else
        # Append a unique counter to gene_id
        modified_line=$(echo "$line" | sed -E "s/(gene_id \"L1HS)(\";)/\1_$counter\2/")
        echo "$modified_line" >> "$OUTPUT_GTF"
        counter=$((counter + 1))
    fi
done < "$INPUT_GTF"

echo "Modified GTF file saved to $OUTPUT_GTF"




#FeatureCounts human with intronic L1
mkdir featureCounts_out
cat sample_list.txt | while read sample; do
    featureCounts -p -B -C -M --fraction -s 2 -T 8 -a /media/kilian/OS/References/GTF_files_featureCounts_hg38/GRCh38_GENCODE_rmsk_L1HS_unique.gtf \
 -o featureCounts_out/${sample}_featurecounts_output.txt BAM_files_multi/${sample}Aligned.sortedByCoord.out.bam
done





