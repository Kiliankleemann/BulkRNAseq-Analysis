#!/bin/bash
# Input and output file paths
INPUT_GTF="mm10_L1Md_filtered.gtf"
OUTPUT_GTF="mm10_L1Md_filtered_unique.gtf"

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
        modified_line=$(echo "$line" | sed -E "s/(gene_id \"[^\"]+)(\";)/\1_$counter\2/")
        echo "$modified_line" >> "$OUTPUT_GTF"
        counter=$((counter + 1))
    fi
done < "$INPUT_GTF"

echo "Modified GTF file saved to $OUTPUT_GTF"
