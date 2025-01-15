def bed_to_gtf(input_bed, output_gtf):
    with open(input_bed, 'r') as bed, open(output_gtf, 'w') as gtf:
        for line in bed:
            if line.startswith("track") or line.strip() == "":
                continue  # Skip track lines or empty lines
            fields = line.strip().split('\t')
            chrom, start, end, name, score, strand = fields[:6]
            attributes = f'gene_id "{name}"; transcript_id "{name}";'
            gtf_line = f"{chrom}\tBED\tgene\t{int(start)+1}\t{end}\t{score}\t{strand}\t.\t{attributes}\n"
            gtf.write(gtf_line)

# Specify the input BED file and output GTF file paths
input_bed = "input.bed"
output_gtf = "output.gtf"

# Perform the conversion
bed_to_gtf(input_bed, output_gtf)
print(f"Converted {input_bed} to {output_gtf}")
