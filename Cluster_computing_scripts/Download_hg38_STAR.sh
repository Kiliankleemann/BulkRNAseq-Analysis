#!/bin/bash
#SBATCH --account=ag_ukikckp_behrendt
#SBATCH --job-name=Make_hg38_genome_STAR
#SBATCH --output=Make_hg38_genome_STAR.out
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3900M
#SBATCH --time=2:00

cd reference
Echo 'Download references #download reference from https://www.gencodegenes.org/human/ Genome sequence (GRCh38.p14)'
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d *.gz
