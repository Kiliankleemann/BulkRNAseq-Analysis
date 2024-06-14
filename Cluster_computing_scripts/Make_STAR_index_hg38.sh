#!/bin/bash
#SBATCH --account=ag_ukikckp_behrendt
#SBATCH --job-name=Index_STAR_genome
#SBATCH --output=Index_STAR_genome.out
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=2:50:00

module load STAR/2.7.10b-GCC-11.3.0
mkdir indexes
STAR --runMode genomeGenerate --genomeDir indexes \
            --genomeFastaFiles reference/hg38.p13.fa \
            --sjdbGTFfile reference/hg38.refGene.gtf \
            --sjdbOverhang 50 --outFileNamePrefix hg38