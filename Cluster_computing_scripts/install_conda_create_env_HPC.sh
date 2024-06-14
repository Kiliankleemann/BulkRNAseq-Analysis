#!/bin/bash
#SBATCH --account=ag_ukikckp_behrendt
#SBATCH --job-name=install_conda_packages
#SBATCH --output=install_conda_packages.out
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5400M
#SBATCH --time=0:10:00

module load Anaconda3
conda create -n Bulk_RNA_TE_pipeline python=3.8
conda activate Bulk_RNA_TE_pipeline
conda install multiqc
conda install samtools
conda install salmon
conda install tetranscripts

