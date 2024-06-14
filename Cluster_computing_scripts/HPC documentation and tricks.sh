HPC documentation and tricks


Login to cluster:
ssh kkleemann@bonna.hpc.uni-bonn.de

sshare
 ag_ukikckp_behrendt                   
  ag_ukikckp_behren+  kkleemann

#Partition otpions:
PARTITION    AVAIL  TIMELIMIT   JOB_SIZE ROOT OVERSUBS     GROUPS  NODES       STATE NODELIST
short           up    8:00:00       1-32   no       NO        all      8     planned node[12-13,17,23,28,30,33,36]
short           up    8:00:00       1-32   no       NO        all     39       mixed node[01-09,11,15-16,18-20,22,25-26,29,31,34-35,38-41,43,50-52,54-60,69-70]
short           up    8:00:00       1-32   no       NO        all     23   allocated node[10,14,21,24,27,32,37,42,44-49,53,61-68]
medium*         up 2-00:00:00       1-16   no       NO        all      8     planned node[12-13,17,23,28,30,33,36]
medium*         up 2-00:00:00       1-16   no       NO        all     39       mixed node[01-09,11,15-16,18-20,22,25-26,29,31,34-35,38-41,43,50-52,54-60,69-70]
medium*         up 2-00:00:00       1-16   no       NO        all     23   allocated node[10,14,21,24,27,32,37,42,44-49,53,61-68]
long            up 7-00:00:00        1-8   no       NO        all      8     planned node[12-13,17,23,28,30,33,36]
long            up 7-00:00:00        1-8   no       NO        all     39       mixed node[01-09,11,15-16,18-20,22,25-26,29,31,34-35,38-41,43,50-52,54-60,69-70]
long            up 7-00:00:00        1-8   no       NO        all     23   allocated node[10,14,21,24,27,32,37,42,44-49,53,61-68]
reservations    up   infinite       1-70   no       NO        all      8     planned node[12-13,17,23,28,30,33,36]
reservations    up   infinite       1-70   no       NO        all     39       mixed node[01-09,11,15-16,18-20,22,25-26,29,31,34-35,38-41,43,50-52,54-60,69-70]
reservations    up   infinite       1-70   no       NO        all     23   allocated node[10,14,21,24,27,32,37,42,44-49,53,61-68]



#Loading modules (packages)
module  avail                 # list available software (modules)
module  load [module]         # set up the environment to use the software
module  list                  # list currently loaded software
module  purge                 # clears the environment
module  help                  #get help

#Loading packages
module load Anaconda3
module load STAR/2.7.10b-GCC-11.3.0

#Directories aufsetzten
mkdir Line1_primer_STAR_alignment

#From local Terminal directory !!!
scp -r /Users/kiliankleemann/Dropbox/Anne_LINE1_primer_mouse_Aug_2023/fastq_files/* kkleemann@bonna.hpc.uni-bonn.de:~/Line1_primer_STAR_alignment/raw_fastq


#Download index from cluster
scp -r kkleemann@bonna.hpc.uni-bonn.de:~/Aref_DNA_damage_Oct_2023/BAM_files /Volumes/NAI3/P2023026RNALB/STAR_output

#Download BAM files from cluster
scp -r kkleemann@bonna.hpc.uni-bonn.de:~/Line1_primer_STAR_alignment/BAM_files /Volumes/OhneTitel/Anne_LINE1_primer_mouse_Aug_2023/BAM_files_STAR_multi
scp -r kkleemann@bonna.hpc.uni-bonn.de:~/Anne_BMDM_try_August_2023/BAM_files /Volumes/OhneTitel/Anne_BMDM_try_August_2023/BAM_files

scp -r kkleemann@bonna.hpc.uni-bonn.de:~/Anne_BMDM_try_August_2023/Aligned.sortedByCoord.out.bam /Volumes/OhneTitel/Anne_BMDM_try_August_2023/Aligned.sortedByCoord.out.bam



#copy between directory
cp Line1_primer_STAR_alignment/mm10_index/* Anne_BMDM_try_August_2023/mm10_index/


#HUMAN
Echo 'Download references #download reference from https://www.gencodegenes.org/human/ Genome sequence (GRCh38.p14)'
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d *.gz

#MOUSE
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.out.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gzip -d *.gz


#Copying files from one directory to another 
scp kkleemann@bonna.hpc.uni-bonn.de:~/Line1_primer_STAR_alignment/mm10_index/* kkleemann@bonna.hpc.uni-bonn.de:~/Anne_BMDM_try_August_2023/index_mm10/



