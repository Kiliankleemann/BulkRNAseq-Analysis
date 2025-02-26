#Squire installation

conda create -n squire
conda activate squire
#conda install r-base=3.4.1
conda install python=2.7
conda install gcc=10

# circumvent pip error
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o
get-pip.py
python get-pip.py --force-reinstall

git clone https://github.com/wyang17/SQuIRE; cd SQuIRE;
pip install -e .

squire Build -s all

conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig ucsc-genepredtobed
conda install -c bioconda star=2.5.3a

squire Fetch -b hg38 -f -c -r -g -x -p 8 -v
squire Fetch -b mm10 -f -c -r -g -x -p 8 -v
