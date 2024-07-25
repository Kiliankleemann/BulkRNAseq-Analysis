list1=(
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/009/ERR6786539/ERR6786539_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/009/ERR6786539/ERR6786539_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/000/ERR6786540/ERR6786540_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/000/ERR6786540/ERR6786540_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/001/ERR6786541/ERR6786541_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/001/ERR6786541/ERR6786541_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/002/ERR6786542/ERR6786542_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/002/ERR6786542/ERR6786542_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/003/ERR6786543/ERR6786543_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/003/ERR6786543/ERR6786543_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/004/ERR6786544/ERR6786544_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/004/ERR6786544/ERR6786544_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/005/ERR6786545/ERR6786545_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/005/ERR6786545/ERR6786545_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/006/ERR6786546/ERR6786546_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/006/ERR6786546/ERR6786546_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/007/ERR6786547/ERR6786547_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/007/ERR6786547/ERR6786547_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/008/ERR6786548/ERR6786548_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/008/ERR6786548/ERR6786548_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/009/ERR6786549/ERR6786549_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/009/ERR6786549/ERR6786549_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/000/ERR6786550/ERR6786550_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/000/ERR6786550/ERR6786550_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/001/ERR6786551/ERR6786551_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/001/ERR6786551/ERR6786551_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/002/ERR6786552/ERR6786552_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/002/ERR6786552/ERR6786552_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/003/ERR6786553/ERR6786553_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/003/ERR6786553/ERR6786553_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/004/ERR6786554/ERR6786554_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/004/ERR6786554/ERR6786554_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/005/ERR6786555/ERR6786555_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/005/ERR6786555/ERR6786555_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/006/ERR6786556/ERR6786556_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/006/ERR6786556/ERR6786556_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/007/ERR6786557/ERR6786557_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/007/ERR6786557/ERR6786557_2.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/008/ERR6786558/ERR6786558_1.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR678/008/ERR6786558/ERR6786558_2.fastq.gz"
)
for i in ${list1[@]}; do
	echo ${i}
	wget ${i}
done

exit