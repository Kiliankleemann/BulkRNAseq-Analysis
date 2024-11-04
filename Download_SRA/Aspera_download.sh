#install Aspera
mkdir ~/Desktop/Aspera
cd Aspera
wget https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.12.780_linux_x86_64.tar.gz
tar xvf ibm-aspera-connect_4.2.12.780_linux_x86_64.tar.gz 
bash ibm-aspera-connect_4.2.12.780_linux_x86_64.sh
#add env path of Aspera
echo 'export PATH=~/.aspera/connect/bin/:$PATH' >> ~/.bashrc # ~/.aspera/connect,  ~/.bashrc  (cd /home/xinwei/Desktop/.aspera/connect) 
source  ~/.bashrc  # ensure the env effective
which ascp #check the position
ascp -h #cheak if it works (work)
#make key txt with following content with named "asperaweb_id_dsa.openssh"
-----BEGIN DSA PRIVATE KEY-----
MIIBuwIBAAKBgQDkKQHD6m4yIxgjsey6Pny46acZXERsJHy54p/BqXIyYkVOAkEp
KgvT3qTTNmykWWw4ovOP1+Di1c/2FpYcllcTphkWcS8lA7j012mUEecXavXjPPG0
i3t5vtB8xLy33kQ3e9v9/Lwh0xcRfua0d5UfFwopBIAXvJAr3B6raps8+QIVALws
yeqsx3EolCaCVXJf+61ceJppAoGAPoPtEP4yzHG2XtcxCfXab4u9zE6wPz4ePJt0
UTn3fUvnQmJT7i0KVCRr3g2H2OZMWF12y0jUq8QBuZ2so3CHee7W1VmAdbN7Fxc+
cyV9nE6zURqAaPyt2bE+rgM1pP6LQUYxgD3xKdv1ZG+kDIDEf6U3onjcKbmA6ckx
T6GavoACgYEAobapDv5p2foH+cG5K07sIFD9r0RD7uKJnlqjYAXzFc8U76wXKgu6
WXup2ac0Co+RnZp7Hsa9G+E+iJ6poI9pOR08XTdPly4yDULNST4PwlfrbSFT9FVh
zkWfpOvAUc8fkQAhZqv/PE6VhFQ8w03Z8GpqXx7b3NvBR+EfIx368KoCFEyfl0vH
Ta7g6mGwIMXrdTQQ8fZs
-----END DSA PRIVATE KEY-----
mv '/home/xinwei/Desktop/Aspera/asperaweb_id_dsa.openssh' ~/.aspera/connect/etc
##find the position 
find ~ -name asperaweb_id_dsa.openssh #/home/xinwei/.aspera/connect/etc/asperaweb_id_dsa.openssh #good 
# download single file for two side reads (delete _1 if it is only one)
ascp  -v -Q -T -l 200m -P 33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR134/006/SRR1342456/   ./ 
#for SRR+8 number such as /SRR134/067/SRR13424567/: "SRR134/" from the first 3 num ; "067/" from 0 + the last two number 67
#for SRR+7 number such as /SRR134/006/SRR1342456/ : "SRR134/" from the first 3 num ; "006/" from 00 + the last number 6
#for SRR+6 number such as /SRR134/SRR134245/  : "SRR134/" from the first 3 num without the middle file

#download multiple task
## set the 
mkdir -p  ~/Aspera/raw/fq/ 
cd  ~/Aspera/raw/fq/
pwd

## 密匙路径
openssh=/home/xinwei/.aspera/connect/etc/asperaweb_id_dsa.openssh 

cat ~/Desktop/Aspera/SRR_Acc_List | while read id
do
num=`echo $id | wc -m `      #please attention，wc -m will treat "$" in the end as a character as well
echo "SRRnum+1 is $num"      #so SRRnum + 1 equal $num value
#if the number is SRR+8 #
if [ $num -eq 12 ]
then
        echo "SRR + 8"
        x=$(echo $id | cut -b 1-6)
        y=$(echo $id | cut -b 10-11)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001  -k 1 -i $openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/0$y/$id/   ./ & )
#if the number is SRR+7 #
elif [ $num -eq 11 ]
then
        echo  "SRR + 7"
        x=$(echo $id | cut -b 1-6)
        y=$(echo $id | cut -b 10-10)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001  -k 1 -i $openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/00$y/$id/   ./ & )
#if the number is SRR+6 #
elif [ $num -eq 10 ]
then
        echo  "SRR + 6"
        x=$(echo $id |cut -b 1-6)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001 -k 1 -i  $openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/$id/   ./ & )
fi
done
