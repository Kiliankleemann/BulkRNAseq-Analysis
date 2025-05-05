#Human_F_M_HC_MCI_AD_E33_34
#Moving all raw files into final submission
for barcode in "101148 Sample_1_WT" "101151 Sample_2_WT" "101152 Sample_3_WT" "101154 Sample_4_WT" "101157 Sample_5_WT" "101149 Sample_6_KI" "101150 Sample_7_KI" "101155 Sample_8_KI" "101156 Sample_9_KI"

do
  set `echo $barcode`
  cat $1.fq.gz > $2.fq.gz
done
   

#Moving all processed files into final submission
for barcode in "101148 Sample_1_WT" "101151 Sample_2_WT" "101152 Sample_3_WT" "101154 Sample_4_WT" "101157 Sample_5_WT" "101149 Sample_6_KI" "101150 Sample_7_KI" "101155 Sample_8_KI" "101156 Sample_9_KI"

do
  set -- $barcode
  set `echo $barcode`
  cp -R $1/ quant_$2
done

#Zip all subfolders in selected merged files directory
cd selected_folders
for i in *
do
[ -d "$i" ] && zip -r "$i.zip" "$i"
done

rm -i !"*.zip"


