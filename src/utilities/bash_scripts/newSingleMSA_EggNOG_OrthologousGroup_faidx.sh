#! /bin/bash


currentSpe_TaxID=$1
SGE_TASK_ID=$2
currentSpe_fastaData=$3
origSTRINGBacteriaProSeqPath=$4
currentSpe_OrthologousGroup_Fa_path=$5
newsingleMSA_RBH_OrthologousGroup_fileName=$6


currentSpe_fastaData="$currentSpe_fastaData"

OIFS=$IFS
IFS=$','
ADDR=($(sed -n "${SGE_TASK_ID}p" ${newsingleMSA_RBH_OrthologousGroup_fileName}))
IFS=$OIFS


echo $currentSpe_fastaData
echo ${ADDR[0]}

# here replace >> with > at the begiing  i want to creat new files /overwrite old content when run again 
samtools faidx "$currentSpe_fastaData" "${ADDR[0]}" > "${currentSpe_OrthologousGroup_Fa_path}${ADDR[0]}.fa" 

ProSize=${#ADDR[@]}  
echo "$ProSize"
i=1
while [ $i -le $(($ProSize-1))  ]

do
  IFS='|' OthologArr=(${ADDR[$i]})
  samtools faidx "${origSTRINGBacteriaProSeqPath}${OthologArr[0]}.fa" "${OthologArr[1]}" >> "${currentSpe_OrthologousGroup_Fa_path}${ADDR[0]}.fa" 
  ((i++))
done


