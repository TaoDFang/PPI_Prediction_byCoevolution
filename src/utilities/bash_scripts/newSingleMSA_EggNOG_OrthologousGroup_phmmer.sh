#!/bin/bash

phmmer_path=$1
SGE_TASK_ID=$2
currentSpe_OrthologousGroup_Fa_path=$3
currentSpeProSeqPath_ByProteins=$4
currentSpe_phmmer_outPath=$5


OIFS=$IFS; IFS=$'\n'; array=($(ls ${currentSpe_OrthologousGroup_Fa_path})); IFS=$OIFS; 
echo "${array[0]}" # start from 0 ;
echo "${#array[@]}"
i=$((SGE_TASK_ID-1))

name=${array[$i]%.fa}
"${phmmer_path}" -o "${currentSpe_phmmer_outPath}${name}" --tblout "${currentSpe_phmmer_outPath}${name}_tblout" --domtblout "${currentSpe_phmmer_outPath}${name}_domtblout"  --notextw  "${currentSpeProSeqPath_ByProteins}${array[$i]}"  "${currentSpe_OrthologousGroup_Fa_path}${array[$i]}" 

