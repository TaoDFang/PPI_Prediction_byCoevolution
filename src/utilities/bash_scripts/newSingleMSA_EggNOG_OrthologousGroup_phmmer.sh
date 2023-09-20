#!/bin/bash

# better run in ipykernel_py3 enviroment 
# for this script , has positional augumente $1, $2, $3 ...
# represent :  current spe name ,idx start from 1 to total number of input files , Path to fasta file of eggnog orthoglos of current spe,
#           :  path to fasta files of current spe which each file contain only one protein seq, output,


#echo "test"
#echo $1, $2, $3, $4, $5, $6,$6
phmmer_path=$1
SGE_TASK_ID=$2
currentSpe_OrthologousGroup_Fa_path=$3
currentSpeProSeqPath_ByProteins=$4
currentSpe_phmmer_outPath=$5


OIFS=$IFS; IFS=$'\n'; array=($(ls ${currentSpe_OrthologousGroup_Fa_path})); IFS=$OIFS; 
echo "${array[0]}" # start from 0 ; notice ls -ls is not exactly the same as ls. it has one more addition  line
echo "${#array[@]}"
i=$((SGE_TASK_ID-1))

name=${array[$i]%.fa}
#/mnt/mnemo5/tao/BeeiveProgram/HMMER3/bin/phmmer
"${phmmer_path}" -o "${currentSpe_phmmer_outPath}${name}" --tblout "${currentSpe_phmmer_outPath}${name}_tblout" --domtblout "${currentSpe_phmmer_outPath}${name}_domtblout"  --notextw  "${currentSpeProSeqPath_ByProteins}${array[$i]}"  "${currentSpe_OrthologousGroup_Fa_path}${array[$i]}" 

