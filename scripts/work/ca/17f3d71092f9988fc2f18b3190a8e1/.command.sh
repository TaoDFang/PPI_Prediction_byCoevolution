#!/bin/bash -ue
mkdir -p ${Base_Folder}nextflow

# download all fasta seq data from new string website 
mkdir -p ${Base_Folder}nextflow/STRING_Data_11.5
wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz .
gzip -d protein.sequences.v11.5.fa.gz
