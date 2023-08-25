#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import numpy as np
import os 


parser = argparse.ArgumentParser(description='prepareFastaDataBySpecies python script')
parser.add_argument('-r','--rawFasta_file', type=str, help='raw fasta file')
parser.add_argument('-s','--STRING_fastaBySpecies', type=str, help='new folder to save parsed fasta file by species')

args = parser.parse_args()
rawFasta_file=args.rawFasta_file
STRING_fastaBySpecies=args.STRING_fastaBySpecies




current_speFasta=""
current_speId=np.inf
species_Id = np.inf
with open(rawFasta_file) as input_handle:
    for line in input_handle:


        if  line.startswith(">"):
            species_Id=int(line[1:].split(".")[0])
            if current_speId !=np.inf and current_speId!=species_Id:
                with open(STRING_fastaBySpecies+str(current_speId) + ".fa", 'w') as output_handle:
                    output_handle.write(current_speFasta)
                current_speId=species_Id
                current_speFasta=""
                current_speFasta +=line

            else:
                current_speId=species_Id
                current_speFasta +=line
            
        else: 
            current_speFasta +=line
            
with open(os.path.join(STRING_fastaBySpecies,str(current_speId) + ".fa",), 'w') as output_handle:
    output_handle.write(current_speFasta)
    