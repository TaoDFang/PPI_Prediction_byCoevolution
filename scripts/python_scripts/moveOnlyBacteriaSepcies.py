#!/usr/bin/env python
# coding: utf-8


import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='moveOnlyBacteriaSepcies python script')
parser.add_argument('-s','--STR_species_mem_file', type=str, help='STR_species_mem_file')
parser.add_argument('-s','--STRING_fastaBySpecies', type=str, help='folder that save parsed fasta file by species')
parser.add_argument('-s','--STRING_fastaByBacteriaSpecies', type=str, help='folder that save only fasta file for bacteria')


args = parser.parse_args()
STR_species_mem_file=args.STR_species_mem_file
STRING_fastaBySpecies=args.STRING_fastaBySpecies
STRING_fastaByBacteriaSpecies=args.STRING_fastaByBacteriaSpecies



STR_species_mem=pd.read_csv(STR_species_mem_file,sep="\t",index_col=None,header=0)
print(STR_species_mem.shape)

STR_bacteria=STR_species_mem.loc[STR_species_mem.iloc[:,4]=="Bacteria",:]
print(STR_bacteria.shape)
STR_backteria_id_list=list(STR_bacteria.iloc[:,0])
STR_backteria_id_list[1:10]


rePureName=re.compile(".*/(.*).fa")
for bid in STR_backteria_id_list:

    print(STRING_fastaBySpecies,STRING_fastaByBacteriaSpecies)
    #copyfile(STRING_fastaBySpecies+str(bid)+".fa", STRING_fastaByBacteriaSpecies+str(bid)+".fa")

        

