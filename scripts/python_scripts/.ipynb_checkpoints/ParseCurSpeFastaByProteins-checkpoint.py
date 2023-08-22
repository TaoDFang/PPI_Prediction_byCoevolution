#!/usr/bin/env python
# coding: utf-8


import argparse

sys.path.append("../../src/utilities")

from create_singleMSA import ParseCurSpeFastaByProteins_STRNG1105


parser = argparse.ArgumentParser(description='ParseCurSpeFastaByProteins_STRNG1105')
parser.add_argument('-f','--currentSpe_fastaData', type=str, help='currentSpe_fastaData)
parser.add_argument('-b','--currentSpeProSeqPath_ByProteins', type=str, help='currentSpeProSeqPath_ByProteins')
parser.add_argument('-n','--currentSpe_protein_info_filename', type=str, help='currentSpe_protein_info_filename')


args = parser.parse_args()
currentSpe_fastaData=args.currentSpe_fastaData
currentSpeProSeqPath_ByProteins=args.currentSpeProSeqPath_ByProteins
currentSpe_protein_info_filename=args.currentSpe_protein_info_filename

ParseCurSpeFastaByProteins_STRNG1105(currentSpe_fastaData,currentSpeProSeqPath_ByProteins)
                    
print("currentSpeProSeqPath_ByProteins:",currentSpeProSeqPath_ByProteins)
pro_num=glob.glob(currentSpeProSeqPath_ByProteins+"*")
print("len(pro_num):",len(pro_num))

# ParseCurSpeFastaByProteins_STRNG1105 function sometimes dont correctly sperate all single proteins 
# here to makesure they do  , otherwise throw error 
currentSpe_protein_info_frame=pd.read_csv(currentSpe_protein_info_filename,skiprows=1,
                                          header=None,index_col=None,sep="\t")
print("currentSpe_protein_info_frame.shape",currentSpe_protein_info_frame.shape)


assert currentSpe_protein_info_frame.shape[0]== len(pro_num)